"""Reversible deletion/restore checks; all labels are isolated software fixtures."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from annotation_store import CLASSES, save_annotation, load_annotation, sha
from annotation_library_store import saved_inventory, delete_saved, restore_saved


class LibraryStoreTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.record = {'id': 'fixture', 'resolution_wh': [12, 10], 'image_sha256': 'image-fixture',
                       'source_sha256': 'source-fixture', 'split': 'train'}
        self.masks = {k: np.zeros((10, 12), np.uint32) for k in CLASSES}
        self.masks['SINGLE_BUBBLE'][2:5, 2:5] = 9
        self.meta = save_annotation(self.root, self.record, self.masks, 'SOFTWARE_TEST', False,
                                    {'SINGLE_BUBBLE:9': 'track'}, 'software fixture',
                                    review_scope='REGIONS', workspace_roi=[1, 1, 8, 10])
        self.folder = self.root / 'fixture'

    def tearDown(self): self.temp.cleanup()

    def delete(self): return delete_saved(self.root, self.record, self.meta['revision'])

    def test_delete_restores_exact_masks_as_new_unreviewed_revision(self):
        history = {p.name: sha(p) for p in self.folder.glob('*') if p.name != 'latest.json'}
        deleted = self.delete()
        self.assertEqual(load_annotation(self.root, self.record), (None, None))
        self.assertEqual(saved_inventory(self.root, [self.record]), [])
        self.assertEqual(len(saved_inventory(self.root, [self.record], True)), 1)
        self.assertEqual(json.loads(deleted.read_text()), self.meta)
        restored = restore_saved(self.root, self.record, deleted)
        meta, masks = load_annotation(self.root, self.record)
        self.assertNotEqual(restored['revision'], self.meta['revision'])
        self.assertEqual(meta['status'], 'DRAFT_UNREVIEWED'); self.assertEqual(meta['review_regions'], [])
        for key in CLASSES: np.testing.assert_array_equal(masks[key], self.masks[key])
        for key in ('track_ids', 'notes', 'workspace_roi_yxyx'):
            self.assertEqual(meta[key], self.meta[key])
        self.assertEqual(saved_inventory(self.root, [self.record], True), [])
        self.assertTrue(deleted.with_suffix('.restored.json').exists())
        for name, digest in history.items(): self.assertEqual(sha(self.folder / name), digest)

    def test_stale_delete_rejected(self):
        with self.assertRaisesRegex(ValueError, '变化'): delete_saved(self.root, self.record, 'old')
        self.assertEqual(load_annotation(self.root, self.record)[0], self.meta)

    def test_failed_archive_preserves_active_pointer(self):
        original = sha(self.folder / 'latest.json')
        with patch.object(Path, 'replace', side_effect=OSError('fixture write failure')):
            with self.assertRaises(OSError): self.delete()
        self.assertEqual(sha(self.folder / 'latest.json'), original)

    def test_corrupt_masks_block_delete_and_restore(self):
        path = self.folder / self.meta['mask_file']; original = path.read_bytes(); path.write_bytes(b'broken')
        with self.assertRaisesRegex(ValueError, 'SHA-256'): self.delete()
        self.assertTrue((self.folder / 'latest.json').exists())
        path.write_bytes(original); deleted = self.delete(); path.write_bytes(b'broken')
        with self.assertRaisesRegex(ValueError, 'SHA-256'): restore_saved(self.root, self.record, deleted)
        self.assertFalse((self.folder / 'latest.json').exists()); self.assertTrue(deleted.exists())

    def test_restore_never_overwrites_new_annotation(self):
        deleted = self.delete()
        fresh = save_annotation(self.root, self.record, self.masks, notes='new work')
        with self.assertRaisesRegex(ValueError, '已有保存'): restore_saved(self.root, self.record, deleted)
        self.assertEqual(load_annotation(self.root, self.record)[0], fresh); self.assertTrue(deleted.exists())

    def test_inventory_keeps_empty_drafts_and_reports_invalid_metadata(self):
        for mask in self.masks.values(): mask.fill(0)
        save_annotation(self.root, self.record, self.masks)
        row = saved_inventory(self.root, [self.record])[0]
        self.assertEqual(row['meta']['instance_counts']['SINGLE_BUBBLE'], 0)
        (self.folder / 'latest.json').write_text('{bad json')
        self.assertTrue(saved_inventory(self.root, [self.record])[0]['error'])

    def test_restoration_drops_approval(self):
        # This fixture is never in the real annotation tree or model training.
        self.meta = save_annotation(self.root, self.record, self.masks, 'Isolated test reviewer', True)
        deleted = self.delete(); restored = restore_saved(self.root, self.record, deleted)
        self.assertFalse(restored['expert_confirmed']); self.assertEqual(restored['status'], 'DRAFT_UNREVIEWED')

    def test_deleted_labels_excluded_from_training_and_inflight_admission(self):
        from sam_auto_training import snapshot, snapshot_still_valid
        project = self.root / 'isolated_project'; (project / 'data').mkdir(parents=True)
        image = project / 'fixture.bin'; image.write_bytes(b'isolated image identity fixture')
        record = dict(self.record, image='fixture.bin', image_sha256=sha(image))
        (project / 'data/candidate_manifest.json').write_text(json.dumps({'rows': [record]}))
        directory = project / 'data/annotations'
        meta = save_annotation(directory, record, self.masks, 'Isolated test reviewer', True)
        snap = snapshot(project)
        self.assertEqual(len(snap['rows']), 1); self.assertTrue(snapshot_still_valid(snap, project))
        archived = delete_saved(directory, record, meta['revision'])
        self.assertEqual(snapshot(project)['rows'], []); self.assertFalse(snapshot_still_valid(snap, project))
        restore_saved(directory, record, archived)
        self.assertEqual(snapshot(project)['rows'], []); self.assertFalse(snapshot_still_valid(snap, project))


if __name__ == '__main__': unittest.main()
