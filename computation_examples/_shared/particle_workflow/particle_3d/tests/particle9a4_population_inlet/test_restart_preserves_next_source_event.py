import pytest,json
from particle_3d.population_inlet_p9a4 import generate
from particle_3d.continuous_infusion import PopulationLedger,canonical_bytes

def test_rejected_events_restore_without_redraw(simple_source,tmp_path):
 whole=generate(simple_source,128);first=generate(simple_source,53)
 assert any(r['particle_id'] is None for r in first.rows)
 first.checkpoint(tmp_path/'cp');restored=PopulationLedger.restore(tmp_path/'cp',simple_source.identity)
 assert restored.next_source_event_id==54
 generate(simple_source,75,ledger=restored,workers=3)
 assert canonical_bytes(restored.rows)==canonical_bytes(whole.rows)

def test_overwrite_tamper_and_partial_refused(simple_source,tmp_path):
 ledger=generate(simple_source,15);folder=tmp_path/'cp';ledger.checkpoint(folder)
 with pytest.raises(FileExistsError):ledger.checkpoint(folder)
 p=folder/'proposal_ledger.jsonl';p.write_text(p.read_text().splitlines()[0]+'\n')
 with pytest.raises(ValueError):PopulationLedger.restore(folder,simple_source.identity)
 partial=tmp_path/'partial';partial.mkdir()
 with pytest.raises(FileNotFoundError):PopulationLedger.restore(partial,simple_source.identity)
