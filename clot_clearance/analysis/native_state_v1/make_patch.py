"""Apply minimal reference-state storage repair ONLY to isolated source copy."""
from pathlib import Path
import shutil,json,difflib,hashlib,subprocess
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run'];src=R/'tools/lammps-contact-onset-v1';dst=O/'source'
# Required evidence exists BEFORE C++ mutation.
assert json.loads((O/'original_sort/report.json').read_text())['T1_sort']['wrong_m_node_count']>0
assert json.loads((O/'original_migration/report.json').read_text())['T2_migration']['wrong_m_node_count']>0
dst.mkdir(exist_ok=True)
subprocess.run(['rsync','-a','--exclude=.git','--exclude=build','--exclude=__pycache__',str(src)+'/',str(dst)+'/'],check=True)
p=dst/'src/BPM/bond_bpm_peri.cpp';s=p.read_text();s=s.replace('index_vinter(-1), index_dtheta(-1), index_ddamage(-1),','index_vinter(-1), index_dtheta(-1), index_ddamage(-1),\n    id_fix_reference_m(nullptr), index_reference_m(-1), index_reference_ready(-1),')
s=s.replace('  memory->destroy(wvolume);\n','')
s=s.replace('      memory->create(wvolume, nmax, "bond:wvolume");\n','')
s=s.replace('  if (id_fix_property_peri) {','  if (id_fix_reference_m) {\n    modify->delete_fix(id_fix_reference_m);\n    delete[] id_fix_reference_m;\n  }\n\n  if (id_fix_property_peri) {',1)
s=s.replace('  // grow the per-step break scratch', '''  // Reference m belongs to atom identity, not a local-index scratch slot.
  // Refresh this alias after sorting/exchange/grow; property/atom owns storage.
  if (state_based) wvolume = atom->dvector[index_reference_m];

  // grow the per-step break scratch''',1)
s=s.replace('      wvolume_setup = 0;    // arrays reallocated -> rebuild the weighted volume','      wvolume_setup = 0;    // refresh ghosts; initialized reference m is retained')
s=s.replace('  for (int i = 0; i < ntotal; i++) wvolume[i] = 0.0;', '''  double *ready = atom->dvector[index_reference_ready];
  // Only virgin atoms receive a reference sum. A zero-neighbor atom legitimately
  // has m=0, hence a separate restart-persistent initialized flag is necessary.
  // Never reconstruct initialized m from the surviving post-fracture topology.
  for (int i = 0; i < nlocal; i++)
    if (ready[i] == 0.0) wvolume[i] = 0.0;''')
s=s.replace('if (i1 < nlocal) wvolume[i1] +=','if (i1 < nlocal && ready[i1] == 0.0) wvolume[i1] +=').replace('if (i2 < nlocal) wvolume[i2] +=','if (i2 < nlocal && ready[i2] == 0.0) wvolume[i2] +=')
s=s.replace('  wvolume_setup = 1;','  for (int i = 0; i < nlocal; i++) ready[i] = 1.0;\n  wvolume_setup = 1;')
marker='  // optional dilatation diagnostic:'
insert='''  // A dedicated fix gives reference m formal copy/grow/exchange/border/restart
  // callbacks. Do not piggyback on s0/smin/vfrac or reset restored values here.
  // On add_fix(), Modify restores matching per-atom restart data by fix ID.
  if (state_based && !id_fix_reference_m) {
    if (atom->find_custom("bpm_peri_reference_m", flag, cols) >= 0 ||
        atom->find_custom("bpm_peri_reference_ready", flag, cols) >= 0 ||
        modify->get_fix_by_id("BPM_PERI_REFERENCE_M"))
      error->all(FLERR, "Reserved BPM_PERI_REFERENCE_M storage already exists");
    id_fix_reference_m = utils::strdup("BPM_PERI_REFERENCE_M");
    modify->add_fix(fmt::format("{} all property/atom d_bpm_peri_reference_m "
                                "d_bpm_peri_reference_ready ghost yes writedata no",
                                id_fix_reference_m));
    index_reference_m = atom->find_custom("bpm_peri_reference_m", flag, cols);
    index_reference_ready = atom->find_custom("bpm_peri_reference_ready", flag, cols);
  }

'''
assert marker in s;s=s.replace(marker,insert+marker);p.write_text(s)
p=dst/'src/BPM/bond_bpm_peri.h';s=p.read_text().replace('  // per-step scratch for the break bookkeeping', '''  // Owned by a dedicated property/atom fix, including initialized-state marker.
  char *id_fix_reference_m;
  int index_reference_m, index_reference_ready;

  // per-step scratch for the break bookkeeping''');p.write_text(s)
diff='';ids={}
for rel in ['src/BPM/bond_bpm_peri.cpp','src/BPM/bond_bpm_peri.h']:
 orig=(src/rel).read_text();new=(dst/rel).read_text();diff+=''.join(difflib.unified_diff(orig.splitlines(True),new.splitlines(True),'a/'+rel,'b/'+rel));ids[rel]={'original':hashlib.sha256(orig.encode()).hexdigest(),'fixed':hashlib.sha256(new.encode()).hexdigest()}
(O/'BUGFIX_SOURCE_DIFF.patch').write_text(diff);(O/'PATCH_IDENTITY.json').write_text(json.dumps(ids,indent=2));print(ids)
