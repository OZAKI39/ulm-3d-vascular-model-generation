"""Second evidence-based guard: do not access owned-only bond arrays for ghosts."""
from pathlib import Path
import json,hashlib,difflib,subprocess
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path.cwd();O=R/'verification'/C['task']/C['run']
log=(O/'mpi_fracture_original_launcher.log').read_text();assert 'Segmentation fault' in log and 'process_broken' in log
src=O/'source';dst=O/'source_owner_fix';dst.mkdir(exist_ok=True);subprocess.run(['rsync','-a',str(src)+'/',str(dst)+'/'],check=True)
p=dst/'src/BPM/bond_bpm.cpp';s=p.read_text();i=s.index('void BondBPM::process_broken(');j=s.index('void BondBPM::process_new(',i);part=s[i:j]
assert part.count('if (id != -1)')==1 and part.count('if (jd != -1)')==1
part=part.replace('  if (id != -1) {','  // Ghosts have no valid owned bond/history arrays. With newton bond off,\n  // the remote owner removes its endpoint; a periodic mapped LOCAL copy is OK.\n  if (id >= 0 && id < nlocal) {').replace('  if (jd != -1) {','  if (jd >= 0 && jd < nlocal) {');s=s[:i]+part+s[j:];p.write_text(s)
files=['src/BPM/bond_bpm_peri.cpp','src/BPM/bond_bpm_peri.h','src/BPM/bond_bpm.cpp'];diff='';identity={}
for f in files:
 a=(R/'tools/lammps-contact-onset-v1'/f).read_bytes();b=(dst/f).read_bytes();identity[f]={'original_sha256':hashlib.sha256(a).hexdigest(),'final_sha256':hashlib.sha256(b).hexdigest()};diff+=''.join(difflib.unified_diff(a.decode().splitlines(True),b.decode().splitlines(True),'a/'+f,'b/'+f))
(O/'BUGFIX_WVOLUME_ONLY.patch').write_bytes((O/'BUGFIX_SOURCE_DIFF.patch').read_bytes());(O/'BUGFIX_SOURCE_DIFF.patch').write_text(diff)
identity['failed_V3_wvolume_only_library_sha256']=hashlib.sha256((O/'build/liblammps.so').read_bytes()).hexdigest();(O/'FINAL_PATCH_IDENTITY.json').write_text(json.dumps(identity,indent=2))
(O/'MPI_OWNED_BOND_GUARD_EVIDENCE.json').write_text(json.dumps({'status':'CONFIRMED_GHOST_OWNERSHIP_DEFECT','production_original_minimal_test':'original_mpi_fracture','original_exit_code':139,'V3_exit':'rank1 SIGSEGV process_broken+0x267','source_lines_original':[460,477],'defect':'id/jd ghost mapping may remain >=nlocal, but old condition only checked !=-1 before reading num_bond/bond_atom/bond_type and bond history. Ghost bond/history arrays are not valid. Assembly at +0x267 dereferences setflag[bond_type[jd][m]].','repair':'Guard both endpoint edits by 0<=index<nlocal; preserve event ownership count, special-bond scheduling, bond history shifts and criterion.','V3_retry_this_round':False,'scope_amendment':'One targeted two-rank fracture on existing341-node geometry justified by actual V3 crash; original production also crashes. No parameter sweep.'},indent=2))
