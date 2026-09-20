import difflib,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import sha256
base=ROOT/'external/sv13p/svMultiPhysics-reuse';target=ROOT/'external/sv13q/svMultiPhysics-reuse';name='Code/Source/solver/main.cpp';old=(base/name).read_text()
new=old.replace('#include <fstream>','#include <fstream>\n#include <chrono>')
needle='    const bool reached_stop_time_step = cTS >= stopTS;'
new=new.replace(needle,needle+'''
    // Stage Q: timestamp the existing native acknowledgement; stop logic is unchanged.
    if (reached_stop_time_step && stopTS < nTS && cm.mas(cm_mod)) {
      const double epoch = std::chrono::duration<double>(
          std::chrono::system_clock::now().time_since_epoch()).count();
      std::ofstream stamp("sv13q_stop_ack.json");
      stamp << std::setprecision(17) << "{\\\"solver_ack_timestamp\\\":" << epoch
            << ",\\\"actual_stop_step\\\":" << cTS << "}\\n";
    }
''')
assert new!=old;(target/name).write_text(new)
r=ROOT/'reports/sv1_3q/source_patch.json';d=json.loads(r.read_text());d['after'][name]=sha256(target/name);d['changed_files'].append(name);d['stop_instrumentation']='One native timestamp written on existing STOP_SIM acknowledgement; no change to stop condition, cadence, or steady criteria.'
patch=''
for n in d['changed_files']:
 patch+=''.join(difflib.unified_diff((base/n).read_text().splitlines(True),(target/n).read_text().splitlines(True),fromfile='a/'+n,tofile='b/'+n))
for n in d['added_files']:
 patch+=''.join(difflib.unified_diff([],(target/n).read_text().splitlines(True),fromfile='/dev/null',tofile='b/'+n))
p=ROOT/d['patch'];p.write_text(patch);d['patch_sha256']=sha256(p);r.write_text(json.dumps(d,indent=2)+'\n')
# Incremental remote patch applies only to the initial Q build, before any CFD.
p=ROOT/'patches/sv1_3q/stop_ack_timestamp.patch';p.write_text(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile='a/'+name,tofile='b/'+name)))
print('Native STOP acknowledgement timestamp added without changing stop semantics.')
