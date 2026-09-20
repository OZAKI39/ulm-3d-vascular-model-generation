"""Upload only WSL-authored Stage O instructions to the native build host."""
from remote import ROOT,upload
for n in ['environment_remote.py','runner_remote.py','flow_parser.py']:
 upload(ROOT/'scripts/sv13o'/n,n)
for local,dst in [('reports/sv1_3o/reference_freeze.json','reference_freeze.json'),('reports/sv1_3o/source_patch.json','source_patch.json'),('configs/sv1_3o/policy.json','policy.json')]:upload(ROOT/local,'configs/'+dst)
upload(ROOT/'patches/sv1_3o/final_output_on_stop.patch','patches/final_output_on_stop.patch')
print('Stage O setup inputs uploaded.')
