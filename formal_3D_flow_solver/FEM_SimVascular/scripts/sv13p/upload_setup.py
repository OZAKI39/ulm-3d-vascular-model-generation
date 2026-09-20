from remote import ROOT,upload
for n in ['environment_remote.py','runner_remote.py','flow_parser.py']:upload(ROOT/'scripts/sv13p'/n,n)
for src,name in [('reports/sv1_3p/reference_freeze.json','reference_freeze.json'),('reports/sv1_3p/source_patch.json','source_patch.json'),('configs/sv1_3p/policy.json','policy.json')]:upload(ROOT/src,'configs/'+name)
upload(ROOT/'patches/sv1_3p/reuse_within_timestep.patch','patches/reuse_within_timestep.patch')
print('Stage P configuration and isolated patch uploaded.')
