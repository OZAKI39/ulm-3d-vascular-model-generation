from remote import ROOT,upload
for n in ['environment_remote.py','runner_remote.py','flow_parser.py']:upload(ROOT/'scripts/sv13q'/n,n)
for src,name in [('reports/sv1_3q/reference_freeze.json','reference_freeze.json'),('reports/sv1_3q/source_patch.json','source_patch.json'),('configs/sv1_3q/policy.json','policy.json')]:upload(ROOT/src,'configs/'+name)
upload(ROOT/'patches/sv1_3q/ilu_rebuild_policy.patch','patches/ilu_rebuild_policy.patch')
print('Stage Q configuration and isolated patch uploaded.')
