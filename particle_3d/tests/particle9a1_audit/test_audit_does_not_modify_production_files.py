from particle_3d.routing_stationary_audit import read,sha

def test_every_protected_original_file_has_identical_bytes(audit,root):
    for path,digest in read(audit/'data/protected_before.json').items():assert sha(root/path)==digest,path
