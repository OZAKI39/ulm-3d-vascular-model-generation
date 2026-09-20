from sv13n_support import *
def test_official_source():
 d=accepted('petsc325_source');source_gate(d)
 assert hashlib.sha256((ROOT/d['source_archive']).read_bytes()).hexdigest()==d['source_archive_sha256']
 assert d['commit'] in (R/'official_tag_lookup.txt').read_text()
 accepted('archive_tag_crosscheck')
def test_floating_release_rejected():
 d=read('petsc325_source');d.update(candidate='Candidate R',branch='release',commit='HEAD')
 with pytest.raises(GateError):source_gate(d)
