from sv13p_support import *
def test_local_and_native_history_untouched():
 d=accepted('preservation_audit');assert not d['historical_changes'] and d['official_source_unchanged']
 d=accepted('remote/remote_preservation');assert all(d['checks'].values())
 for f in json.loads((ROOT/'reports/sv1_3o/delivery_manifest.json').read_text())['files']:assert hashlib.sha256((ROOT/f['path']).read_bytes()).hexdigest()==f['sha256']
