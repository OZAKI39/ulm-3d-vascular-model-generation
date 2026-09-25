"""Verify this published snapshot, compressed evidence and frozen scientific inputs."""
from pathlib import Path
import argparse,hashlib,importlib.util,json
ROOT=Path(__file__).resolve().parents[2];META=Path(__file__).resolve().parent

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--scientific-inputs',action='store_true');a=p.parse_args();count=0
 for line in (META/'SNAPSHOT_SHA256.txt').read_text().splitlines():
  expected,relative=line.split('  ',1);path=(ROOT/relative).resolve()
  if not path.is_relative_to(ROOT):raise ValueError('Path outside snapshot')
  with path.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
  if actual!=expected:raise ValueError('SHA mismatch: '+relative)
  count+=1
 restored=load('p9a4_packed_verifier',META/'restore_large_artifacts.py').restore(ROOT,True)
 print(json.dumps(dict(published_file_hashes=count,compressed_originals_verified=len(restored),status='PASS'),indent=2))
 if a.scientific_inputs:load('inherited_scientific_verifier',ROOT/'sync_metadata/network_h0_particle_20260925/verify_snapshot.py').scientific_inputs()
if __name__=='__main__':main()
