#!/usr/bin/env python3
"""Resolve official registry metadata only; this is not a successful image pull."""
import hashlib
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.provenance import now,write_json
repo='simvascular/solver';base='https://registry-1.docker.io/v2/'+repo
scope=urllib.parse.urlencode({'service':'registry.docker.io','scope':'repository:'+repo+':pull'})
token=json.load(urllib.request.urlopen('https://auth.docker.io/token?'+scope,timeout=30))['token']
accept=', '.join(['application/vnd.oci.image.index.v1+json','application/vnd.docker.distribution.manifest.list.v2+json','application/vnd.oci.image.manifest.v1+json','application/vnd.docker.distribution.manifest.v2+json'])
def request(suffix):
    request=urllib.request.Request(base+suffix,headers={'Authorization':'Bearer '+token,'Accept':accept})
    response=urllib.request.urlopen(request,timeout=45);raw=response.read()
    return json.loads(raw),dict(response.headers),raw
manifest,headers,raw=request('/manifests/latest');tag_digest=headers.get('Docker-Content-Digest') or 'sha256:'+hashlib.sha256(raw).hexdigest()
assert tag_digest=='sha256:'+hashlib.sha256(raw).hexdigest()
index_digest=None
if 'manifests' in manifest:
    index_digest=tag_digest;matches=[m for m in manifest['manifests'] if m['platform']['os']=='linux' and m['platform']['architecture']=='amd64']
    assert len(matches)==1
    tag_digest=matches[0]['digest'];manifest,headers,raw=request('/manifests/'+tag_digest)
    assert tag_digest=='sha256:'+hashlib.sha256(raw).hexdigest()
config_digest=manifest['config']['digest'];config,headers,raw_config=request('/blobs/'+config_digest)
assert config_digest=='sha256:'+hashlib.sha256(raw_config).hexdigest()
output=ROOT/'outputs/sv0/environment';output.mkdir(parents=True,exist_ok=True)
(output/'solver_registry_manifest.json').write_bytes(raw);(output/'solver_registry_config.json').write_bytes(raw_config)
write_json(ROOT/'reports/sv0/solver_registry_resolution.json',{'timestamp':now(),'repository':repo,'requested_tag':'latest',
    'resolved_immutable_reference':repo+'@'+tag_digest,'digest':tag_digest,'index_digest':index_digest,
    'config_digest':config_digest,'created':config.get('created'),'os':config.get('os'),'architecture':config.get('architecture'),
    'labels':config.get('config',{}).get('Labels'),'layer_count':len(manifest['layers']),
    'compressed_layer_bytes':sum(l['size'] for l in manifest['layers']),
    'image_layers_pulled':False,'container_run':False,'evidence_scope':'Official public registry manifest and image configuration only; not docker inspect of a locally installed image'})
print(json.dumps({'immutable_reference':repo+'@'+tag_digest,'created':config.get('created'),'compressed_layer_bytes':sum(l['size'] for l in manifest['layers'])},indent=2))
