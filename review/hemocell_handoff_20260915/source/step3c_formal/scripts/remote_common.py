"""Portable read-only contracts, hashing and result IO; no source-project access."""
import csv
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np

FROZEN_CONVERGENCE_SHA256='f3b02ba2b0607144194af09ec930db8af1161fc89b75ee1dcee313986989d7de'
PORTS=['Qin','Qout01','Qout02','Qout03']
OFFSETS=[2,4,6]
FIELD_DTYPE=np.dtype([('index','<u8'),('rho','<f8'),('u','<f8',(3,))])

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())

def write(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    encoded=json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n'
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name+'.',dir=path.parent)
    try:
        with os.fdopen(fd,'w') as f:f.write(encoded);f.flush();os.fsync(f.fileno())
        os.replace(tmp,path)
    finally:
        if Path(tmp).exists():Path(tmp).unlink()

def relative_path(root,value):
    p=Path(value)
    if p.is_absolute() or '..' in p.parts:raise ValueError('Bundle binding must be relative: '+value)
    result=root/p
    if not result.resolve().is_relative_to(root.resolve()):raise ValueError('Path escapes bundle')
    return result

def contract(root):
    p=Path(root)/'frozen_contracts/convergence_contract.json'
    if sha(p)!=FROZEN_CONVERGENCE_SHA256:raise ValueError('Frozen convergence contract changed')
    return read(p)

def verify_bundle(root):
    root=Path(root);manifest=read(root/'manifest.json');failures=[]
    entries={e['path']:e for e in manifest['files']}
    for name,e in entries.items():
        p=relative_path(root,name)
        if not p.is_file() or p.is_symlink() or p.stat().st_size!=e['bytes'] or sha(p)!=e['sha256']:failures.append(name)
    sums={}
    for line in (root/'SHA256SUMS').read_text().splitlines():
        digest,name=line.split('  ',1);sums[name]=digest
    if set(sums)!=(set(entries)|{'manifest.json'}):failures.append('SHA256SUMS inventory')
    if sums.get('manifest.json')!=sha(root/'manifest.json'):failures.append('manifest.json')
    for name,e in entries.items():
        if sums.get(name)!=e['sha256']:failures.append('SHA256SUMS:'+name)
    c=contract(root)
    linked=[('frozen_contracts/multiplane_measurement_contract.json',c['flux']['geometry_sha256']),
            ('frozen_contracts/control_volume_geometry.json',c['control_volume']['geometry_sha256']),
            ('frozen_contracts/control_volume_nodes.npz',c['control_volume']['mask_sha256'])]
    linked += [('frozen_contracts/'+Path(e['copy']).name,e['sha256']) for e in c['physical_contract_copies']]
    measure=read(root/'frozen_contracts/multiplane_measurement_contract.json')
    linked.append(('frozen_contracts/multiplane_quadrature.tsv',measure['quadrature_sha256']))
    geometry=read(root/'frozen_contracts/geometry_reuse_contract.json');binding=read(root/'runtime_bindings.json')
    linked.append((binding['stl'],geometry['input_stl_sha256']))
    proof=read(root/'source/physics_reuse_proof.json')
    linked += [('source/vascularPureFluidLong.cpp',proof['long_source_sha256']),
               ('frozen_inputs/step3_baseline_source/vascularPureFluid.cpp',proof['baseline_sha256'])]
    for name,digest in linked:
        if sha(relative_path(root,name))!=digest:failures.append('frozen reference:'+name)
    if failures:raise ValueError('Bundle integrity failed: '+repr(failures))
    return {'status':'PASS','file_count':len(entries),'manifest_sha256':sha(root/'manifest.json'),
            'convergence_contract_sha256':FROZEN_CONVERGENCE_SHA256}

def verify_run_inputs(root,run):
    receipt=read(Path(run)/'provenance/input_receipt.json')
    if receipt['manifest_sha256']!=sha(Path(root)/'manifest.json'):raise ValueError('Run belongs to different bundle')
    for rel,digest in receipt['materialized_files'].items():
        if sha(relative_path(Path(run),rel))!=digest:raise ValueError('Materialized input changed: '+rel)
    if sha(Path(run)/'contracts/convergence_contract.json')!=FROZEN_CONVERGENCE_SHA256:raise ValueError('Run convergence contract changed')
    return receipt

def fields(path):
    path=Path(path)
    with path.open('rb') as f:
        count=np.fromfile(f,dtype='<u8',count=1)
        if count.size!=1:raise ValueError('Missing field count')
        data=np.fromfile(f,dtype=FIELD_DTYPE)
    if len(data)!=count[0] or path.stat().st_size!=8+40*len(data):raise ValueError('Truncated field '+str(path))
    if not np.isfinite(data['rho']).all() or not np.isfinite(data['u']).all():raise ValueError('Nonfinite field')
    return data

def history(run):
    data=np.genfromtxt(Path(run)/'diagnostics/flow_history.csv',names=True,delimiter=',',dtype=float)
    return np.atleast_1d(data)

def csv_write(path,rows,columns):
    path=Path(path)
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)

def field_path(run,step,samples=False):
    return Path(run)/f'diagnostics/field_samples/{"samples" if samples else "fields"}_{step}.bin'

# Formal harness bindings only. Frozen convergence algebra above is unchanged.
def verify_run_inputs(root,run):
    root=Path(root);run=Path(run)
    receipt=read(root/'provenance/INPUT_RECEIPT.json')
    for rel,digest in receipt['materialized_files'].items():
        if sha(run/rel)!=digest:raise ValueError('Frozen formal input changed: '+rel)
    if sha(run/'contracts/convergence_contract.json')!=FROZEN_CONVERGENCE_SHA256:raise ValueError('Convergence contract changed')
    return receipt

def verify_bundle(root):
    root=Path(root);freeze=read(root/'FORMAL_STEP3C_GPU_FREEZE.json')
    for path,digest in freeze['execution_file_hashes'].items():
        if sha(Path(path))!=digest:raise ValueError('Frozen execution file changed: '+path)
    return {'status':'PASS','frozen_files':len(freeze['execution_file_hashes'])}
