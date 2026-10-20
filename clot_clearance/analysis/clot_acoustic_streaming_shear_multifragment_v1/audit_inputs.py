"""Read actual upstream files; no bubble trajectory is promoted to fluid velocity."""
from common import *
import h5py,shutil,time,platform,importlib.metadata

def environment():
 packages={}
 for p in ['numpy','scipy','matplotlib','h5py','mpi4py','optuna','botorch']:
  try:packages[p]=importlib.metadata.version(p)
  except importlib.metadata.PackageNotFoundError:packages[p]='NOT_INSTALLED_IN_NATIVE_ENV'
 dump(O/'ENVIRONMENT.json',dict(python=platform.python_version(),packages=packages,library_sha256=sha(LIB),optimizer_note='Optional optimizer uses its previous isolated environment; no installation or trial needed this round.'))

def main():
 start=time.perf_counter();O.mkdir(parents=True,exist_ok=True);F.mkdir(parents=True,exist_ok=True)
 assert not (O/'BASELINE_MANIFEST.json').exists()
 old=json.loads((B/'DELIVERY_MANIFEST.json').read_text())
 assert sha(B/'DELIVERY_MANIFEST.json')=='7471e7d5e0a9b3e8f147de9a625185a5a3f8793ffb551283c01406ef15b39ef7'
 (O/'provenance').mkdir()
 for name in ['reference.npz','reference.native','body.data']:shutil.copy2(B/'provenance'/name,O/'provenance'/name)
 baseline={str(p):dict(sha256=sha(p),bytes=p.stat().st_size) for p in [LIB,B/'DELIVERY_MANIFEST.json',R/'analysis/v3_softening_pilot_v1/native.py',R/'analysis/softening_feasibility_v1/model.py',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py',*(O/'provenance').iterdir()]}
 assert sha(LIB)=='465f8144c6da5c9de3aec4137b9c325f53c26737b199c536fb4a20845511ef77'
 dump(O/'BASELINE_MANIFEST.json',dict(native_unchanged=True,material_unchanged=True,files=baseline,scope='Only new fluid/loading/observer adapters; no indenters, no new precuts, no constitutive/MPI/restart changes'))
 oldaudit=json.loads((INV/'MEVO_ASSET_AUDIT.json').read_text());assets={};schemas={}
 for name,rec in oldaudit['assets'].items():
  p=Path(name);assets[name]=dict(exists=p.exists())
  if p.exists():
   assets[name].update(sha256=sha(p),bytes=p.stat().st_size,identical_to_previous_audit=sha(p)==rec['sha256'])
   if p.suffix=='.h5':
    with h5py.File(p) as h:
     ds={};h.visititems(lambda key,obj:ds.update({key:dict(shape=obj.shape,dtype=str(obj.dtype))}) if isinstance(obj,h5py.Dataset) else None)
     schemas[name]=dict(datasets=ds,attrs={k:str(v)[:20000] for k,v in h.attrs.items()})
 # Actual same files re-read, with current hashes, plus compact values already audited.
 shutil.copy2(INV/'MEVO_ACTUAL_EVIDENCE.npz',O/'MEVO_ACTUAL_EVIDENCE.npz')
 mevo=Path('/workspace/mevo_acoustic_microbubbles');exe=mevo/'f10s_control_validation/.execution/14731fd56f0477c89cbbdea5ae1f13e8e38e3f3c'
 scope=[]
 for root in [exe/'results/straight_tube_exp1/f10n',mevo/'references']:
  for p in root.rglob('*'):
   if p.is_file() and (p.suffix in ['.h5','.npz','.pdf'] or p.name.startswith('RP_')):scope.append(dict(path=str(p),bytes=p.stat().st_size))
 # Read coupling-relevant source lines, preserving exact paths/line numbers.
 excerpts=[]
 for p in [exe/'src/mevo_mb/fast_sim/water_b2.py',exe/'src/mevo_mb/bubbles/apecss.py',exe/'src/mevo_mb/bubbles/adaptive_table.py']:
  lines=p.read_text().splitlines();hits=[i for i,s in enumerate(lines) if any(w in s.lower() for w in ['strict_zero','flow_enabled','harmonic','radius','velocity','r(t)','results'])]
  indices=sorted(set(j for i in hits for j in range(max(0,i-1),min(len(lines),i+2))))
  excerpts.append(dict(path=str(p),sha256=sha(p),lines=[dict(line=i+1,text=lines[i]) for i in indices]))
 audit=dict(status='ACTUAL_FILES_REOPENED',assets=assets,h5_schemas=schemas,current_acoustics=oldaudit['current_acoustics'],historical_kwave=oldaudit['historical_kwave'],apecss=oldaudit['apecss'],trajectory=oldaudit['trajectory'],searched_file_inventory=scope,source_excerpts=excerpts,
  available=['F1 complex pressure/gradient on its own geometry','F10N/F10S pressure basis proxy and phase convention','APECSS volume-harmonic response table over valid adaptive leaves','F10S particle positions/times and sampled forces'],
  unavailable_for_clot=['registered bubble locations in current small-specimen/V3 coordinates','measured bubble-to-clot gap','near-wall no-slip first-order oscillatory velocity/gradient','bubble surface kinematics including translation/radial phase at clot','cycle-averaged Reynolds stress with consistent boundary terms','solved second-order mean flow and pressure on clot geometry'],
  radial_caveat='A harmonic response cache is not a time-resolved near-wall fluid field. Underlying radial solver can generate R(t), but even full R(t) would not fix unknown geometry/boundary conditions.',
  decision='STREAMING_LOAD_NOT_YET_DETERMINED',selected_fallback='PRESCRIBED_STOKES_SHEAR_TEMPLATE',bubble_to_clot_distance_m='UNAVAILABLE',bubble_streaming_stress_Pa='UNAVAILABLE',no_other_project_modified=True)
 dump(O/'ACOUSTIC_AND_BUBBLE_INPUT_AUDIT.json',audit)
 pipeline=[('MEVO_acoustics','IMPLEMENTED','F1 actual acoustics and F10S geometric proxy, separate geometries'),('APECSS_radial','REUSABLE','actual response cache, not near-wall fluid'),('bubble_to_clot_registration','UNAVAILABLE','no common coordinate/boundary definition'),('first_order_near_wall_fluid','UNAVAILABLE','pressure gradients alone omit boundary layers'),('second_order_streaming','INCOMPLETE','no solved source/boundary-consistent flow'),('prescribed_Stokes','REUSABLE','upgrade to one continuous transverse counter-shear mode'),('V3_exposed_surface','REUSABLE',str(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz')),('surface_node_mapping','REUSABLE','audit signed affine map and local nonnegative comparison'),('native_external_force','REUSABLE','existing fix external callback, physical N'),('native_fracture','REUSABLE',str(LIB)),('attachment_topology','REUSABLE','original live bulk/interface graph'),('Optuna_BoTorch','REUSABLE','installed interface, zero new trials'),('rigid_indenters','IMPLEMENTED','historical evidence only; excluded from current loading')]
 dump(O/'ACTIVE_STREAMING_SHEAR_PIPELINE.json',dict(components=[dict(name=n,status=s,scope=d) for n,s,d in pipeline],deleted_files=[],reason_no_deletion='Historical/provenance modules remain referenced. New entrypoint excludes indenters without physical deletion.'))
 dump(O/'INPUT_AUDIT_RESOURCES.json',dict(wall_s=time.perf_counter()-start))
 environment()
 print('INPUT_AUDIT_DONE',len(assets),len(scope))
if __name__=='__main__':main()
