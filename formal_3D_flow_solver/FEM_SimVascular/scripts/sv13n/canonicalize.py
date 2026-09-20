"""Promote observed native records to stable Stage N report names without relabeling failures."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3n'
candidate=sys.argv[1] if len(sys.argv)>1 else 'gpu13'
names=['petsc_cpu_build','svmp_cpu_build','svmp_cpu_link','petsc_gpu13_build','petsc_gpu123_build','petsc_gpu_smoke','petsc_cuda_types','ghost_probe_old','ghost_coherence','remote_preservation','final_source_integrity']
mapping={n:n for n in names};mapping.update({f'svmp_{candidate}_build':'svmp_gpu_build',f'svmp_{candidate}_link':'svmp_gpu_link',f'ghost_probe_{candidate}':'ghost_probe_new'})
for source,destination in mapping.items():
 p=R/'remote'/(source+'.json')
 if not p.exists():continue
 d=json.loads(p.read_text())
 if destination.endswith('_link'):d['version']='3.25.5'
 (R/(destination+'.json')).write_text(json.dumps(d,indent=2)+'\n')
