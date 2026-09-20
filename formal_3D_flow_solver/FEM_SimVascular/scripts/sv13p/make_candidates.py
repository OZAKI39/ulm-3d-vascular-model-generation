"""Choose a single documented PC profile only after reading live help."""
import json,re,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3p';C=ROOT/'configs/sv1_3p/candidates'
name=sys.argv[1];policy=json.loads((C.parent/'policy.json').read_text());outer=policy['baseline_options'].split(' -pc_type asm')[0]
if name=='HYPRE_BASELINE':options=policy['baseline_options'];build='svmp_hypre_build';proof='Mandatory short parity smoke on new PETSc binary.'
elif name=='P2':
 help=(ROOT/'logs/sv1_3p/remote/hypre_ilu_standalone_help.log').read_text()
 profile='-pc_type hypre -pc_hypre_type ilu -pc_hypre_ilu_type Block-Jacobi-ILUk -pc_hypre_ilu_level 0 -pc_hypre_ilu_maxiter 1 -pc_hypre_ilu_tol 0 -pc_hypre_ilu_tri_solve false -pc_hypre_ilu_local_reordering false'
 for option in re.findall(r'(?:^|\s)(-pc_\w+)',profile):assert option in help,option
 assert 'Block-Jacobi-ILUk' in help
 options=outer+' '+profile;build='svmp_hypre_build';proof='Pinned hypre 3.1.0 ILU GPU table: BJ-ILU0 setup+solve supported; iterative triangular solves; no reordering; actual help verified. Default Jacobi sweep counts retained; no profile scan.'
elif name=='P3':
 help=(ROOT/'logs/sv1_3p/remote/hypre_boomeramg_standalone_help.log').read_text()
 profile='-pc_type hypre -pc_hypre_type boomeramg -pc_hypre_boomeramg_coarsen_type PMIS -pc_hypre_boomeramg_interp_type ext+i -pc_hypre_boomeramg_relax_type_all l1scaled-Jacobi -pc_hypre_boomeramg_relax_type_coarse l1scaled-Jacobi -pc_hypre_boomeramg_max_iter 1 -pc_hypre_boomeramg_tol 0'
 for option in re.findall(r'(?:^|\s)(-pc_\w+)',profile):assert option in help,option
 for v in ('PMIS','ext+i','l1scaled-Jacobi'):assert v in help,v
 options=outer+' '+profile;build='svmp_hypre_build';proof='Pinned hypre 3.1.0 GPU-supported PMIS (8), extended+i (6), l1-Jacobi (18), one V cycle. Coarse relaxation also explicitly GPU-supported; no unified memory fallback.'
elif name=='P4':
 assert json.loads((R/'gamg_matrix_audit.json').read_text())['status']=='AUDITED'
 options=outer+' -pc_type gamg';build='svmp_baseline_build';proof='Official default PCGAMG aggregation, original mixed block-size-4 matrix, no layout or physical nullspace changes.'
else:raise ValueError(name)
d=dict(PETSC_OPTIONS=options,build_report=build,profile_rationale=proof)
p=C/(name+'.json');assert not p.exists();p.write_text(json.dumps(d,indent=2)+'\n')
print(name,options)
