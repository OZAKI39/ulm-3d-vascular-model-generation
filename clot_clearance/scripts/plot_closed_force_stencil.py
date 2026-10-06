"""Four scientific audit figures and the focused review; no presentation work."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image

BLUE='#0072B2';ORANGE='#D55E00';GREEN='#009E73';INK='#233746';GRAY='#8495a3'


def make_figures(summary,destination):
    destination.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':10,
                         'axes.titlesize':12,'legend.fontsize':9,'pdf.fonttype':42,
                         'axes.spines.top':False,'axes.spines.right':False,'text.color':INK,
                         'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
    centers=[c['observations'][0] for c in summary['cases']];names=[f"N{c['n']}" for c in summary['cases']]
    force_tol=summary['configuration']['primary_force_tolerance'];coef_tol=summary['configuration']['coefficient_match_tolerance']
    oracle=summary['oracle'];scale=oracle['local_force_scale_N'];figures=[]

    def export(fig,name):
        fig.savefig(destination/f'{name}.png',dpi=600,facecolor='white')
        fig.savefig(destination/f'{name}.pdf',facecolor='white')
        with Image.open(destination/f'{name}.png') as im:
            assert min(im.info['dpi'])>599
            dimensions=list(im.size);im.thumbnail((1500,1100));im.save(destination/f'{name}_preview.png')
        figures.append(dict(name=name,formats=['png','pdf'],png_dpi=600,dimensions_px=dimensions))
        plt.close(fig)

    fig,ax=plt.subplots(figsize=(8,5.2));fig.subplots_adjust(left=.13,right=.95,bottom=.23,top=.83)
    y=[p['force_local_normalized'] for p in centers]
    ax.semilogy(range(3),y,'o-',color=BLUE,lw=1.4,label='Center force (reconstructed affine state)')
    ax.axhline(force_tol,color=ORANGE,ls='--',lw=1.2,label='Predeclared force tolerance: 1e-10')
    ax.axhline(oracle['center_force_local_normalized'],color=GRAY,ls=':',lw=1.5,label='Infinite-lattice oracle noise')
    for i,v in enumerate(y):ax.annotate(f'{v:.2e}',(i,v),xytext=(0,10),textcoords='offset points',ha='center',fontsize=9)
    ax.set(xticks=range(3),xticklabels=[f'{n}\n{j+1}-hop center closure' for j,n in enumerate(names)],
           ylabel=r'$\|\mathbf{f}_{center}\| / (|P_{zz}|h^2)$',ylim=(1e-18,1e-9),xlim=(-.3,2.3))
    ax.grid(axis='y',alpha=.2);ax.legend(loc='upper center',fontsize=8.4,framealpha=1)
    fig.suptitle('Center forces vanish before coefficient closure',y=.96,fontsize=15,weight='bold')
    fig.text(.5,.055,'All centers pass the force-only gate. Symmetry prevents identifying a minimum sufficient closure.',ha='center',fontsize=9)
    export(fig,'closure_depth_force_residual')

    fig,ax=plt.subplots(figsize=(8,5.4));fig.subplots_adjust(left=.13,right=.94,bottom=.24,top=.79)
    for key,label,color,marker in [('omega0',r'$\omega_0$: required D0–D2',BLUE,'o'),
                                   ('K','K: required D0–D2',GREEN,'s'),('Z','Z: required D0–D1',ORANGE,'^')]:
        values=[p['coefficients']['max_relative_error'][key] for p in centers]
        ax.semilogy(range(3),values,marker+'-',color=color,lw=1.4,ms=6,label=label)
    ax.axhline(coef_tol,color=GRAY,ls='--',lw=1.2,label='Predeclared coefficient tolerance: 1e-12')
    zerr=centers[-1]['coefficients']['max_relative_error']['Z']
    ax.annotate(f'N27 Z: {zerr:.3e}',xy=(2,zerr),xytext=(.50,3e-10),textcoords='data',
                arrowprops=dict(arrowstyle='-',color=ORANGE),fontsize=9,color=ORANGE)
    ax.set(xticks=range(3),xticklabels=names,ylabel='Maximum relative coefficient error',ylim=(1e-16,3),xlim=(-.2,2.25))
    ax.grid(axis='y',alpha=.2);ax.legend(loc='lower left',bbox_to_anchor=(-.01,1.025),ncol=2,frameon=False,fontsize=8.6)
    fig.suptitle('Coefficient convergence to the infinite lattice',y=.965,fontsize=15,weight='bold')
    fig.text(.5,.12,'Matrix errors use Frobenius norms. K includes the local F reconstruction dependencies.',ha='center',fontsize=9)
    fig.text(.5,.067,'The raw reconstructed-stress Z comparison retains its small exceedance of the fixed tolerance.',ha='center',fontsize=9,color=ORANGE)
    export(fig,'coefficient_convergence_to_bulk')

    n27=summary['cases'][-1]
    rows=sorted([p for p in n27['observations'] if p['label'].endswith('_bottom') or p['label']=='center'],key=lambda p:p['boundary_distance_over_delta'])
    x=np.array([p['boundary_distance_over_delta'] for p in rows])
    fig,axes=plt.subplots(1,2,figsize=(12,5.4));fig.subplots_adjust(left=.075,right=.98,bottom=.25,top=.74,wspace=.28)
    for key,label,color,marker in [('force','Total',INK,'o'),('direct','Direct',BLUE,'s'),('Z_correction','Z correction',ORANGE,'^')]:
        values=np.array([np.linalg.norm(p['channels']['original'][key])/scale for p in rows])
        axes[0].semilogy(x,values,marker+'-',label=label,color=color,lw=1.25,ms=5)
    axes[0].axhline(force_tol,color=GRAY,ls='--',lw=1,label='1e-10 tolerance')
    axes[0].set(ylabel=r'Force norm / $(|P_{zz}|h^2)$',title='(a) Direct and Z channels')
    for weight,label,color,marker in [('original','Original force',INK,'o'),('bulk_weight','Bulk-weight component',BLUE,'s'),('normalization_excess','Normalization-excess component',ORANGE,'^')]:
        values=np.array([p['channels'][weight]['force'][2]/scale for p in rows])
        axes[1].plot(x,values,marker+'-',label=label,color=color,lw=1.25,ms=5)
    axes[1].set_yscale('symlog',linthresh=force_tol,linscale=.5)
    axes[1].axhline(0,color=GRAY,lw=.6)
    axes[1].set(ylabel=r'Signed $f_z / (|P_{zz}|h^2)$',title='(b) Exact weight decomposition')
    for ax in axes:
        ax.set_xlabel(r'Nearest boundary distance / $\delta$');ax.grid(axis='y',alpha=.2)
        ax.legend(loc='lower left',bbox_to_anchor=(-.01,1.14),frameon=False,ncol=2,fontsize=8)
        ax.set_xticks([1,2,3]);ax.set_xlim(.85,3.38)
    fig.suptitle('N27: off-center residuals and the fully closed center',fontsize=16,weight='bold',y=.97)
    fig.text(.5,.12,'Actual sampled depths: 4h, 5h, 8h, 9h, 12h and center (13h). Reflected pairs are retained in the data.',ha='center',fontsize=9)
    fig.text(.5,.065,'All components reconstruct the original model; no weight is replaced and no boundary fix is applied.',ha='center',fontsize=9)
    export(fig,'n27_force_vs_boundary_distance')

    fig=plt.figure(figsize=(11,5.4));ax=fig.add_axes([.055,.42,.89,.34]);ax.set_axis_off()
    levels=centers[-1]['closure']['levels'];positions=np.linspace(.105,.895,4)
    for i,(position,level) in enumerate(zip(positions,levels)):
        ax.add_patch(FancyBboxPatch((position-.095,.06),.19,.84,boxstyle='round,pad=.009',
                                   transform=ax.transAxes,facecolor='#f5f8fa',edgecolor=GRAY,lw=.8))
        ax.text(position,.77,f'D{i}',transform=ax.transAxes,ha='center',fontsize=16,weight='bold')
        count=level['required_unique_nodes']
        ax.text(position,.54,f"{count:,} {'node' if count==1 else 'nodes'}",transform=ax.transAxes,ha='center',fontsize=12)
        ax.text(position,.32,f"max radius: {level['max_required_radius_m']*1e6:g} µm",transform=ax.transAxes,ha='center',fontsize=9.5)
        ax.text(position,.14,f"max coordinate offset: {level['max_required_coordinate_offset_cells']}h",transform=ax.transAxes,ha='center',fontsize=8.5)
        if i<3:ax.annotate('',xy=(positions[i+1]-.105,.5),xytext=(position+.105,.5),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',lw=1.1,color=INK))
    fig.suptitle('Actual three-hop coefficient dependency graph',fontsize=17,weight='bold',y=.965)
    fig.text(.5,.843,'N27 center: every required node exists and has interior volume; every D2 family is complete.',ha='center',fontsize=10.5)
    ax=fig.add_axes([.115,.16,.77,.22]);ax.axis('off')
    rows=[]
    for name,p in zip(names,centers):rows.append([name]+['CLOSED' if lev['closed_with_bulk_volumes'] else 'OPEN' for lev in p['closure']['levels'][1:]])
    tab=ax.table(cellText=rows,colLabels=['Domain','D1 + bulk volumes','D2 + bulk volumes','D3 + bulk volumes'],cellLoc='center',bbox=[0,0,1,1])
    tab.auto_set_font_size(False);tab.set_fontsize(9.5)
    for (r,c),cell in tab.get_celld().items():
        cell.set_edgecolor('#d8e0e6');cell.set_linewidth(.6)
        if r==0:cell.set_facecolor('#edf2f5')
        elif c>0:cell.get_text().set_color(GREEN if cell.get_text().get_text()=='CLOSED' else ORANGE)
    fig.text(.5,.071,'Dq = neighbors of D(q−1); sets overlap. Counts are unique nodes at each walk depth, not disjoint shells.',ha='center',fontsize=9)
    fig.text(.5,.029,'Geometric closure is verified independently of zero center force. The conservative coefficient depth is three.',ha='center',fontsize=9)
    export(fig,'force_stencil_dependency_graph')
    return figures


def write_review(summary,path,figures):
    s=summary['status'];o=summary['oracle'];scale=o['local_force_scale_N'];cases=summary['cases'];c27=cases[-1];p27=c27['observations'][0]
    e=p27['coefficients']['max_relative_error'];diag=c27['center_stress_roundoff_diagnostic'];r=summary['resources']
    center_rows='\n'.join(f"| N{c['n']} | {c['L_m']*1000:g} | {c['observations'][0]['boundary_distance_over_delta']:.6f} | {c['observations'][0]['force_norm_N']:.9e} | {c['observations'][0]['force_local_normalized']:.9e} | {s['N'+str(c['n'])+'_CLOSURE']} |" for c in cases)
    coefficient_rows='\n'.join(f"| N{c['n']} | "+' | '.join(f"{c['observations'][0]['coefficients']['max_relative_error'][k]:.9e}" for k in ('omega0','K','Z','wbar'))+' |' for c in cases)
    off=[p for p in c27['observations'] if p['label'].endswith('_bottom') or p['label']=='center'];off.sort(key=lambda p:p['boundary_distance_cells'])
    off_rows='\n'.join(f"| {p['label']} | {p['boundary_distance_cells']} | {p['boundary_distance_over_delta']:.6f} | {p['force_local_normalized']:.9e} | {p['channels']['original']['direct'][2]/scale:.9e} | {p['channels']['original']['Z_correction'][2]/scale:.9e} | {p['channels']['bulk_weight']['force'][2]/scale:.9e} | {p['channels']['normalization_excess']['force'][2]/scale:.9e} |" for p in off)
    validity='\n'.join(f"| N{v['n']} | {v['force_max_local_scaled_error']:.9e} | {v['coefficient_max_scaled_errors']['omega0']:.9e} | {v['coefficient_max_scaled_errors']['K']:.9e} | {v['coefficient_max_scaled_errors']['Z']:.9e} | {'PASS' if v['PASS'] else 'FAIL'} |" for v in summary['selective_validation'])
    maximum_decomposition=max(max(p['decomposition_errors_local_scaled'].values()) for c in cases for p in c['observations'])
    symmetry=max(p['pair_defect_local_normalized'] for c in cases for p in c['symmetry_pairs'])
    figure_links='\n'.join(f"- [{f['name']} — PNG](../../results/{figures.parent.name}/figures/{f['name']}.png) · [PDF](../../results/{figures.parent.name}/figures/{f['name']}.pdf)" for f in summary['figures'])
    lines=[]
    for k,v in s.items():
        if k not in ('minimum_depth_reason','center_force_only_first_observed_closed_hop'):lines.append(f'{k}: {v}')
    labels='\n'.join(lines)
    verdict=(f"**The algorithmically closed N27 center recovers zero nodal force at the predeclared force tolerance:** "
             f"{s['N27_CENTER_FORCE_LOCAL_NORMALIZED']:.9e} ≤ 1e-10. "
             f"**The combined scientific acceptance status is `{s['BULK_AFFINE_FORCE_PATCH_TEST']}`**, with primary root cause `{s['ROOT_CAUSE_PRIMARY']}`.")
    limitation=(f"The raw Z bulk mismatch is {e['Z']:.9e}, exceeding the fixed 1e-12 coefficient gate by {(e['Z']/1e-12-1)*100:.3f}%. "
                "This exceedance is retained. It is not relabeled PASS or replaced by a constant-stress diagnostic. The fully closed force result does not support a nonzero bulk-force defect, but the user-specified conjunction needed to upgrade the formal root-cause label has not been met.") if s['N27_COEFFICIENT_GATE']=='FAIL' else 'All predeclared coefficient and force gates pass.'
    text=rf'''# PFPD_BAQP_FULLY_CLOSED_FORCE_STENCIL_AUDIT_V1

{verdict}

{limitation}

## 1. Why d ≥ δ was not enough

An intact nodal force reads incident-neighbor K inverse and Z. Neighbor Z reads endpoint-normalized weights on its own bonds; those weights read omega0 at a second neighbor, whose family reaches a third hop. An isolated own family is not an isolated force stencil. No fracture calculation was used to infer this chain.

## 2. Infinite-lattice bulk oracle

Fixed h={summary['configuration']['h_m']:.12g} m, δ={summary['configuration']['delta_m']:.12g} m, δ/h=4.02, ε=0.003, SVK μ=1000 Pa, λ=0, quadratic kernel. The integer-offset set uses the same `0 < |offset| h ≤ δ` inclusion as the existing neighbor construction. The kernel is zero at exactly δ. There are **{o['neighbor_offset_count']}** offsets and every volume is h³.

The oracle computes omega0, K, K inverse, wbar, Z, and the incident center force using the existing r² direct formula. Pzz={o['Pzz']:.16g} Pa and the fixed local scale is **{scale:.16e} N**. Oracle force is `{np.asarray(o['center_force_N']).tolist()}` N; norm/local={o['center_force_local_normalized']:.9e}, PASS. Tensor entries, eigenvalues, offsets and weights are in [BULK_LATTICE_ORACLE.json](BULK_LATTICE_ORACLE.json).

The force scale remains |Pzz|h² for every domain and the zero-strain baseline. No local residual is normalized by Pzz L².

## 3. How the BAQP force dependency propagates

`D0={{i}}; D1=N(D0); D2=N(D1); D3=N(D2)` are explicit graph-walk sets, not disjoint shells. Infinite-lattice counts are **1, 256, 1,797, 5,821**, with maximum coordinate offsets 0, 4, 8, 12 cells and radii 0, 250, 500, 750 µm. Each finite graph is built separately and compared to its required untruncated counterpart; missing nodes and non-interior volumes are recorded.

The selective evaluator retains local Fbar, Ftilde and P reconstruction so its numerical operation follows `baqp.evaluate()`. This additionally requires K/Fbar on D2; those families still end in D3. Only incident batches of Ftilde/P exist in memory. No full N19/N27 bond state, dense matrix, phase state or equilibrium solve is allocated.

For each unordered bond (a,b), with ξ=Xb−Xa, Bi=ω Ka⁻¹ξ and Bj=−ω Kb⁻¹ξ, the implementation uses

$$
t_a=\bar\omega P\xi/|\xi|^2+Z_aB_i,\qquad
t_b=-\bar\omega P\xi/|\xi|^2+Z_bB_j,
$$
$$f_{{ab}}=(t_a-t_b)V_aV_b.$$

The lower endpoint receives +fab and the upper endpoint −fab. Both selective and full validation use lexicographic unordered bonds, including the original two endpoint-accumulation passes. No physical bond is duplicated.

## 4. Why N11 / N19 / N27 were chosen

All cubes use fixed h, L=(n−1)h, and tensor-product boundary-node volumes h³, h³/2, h³/4, h³/8. Their center distances are 5h, 9h, 13h, respectively. The largest volume-sum relative error is {max(c['volume_relative_error'] for c in cases):.3e}.

| Domain | L [mm] | Center d/δ | Center force norm [N] | Norm / local scale | Verified closure with bulk volumes |
|---|---:|---:|---:|---:|---|
{center_rows}

## 5. Selective evaluator validation

Eight center/near-center points are selected by lattice offsets, including both signs of all axes and (1,−1,−1). Force, direct and Z correction are compared to `baqp.evaluate()` and its full diagnostic decomposition. Every needed omega0 and every D0/D1 K/Z coefficient is checked, not only the central coefficient. The force error is scaled by the fixed local force scale; coefficient errors by the infinite-lattice scalar or Frobenius norm.

| Grid | Maximum force/direct/Z error | omega0 error | K error | Z error | Gate |
|---|---:|---:|---:|---:|---|
{validity}

An initial analytic-constant-P prototype exceeded the N13 Z comparison threshold (1.0895332849e-12). It was retired before N19/N27 use. Its records remain as `INITIAL_CONSTANT_P_VALIDATION_N11.json` and `INITIAL_CONSTANT_P_VALIDATION_N13.json`. The final evaluator reconstructs local F and P and passes both grids without changing the 1e-12 threshold. The initial mismatch is not hidden or treated as a physics modification.

## 6. One-hop closure result

N11 has closed D1 with interior volumes but open D2/D3. Its center passes the force-only gate while force-relevant coefficients remain strongly nonuniform. Inversion symmetry at the center can cancel a boundary-contaminated force. This is not evidence that one hop generally suffices.

## 7. Two-hop closure result

N19 closes D2 with bulk volumes, but required D3 dependencies still encounter the finite boundary. Its center again passes the force-only gate. Near-center and geometrically selected off-center results are retained in the observation table rather than using center cancellation to declare translation invariance.

## 8. Three-hop fully closed result

For the N27 center all required D0–D3 nodes exist; every required volume is h³; every D2 neighbor offset set equals the bulk oracle set. `FULL_FORCE_STENCIL_CLOSED=YES` follows from these explicit checks, not just d>3δ. Its force vector is `{np.asarray(p27['force_N']).tolist()}` N, norm **{p27['force_norm_N']:.9e} N**, norm/local **{p27['force_local_normalized']:.9e}**.

The nearest 3δ lattice point is at depth 12h, not the center at 13h. At 12h some required D3 nodes lie exactly on the boundary and have half volume; mere node existence does not satisfy the bulk-volume closure test. [DEPENDENCY_CLOSURE.json](DEPENDENCY_CLOSURE.json) records both kinds of closure.

## 9. Coefficient convergence to bulk values

omega0 and K maxima cover D0–D2 (K includes numerical F reconstruction); Z covers D0–D1; wbar covers all bonds incident to D0/D1. Direct-B K restricted to D0/D1 is separately retained in JSON. Matrices use relative Frobenius errors, avoiding division by exact zero entries. Center and representative D1 tensors are also saved.

| Domain | omega0 maximum relative error | K maximum relative error | Z maximum relative error | wbar maximum relative error |
|---|---:|---:|---:|---:|
{coefficient_rows}

At the closed N27 center the supplementary diagnostic holds the existing weights and geometry fixed and evaluates Z at the exact analytic affine P. Its maximum Z bulk error is **{diag['constant_P_Z_max_bulk_relative_error']:.9e}**. The reconstructed-P minus constant-P Z difference is **{diag['reconstructed_P_minus_constant_P_Z_max_relative_error']:.9e}**, while the maximum Fbar component error is **{diag['Fbar_D2_max_absolute_error']:.9e}**. These observations locate the small raw Z mismatch in floating-point stress reconstruction; they do not waive the original raw coefficient gate. See [NOISE_BASELINE.json](NOISE_BASELINE.json).

## 10. Does affine bulk force vanish?

**Yes, at the prescribed local force tolerance for the fully closed N27 center.** The combined acceptance rule also requires coefficients ≤1e-12 and is therefore **{s['BULK_AFFINE_FORCE_PATCH_TEST']}** in this run. These are distinct statements.

`MINIMUM_OBSERVED_CLOSURE_DEPTH_FOR_BULK_FORCE=NOT_ESTABLISHED`. The first center with zero force is already N11, but all centers possess inversion symmetry. Some off-center points at the same nominal closure depth have different residuals. This audit establishes the measured fully closed center force; it does not prove a universal minimum closure depth or prescribe a collar thickness.

Noise references: infinite-lattice norm/local={o['center_force_local_normalized']:.9e}; zero-strain maximum force={summary['noise']['zero_strain_max_force_N']:.9e} N; largest top/bottom paired vector defect/local={symmetry:.9e}. The 1e-10 force threshold was fixed before these evaluations.

## 11. Which channel transmits finite-domain influence?

The original force is decomposed as direct + Z correction. Separately, wbar = omega/omega0_inf + (wbar−omega/omega0_inf); Z is independently accumulated for each weight component using the same reconstructed P and K. These are exact summands, not alternative models. Maximum decomposition error/local={maximum_decomposition:.9e} ≤1e-12.

N27 observations below use the bottom axis; matching reflected points were also evaluated. The three approximate horizon depths use nearest lattice layers (4h, 8h, 12h). Two extra samples at 5h and 9h provide noncentral one/two-hop closure witnesses. All entries except distance and force norm are signed axial contributions normalized by the same local force scale.

| Observation | Depth [h] | d/δ | Force norm / local | Direct fz / local | Z fz / local | Bulk-weight fz / local | Normalization-excess fz / local |
|---|---:|---:|---:|---:|---:|---:|---:|
{off_rows}

Secondary mechanism label: `{s['ROOT_CAUSE_SECONDARY']}`. Nonzero Z and normalization-excess components identify transmission channels; their size does not make either term erroneous.

## 12. Final root cause classification

Primary: **`{s['ROOT_CAUSE_PRIMARY']}`**. Secondary: **`{s['ROOT_CAUSE_SECONDARY']}`**. The strict reason the primary label was not upgraded is the recorded raw Z coefficient exceedance, not a nonzero fully closed bulk force. No tolerance, old scientific result, production formula or boundary condition was changed.

Seven focused tests were run: deterministic lattice stencil; selective/full N11 equivalence; deterministic D1/D2/D3; both force decompositions; zero strain; N27 geometry/volume closure; infinite-lattice symmetry. The [test log](TARGETED_TESTS.log) contains the actual results. Full pytest and historical cleanup QA were not run. An initial new symmetry-test assertion incorrectly used a dimensional 1e-30 bound for exact-zero matrix entries; the final assertion uses the predeclared 1e-12 relative tensor norm. The initial log is retained as `TARGETED_TESTS_INITIAL.log`; no historical test or scientific acceptance limit was changed.

## 13. Implication for existing C5/C7 activation drift

The intact mechanical evidence is consistent with finite-domain multi-hop coefficients contaminating a mechanical driver field. Under the full predeclared acceptance rule, the formal root cause remains unresolved pending the coefficient-reference issue. Existing C5/C7 activation values and accepted states remain **UNCHANGED**. No claim of fracture mesh objectivity, a correct physical nonlocal BC, or an incorrect surface Dirichlet condition follows from this audit. No activation, phase, stability, Gc/Yc or CFD calculation was performed.

## 14. What should be developed next

First obtain human review of the exact-affine versus reconstructed-stress coefficient reference and its fixed tolerance; keep this run's exceedance intact. A focused affine-coefficient-reference audit can settle that numerical acceptance issue. A subsequent, separately approved boundary-development task should compare a physically specified nonlocal loading/interaction-domain formulation against its intended traction and virtual work, at fixed h and δ. Do not select a collar thickness by forcing C5/C7 activation thresholds to agree. No collar, ghost region, force-density loading or volume constraint was implemented here.

CPU workers=1, BLAS threads=1. Recorded audit computation, targeted tests and scientific figure exports: {r['wall_seconds']:.3f} s; sampled process-tree peak RSS={r['process_tree_peak_RSS_sampled_bytes']/2**30:.3f} GiB; minimum sampled available RAM fraction={r['minimum_available_fraction']:.3%}. GPU=`NOT_USED_BY_DESIGN`. The monitor samples at 50 ms; peaks are reported as sampled measurements, not exact continuous-time maxima.

Reproduce from the `clot_clearance` root using the existing environment:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/clot_closed_stencil_mpl /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m scripts.fully_closed_force_stencil_audit
```

Scientific figures only (English; PNG 600 dpi and PDF):

{figure_links}

```text
{labels}
```

Stop for human review. No next development task was started.
'''
    path.write_text(text)
