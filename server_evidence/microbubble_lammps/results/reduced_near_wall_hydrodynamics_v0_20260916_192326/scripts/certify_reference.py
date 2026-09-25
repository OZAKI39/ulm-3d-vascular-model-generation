"""Source-gated CF2003 audit and bounded reference freeze. No flow simulation."""
import csv
import hashlib
import json
import math
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def horner(c,x):
    p=c[-1]
    for a in c[-2::-1]: p=p*x+a
    return p
def bernstein(c,lo,hi):
    n=len(c)-1
    power=[sum(c[j]*Decimal(math.comb(j,k))*lo**(j-k)*(hi-lo)**k
               for j in range(k,n+1)) for k in range(n+1)]
    return [sum(power[k]*Decimal(math.comb(i,k))/Decimal(math.comb(n,k))
                for k in range(i+1)) for i in range(n+1)]

def main():
    source=json.loads((ROOT/'validation/CF2003_SOURCE_VERIFICATION.json').read_text())
    assert all(source[k] for k in ['table17_verified','corrigendum_checked','normalization_visually_verified'])
    rows=list(csv.DictReader((ROOT/'reference/CF2003_TABLE17_VERIFICATION.csv').open()))
    assert len(rows)==20 and all(r['match']=='YES' for r in rows)
    literals=[[r[key] for r in rows] for key in ['u_final','omega_final']]
    coefs=[list(map(float,c)) for c in literals]
    required=np.array([.001,.002,.003202,.005,.01,.02,.05,.08,.1,.2])
    grid=np.unique(np.r_[np.geomspace(.001,.2,10001),required])
    vals=np.array([[1/horner(c,math.log(e)) for c in coefs] for e in grid])
    assert np.isfinite(vals).all() and (vals>0).all() and (vals<1).all()
    assert (np.diff(vals,axis=0)>0).all()
    ratio=vals[:,1]/(2*(1+grid)*vals[:,0])
    assert np.all(np.diff(ratio)<0) and ratio[-1]<ratio[0]<.5676
    corrected=1/(.6600-.2693*np.log(grid))
    old_u=.7431/(.6376-.200*np.log(grid))
    old_o=.8436/(.6376-.200*np.log(grid))
    bmask=grid<.005
    corrected_error=np.abs(corrected-vals[:,0])/vals[:,0]
    # The user specified approximate 1e-3 precision, not a rigorous uniform
    # bound. Report its small exceedance near .005; never tune primary data.
    required_small_errors=[abs(1/(.6600-.2693*math.log(e))-1/horner(coefs[0],math.log(e)))
                           /(1/horner(coefs[0],math.log(e))) for e in required if e<.005]
    # Strong continuous-interval sanity check: Bernstein convex-hull bounds,
    # evaluated at 70 decimal digits on 64 subdivisions in ln(epsilon).
    bounds=[]; fp_error=0.
    with localcontext() as ctx:
        ctx.prec=70
        left=Decimal('.001').ln();right=Decimal('.2').ln()
        for cs in literals:
            c=list(map(Decimal,cs));der=[Decimal(i)*c[i] for i in range(1,len(c))]
            pmin=Decimal('Infinity');dmax=Decimal('-Infinity')
            for i in range(64):
                lo=left+(right-left)*Decimal(i)/64;hi=left+(right-left)*Decimal(i+1)/64
                pmin=min(pmin,min(bernstein(c,lo,hi)))
                dmax=max(dmax,max(bernstein(der,lo,hi)))
            assert pmin>1 and dmax<0
            bounds.append({'inverse_polynomial_lower_bound':str(pmin),'derivative_upper_bound':str(dmax)})
        for e in grid[::10]:
            x=Decimal(str(float(e))).ln()
            for cf,cs in zip(coefs,literals):
                exact=Decimal(1)/horner(list(map(Decimal,cs)),x)
                fp_error=max(fp_error,abs(1/horner(cf,math.log(e))-float(exact)))
    assert fp_error<1e-11
    gcb_table=[.3966,.4268]
    history_rel=[abs(float(v)-t)/t for v,t in zip(vals[0],gcb_table)]
    # Table 3 is a low-order historical diagnostic, not an acceptance oracle.
    historical={'epsilon':.001,'GCB_table3_transcribed':gcb_table,'CF2003':vals[0].tolist(),
                'relative_differences':history_rel,'classification':'HISTORICAL_DIAGNOSTIC',
                'rounding_half_unit':.00005,'explanation':'Residual includes old asymptotic error, not rounding alone; no coefficient tuning',
                'old_rotation_max_relative_error':float(np.max(abs(old_o-vals[:,1])/vals[:,1]))}
    with (ROOT/'raw/CF2003_REFERENCE_DENSE_GRID.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['epsilon','FU','FOMEGA','Omega_a_over_U','GCB_corrected_FU','GCB_old_FU','GCB_old_FOMEGA'])
        w.writerows(zip(grid,*vals.T,ratio,corrected,old_u,old_o))
    gridpath=ROOT/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.csv'
    with gridpath.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['epsilon','FU','FOMEGA','Omega_a_over_U'])
        for e in required:
            u,o=[1/horner(c,math.log(e)) for c in coefs];w.writerow([e,u,o,o/(2*(1+e)*u)])
    record={'format':'CF2003_FREE_SHEAR_REFERENCE_V0','paper_DOI':'10.1093/qjmam/56.3.381',
            'corrigendum_DOI':'10.1093/qjmam/hbs012','paper_sha256':sha(ROOT/'reference/materials/CF2003_user_supplied_published.pdf'),
            'corrigendum_pdf_sha256':None,'corrigendum_sha256':sha(ROOT/'reference/materials/corrigendum_preview.gif'),
            'corrigendum_material':'Official complete single-page publisher GIF; PDF unavailable',
            'correction_status':'DOES_NOT_AFFECT_CLOSURE','corrigendum_checked':True,
            'normalization':{'FU':'Ux/(kappa*l)','FOMEGA':'Omega_y/(kappa/2)','far_field_convention':[1,1],
                             'rotation_sign':'n=+z, g=+kappa*x => n cross g=+kappa*y; right-handed'},
            'gap_definition':'h=l-a; epsilon=h/a; l/a=1+epsilon; d=l',
            'coefficient_order':'ascending powers j=0..19','u_coefficients':coefs[0],'omega_coefficients':coefs[1],
            'u_decimal_literals':literals[0],'omega_decimal_literals':literals[1],
            'certified_epsilon_domain':[.001,.2],'formula':'FU=1/P_U(ln(epsilon)); FOMEGA=1/P_Omega(ln(epsilon))',
            'numeric_evaluator_method':'Horner; natural logarithm; no extrapolation',
            'source_pages':{'definitions':383,'equations_6_1_to_6_8':404,'table17_and_eq6_9':407,'figure7':408},
            'source_accuracy':'Paper reports fitted small-gap approximation error below 1e-11, pp404/407; not an all-gap exact solution',
            'domain_basis':'p407 small/large approximation crossover 0.4; project restricts to [0.001,0.2]',
            'REFERENCE_RESOLUTION_PASS':True,'REFERENCE_SOURCE_VALIDATED':True,
            'verification_csv_sha256':sha(ROOT/'reference/CF2003_TABLE17_VERIFICATION.csv')}
    path=ROOT/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.json';path.write_text(json.dumps(record,indent=2)+'\n')
    (ROOT/'reference/CF2003_FREE_SHEAR_REFERENCE_V0.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in [path,gridpath]))
    validation={'REFERENCE_RESOLUTION_PASS':True,'REFERENCE_SOURCE_VALIDATED':True,'certified_epsilon_domain':[.001,.2],
                'dense_grid_points':len(grid),'finite_positive_below_one':True,'FU_FOMEGA_strictly_increasing':True,
                'rolling_ratio_decreases_with_gap':True,'rolling_ratio_endpoints':[float(ratio[0]),float(ratio[-1])],
                'contact_limit_used_as_exact_at_min_gap':False,'Horner_vs_70_digit_max_abs':fp_error,
                'continuous_domain_sanity':{'method':'70-digit Bernstein coefficient hull on 64 log-gap subintervals; numeric bounds, not interval-arithmetic proof', 'bounds':bounds},
                'corrected_first_order_max_relative_error_below_005':float(corrected_error[bmask].max()),
                'corrected_first_order_published_relative_precision':1e-3,
                'corrected_first_order_required_small_grid_relative_errors':required_small_errors,
                'strict_uniform_1e_3_diagnostic_pass':bool(corrected_error[bmask].max()<=1e-3),
                'corrected_first_order_interpretation':'Order 1e-3 historical consistency; strict 1e-3 not established near .005. No coefficient adjustment and no physical tolerance invented.',
                'historical_cross_checks':historical}
    (ROOT/'validation/CF2003_REFERENCE_VALIDATION.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2));print(gridpath.read_text())

if __name__=='__main__':main()
