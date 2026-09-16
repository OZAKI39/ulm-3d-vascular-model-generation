"""Independent SI/convention/matrix audit of frozen first package outputs."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import h5py
from classical_wall_reference import normal_brenner,published_wall_asymptotics
R=Path(__file__).resolve().parents[1]
P=json.loads((R/'WALL_REFERENCE_CONVENTION.json').read_text())
LABELS=['RMBW_SPHERE','RMBW_LUBRICATION','PYSTOKES']
MODES={'normal_TT':(2,2),'tangential_TT':(0,0),'TR_x_Ty':(0,4),'RT_y_Fx':(4,0),'TR_y_Tx':(1,3),'RT_x_Fy':(3,1),'RR_parallel':(3,3),'RR_normal':(5,5)}
def dump(name,obj): (R/name).write_text(json.dumps(obj,indent=2,allow_nan=True)+'\n')
def csvout(name,rows):
    with (R/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def rel(a,b):return float(np.max(np.abs(a-b)/np.maximum(np.maximum(np.abs(a),np.abs(b)),1e-30)))
def classify(x):return 'CLOSE_AGREEMENT' if x<.02 else 'MODERATE_DIFFERENCE' if x<=.1 else 'LARGE_DIFFERENCE'

def main():
    records=[];matrows=[];resrows=[];lookup={};Mlist=[];Rlist=[];normM=[];normR=[]
    inputs={label:[json.loads(s) for s in (R/'raw'/(label+'.jsonl')).read_text().splitlines()] for label in LABELS}
    expected=len(P['radii_m'])*len(P['gaps_h_over_a'])*len(P['viscosities_Pa_s'])
    assert all(len(rows)==expected for rows in inputs.values())
    theory=[]
    for e in P['gaps_h_over_a']:
        exact,n=normal_brenner(e,60);check,n2=normal_brenner(e,80);assert abs(exact-check)<=1e-12*exact
        t={'epsilon':e,'normal_resistance_brenner_over_bulk':exact,'normal_mobility_brenner_over_bulk':1/exact,'series_terms_60digits':n,'series_terms_80digits':n2,'precision_repeat_relative':abs(exact-check)/exact}
        t.update(published_wall_asymptotics(e));theory.append(t)
    csvout('validation/CLASSICAL_THEORY_VALUES.csv',theory);tmap={t['epsilon']:t for t in theory}
    allow=np.eye(6,dtype=bool)
    for i,j in [(0,4),(4,0),(1,3),(3,1)]:allow[i,j]=True
    Q=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.] ]);S=np.zeros((6,6));S[:3,:3]=Q;S[3:,3:]=Q
    for label,rawrows in inputs.items():
        for row in rawrows:
            a,mu,e=row['radius_m'],row['viscosity_Pa_s'],row['epsilon'];MT=1/(6*np.pi*mu*a);MR=1/(8*np.pi*mu*a**3)
            d=np.sqrt([MT]*3+[MR]*3);M=np.array(row['matrix_SI']) if row['matrix_SI'] is not None else np.full((6,6),np.nan);B=M/np.outer(d,d);Ri=np.full((6,6),np.nan);C=Ri.copy()
            rec={k:row[k] for k in ['implementation','size','radius_m','gap_m','epsilon','viscosity_Pa_s','status','finite']};rec['warnings']='; '.join(row['warnings']);rec['exception']=row.get('error','')
            if row['finite']:
                TR=M[:3,3:];RT=M[3:,:3].T;absrec=float(np.max(abs(TR-RT)));relrec=absrec/max(np.max(abs(TR)),np.max(abs(RT)),1e-30)
                eig=np.linalg.eigvalsh((B+B.T)/2);condition=float(np.linalg.cond(B));zero=float(np.max(abs(B[~allow])));isot=float(np.max(abs(S@B@S.T-B)));invstatus='CONDITION_REJECTED'
                valid=bool(relrec<=1e-12 and eig.min()>-1e-12 and zero<=1e-13 and isot<=1e-12)
                if condition<1e12:
                    C=np.linalg.solve(B,np.eye(6));Ri=C/np.outer(d,d);invstatus='PHYSICALLY_VALID' if valid else 'ALGEBRAIC_DIAGNOSTIC_ONLY_INVALID_M'
                rec.update(reciprocity_absolute_SI=absrec,reciprocity_relative=float(relrec),reciprocity_status='PASS' if relrec<=1e-12 else 'FAIL',scaled_min_eigenvalue=float(eig.min()),positivity_status='PASS' if eig.min()>-1e-12 else 'FAIL',scaled_condition_number=condition,unexpected_scaled_element_max=zero,planar_rotation_symmetry_error=isot,planar_symmetry_status='PASS' if isot<=1e-12 and zero<=1e-13 else 'FAIL',physical_matrix_valid=valid,resistance_status=invstatus,inverse_residual=float(np.max(abs(B@C-np.eye(6)))) if np.all(np.isfinite(C)) else np.nan)
            else:rec.update(reciprocity_absolute_SI=np.nan,reciprocity_relative=np.nan,reciprocity_status='UNAVAILABLE',scaled_min_eigenvalue=np.nan,positivity_status='UNAVAILABLE',scaled_condition_number=np.nan,unexpected_scaled_element_max=np.nan,planar_rotation_symmetry_error=np.nan,planar_symmetry_status='UNAVAILABLE',physical_matrix_valid=False,resistance_status='UNAVAILABLE',inverse_residual=np.nan)
            for mode,(i,j) in MODES.items():rec['Mhat_'+mode]=float(B[i,j]);rec['Rhat_'+mode]=float(C[i,j])
            rec['normal_vs_Brenner_relative']=abs(rec['Rhat_normal_TT']/tmap[e]['normal_resistance_brenner_over_bulk']-1)
            records.append(rec);lookup[label,row['size'],e,mu]=(M,B,Ri,C,rec);Mlist.append(M);Rlist.append(Ri);normM.append(B);normR.append(C)
            for i in range(6):
                for j in range(6):
                    info={k:rec[k] for k in ['implementation','size','epsilon','viscosity_Pa_s']};info.update(row=i,column=j,physical_matrix_valid=rec['physical_matrix_valid'])
                    matrows.append(dict(info,value_SI=float(M[i,j]),normalized=float(B[i,j])))
                    resrows.append(dict(info,value_SI=float(Ri[i,j]),normalized=float(C[i,j]),resistance_status=rec['resistance_status']))
    csvout('WALL_REFERENCE_RAW_RESULTS.csv',records);csvout('WALL_MOBILITY_MATRIX.csv',matrows);csvout('WALL_RESISTANCE_MATRIX.csv',resrows)
    for name,values,normalized in [('MOBILITY',Mlist,normM),('RESISTANCE',Rlist,normR)]:
        with h5py.File(R/('WALL_REFERENCE_'+name+'_MATRICES.h5'),'w') as h:
            h.create_dataset('matrix_SI',data=values);h.create_dataset('work_conjugate_normalized_matrix',data=normalized)
            for k in ['implementation','size','epsilon','radius_m','gap_m','viscosity_Pa_s','physical_matrix_valid']:
                vals=[r[k] for r in records];h.create_dataset(k,data=vals,dtype=h5py.string_dtype('utf-8') if isinstance(vals[0],str) else None)
            h.attrs['convention_sha256']=hashlib.sha256((R/'WALL_REFERENCE_CONVENTION.json').read_bytes()).hexdigest();h.attrs['invalid_matrix_inverse_policy']='Diagnostic only; inspect physical_matrix_valid';h.attrs['reference_repaired']=False
    comparison=[]
    for label in ['RMBW_LUBRICATION','RMBW_SPHERE']:
        for size in P['radii_m']:
            for e in P['gaps_h_over_a']:
                ar=lookup[label,size,e,.001][4];br=lookup['PYSTOKES',size,e,.001][4]
                for mode in MODES:
                    av=ar['Mhat_'+mode];bv=br['Mhat_'+mode];diff=abs(av-bv)/max(abs(av),abs(bv),1e-30)
                    comparison.append({'reference_A':label,'reference_B':'PYSTOKES','size':size,'epsilon':e,'mode':mode,'A_normalized_mobility':av,'B_normalized_mobility':bv,'symmetric_relative_difference':diff,'agreement_label':classify(diff),'A_physical_matrix_valid':ar['physical_matrix_valid'],'B_physical_matrix_valid':br['physical_matrix_valid'],'comparison_role':'DESCRIPTION_ONLY_NOT_ACCURACY_GATE'})
    csvout('WALL_REFERENCE_CROSS_COMPARISON.csv',comparison);csvout('NEAR_WALL_COMPARISON.csv',[r for r in comparison if r['epsilon']<=.2])
    scaling={}
    for label in LABELS:
        radial=[];visc=[]
        for e in P['gaps_h_over_a']:
            ref=lookup[label,'d50',e,.001]
            for size in ['d10','d90']:
                oth=lookup[label,size,e,.001];radial.append({'epsilon':e,'size':size,'mobility_max_relative':rel(ref[1],oth[1]),'resistance_max_relative':rel(ref[3],oth[3])})
            for size in P['radii_m']:
                ref=lookup[label,size,e,.001]
                for mu in [.0005,.002]:
                    oth=lookup[label,size,e,mu];visc.append({'epsilon':e,'size':size,'mu':mu,'mobility_inverse_mu_error':rel(ref[0],oth[0]*mu/.001),'resistance_mu_error':rel(ref[2],oth[2]*.001/mu)})
        mr=max(max(x['mobility_max_relative'],x['resistance_max_relative']) for x in radial);mv=max(max(x['mobility_inverse_mu_error'],x['resistance_mu_error']) for x in visc)
        scaling[label]={'radius_status':'PASS' if mr<=1e-12 else 'FAIL','viscosity_status':'PASS' if mv<=1e-12 else 'FAIL','max_radius_relative_error':mr,'max_viscosity_relative_error':mv,'radius_details':radial,'viscosity_details':visc}
    dump('WALL_REFERENCE_SCALING_AUDIT.json',scaling)
    near=[]
    for label in LABELS:
        for e in P['gaps_h_over_a']:
            r=lookup[label,'d50',e,.001][4];t=tmap[e]
            near.append({'implementation':label,'epsilon':e,'normal_R_over_bulk':r['Rhat_normal_TT'],'normal_R_times_epsilon':r['Rhat_normal_TT']*e,'normal_vs_Brenner_relative':r['normal_vs_Brenner_relative'],'Ytt_R_over_6pi':r['Rhat_tangential_TT'],'Ytr_R_over_6pi':r['Rhat_TR_x_Ty']*np.sqrt(4/3),'Yrr_R_over_6pi':r['Rhat_RR_parallel']*4/3,'Xrr_R_over_6pi':r['Rhat_RR_normal']*4/3,'published_Xtt':t['Xtt'],'published_Ytt':t['Ytt'],'published_Ytr':t['Ytr'],'published_Yrr':t['Yrr'],'published_Xrr':t['Xrr'],'asymptotic_comparison_domain':'Xtt/Ytr/Yrr <=0.1; Ytt/Xrr <=0.01, outside are unused diagnostics'})
    csvout('validation/ASYMPTOTIC_COMPARISON.csv',near)
    summary={}
    for label in LABELS:
        rr=[r for r in records if r['implementation']==label];main=[r for r in rr if r['size']=='d50' and r['viscosity_Pa_s']==.001];low=[r for r in main if r['epsilon']<=.005];far=[r for r in main if r['epsilon']>=5]
        summary[label]={'formal_rows':len(rr),'finite_rows':sum(r['finite'] for r in rr),'exception_rows':sum(r['status']=='EXCEPTION' for r in rr),'negative_mobility_rows':sum(r['scaled_min_eigenvalue']<-1e-12 for r in rr),'reciprocity_max_relative':max(r['reciprocity_relative'] for r in rr),'reciprocity_max_absolute_SI':max(r['reciprocity_absolute_SI'] for r in rr),'min_scaled_eigenvalue':min(r['scaled_min_eigenvalue'] for r in rr),'max_scaled_condition_number':max(r['scaled_condition_number'] for r in rr),'max_planar_rotation_symmetry_error':max(r['planar_rotation_symmetry_error'] for r in rr),'all_physical_matrix_valid':all(r['physical_matrix_valid'] for r in rr),'max_normal_error_vs_Brenner':max(r['normal_vs_Brenner_relative'] for r in main),'normal_error_at_smallest_gap':main[0]['normal_vs_Brenner_relative'],'normal_monotonic_with_gap':bool(np.all(np.diff([r['Mhat_normal_TT'] for r in main])>=-1e-12)),'normal_R_epsilon_at_smallest_gap':low[0]['Rhat_normal_TT']*low[0]['epsilon'],'far_field_diagonal_distance_to_bulk':[{'epsilon':r['epsilon'],'max_abs_diagonal_minus1':max(abs(r['Mhat_'+m]-1) for m in ['normal_TT','tangential_TT','RR_parallel','RR_normal'])} for r in far],'radius_status':scaling[label]['radius_status'],'viscosity_status':scaling[label]['viscosity_status']}
    dump('validation/MATRIX_AUDIT_SUMMARY.json',summary)
    print(json.dumps(summary,indent=2));print('FIRST_FORMAL_OUTPUTS_PRESERVED',sum(map(len,inputs.values())))
if __name__=='__main__':main()
