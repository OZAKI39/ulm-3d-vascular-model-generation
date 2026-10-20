from common import *
from graph_ops import topology,boundary
import csv
def main():
 validation=json.loads((O/'SMALL_SPECIMEN_VALIDATION.json').read_text());comparisons=[]
 for case in validation['cases']:
  root=O/validation.get('benchmark_root','benchmark')/case['case'];r=load_npz(root/'provenance/reference.npz');b=r['bonds'];bt=r['bond_types'];n=int(r['clot_count']);X=r['X'];cut=r['initial_cut'];mid=.5*(X[b[:,0]]+X[b[:,1]]);tip=np.hypot(mid[:,0]+DX,mid[:,2])<=HORIZON;grip=(X[:n,2].max()-mid[:,2])<=HORIZON;interface=(bt==2);bulk=bt==1
  soft=np.full(len(b),np.inf);failure=np.full(len(b),np.inf)
  for line in (root/'damage_events.jsonl').open():
   e=json.loads(line);ids=np.asarray(e['bond_ids'])-1;soft[ids]=np.minimum(soft[ids],e['time_s']);dd=np.asarray(e['d']);failure[ids[dd==1]]=np.minimum(failure[ids[dd==1]],e['time_s'])
  final=load_npz(sorted((root/'keys').glob('*.npz'))[-1]);natural=(final['d']==1)&~cut;cross=r['cross'];ci=r['cross_ids'];hit=r['cross_hits'];remaining=cross&~cut
  trace=[]
  for path in sorted((root/'keys').glob('*.npz')):
   s=load_npz(path);ll,ss,att=topology(r,s['d']);f=(s['d']==1)&~cut;softened=(s['d']>0)&~cut;crossfail=remaining&f;validhits=hit[(s['d'][ci]==1)&~cut[ci]];regions=[]
   for mask,label in [(tip&bulk,'notch_delta'),(grip&bulk,'top_grip_delta'),(bulk&~tip&~grip,'other_bulk'),(interface,'interface')]:regions.append(dict(region=label,original_natural_eligible=int(np.sum(mask&~cut)),softened=int(np.sum(mask&softened)),failed=int(np.sum(mask&f))))
   trace.append(dict(step=int(s['step']),time_s=float(s['time']),regions=regions,cross_plane_live=int(np.sum(cross&(s['d']<1))),new_cross_plane_failed=int(crossfail.sum()),furthest_failed_cross_intersection_x_m=float(validhits[:,0].max()) if len(validhits) else None,body_components=len(ss),body_component_sizes=ss,anchor_attached=att))
  early=soft[bulk&~cut];on=early.min();win=on+50e-6;earlymask=bulk&~cut&(soft<=win);firstids=np.flatnonzero(bulk&~cut&(soft==on));hist=np.genfromtxt(root/'history.csv',delimiter=',',names=True)
  groupstats=[]
  for comp in case['final_components']:
   if comp['status']!='ACTUAL DETACHED':continue
   nodes=np.asarray(comp['atom_ids'])-1;outside,inside=boundary(r,nodes);bulkinside=inside[bt[inside]==1];groupstats.append(dict(atom_ids=nodes+1,nodes=len(nodes),equivalent_diameter_m=(6*len(nodes)*V/np.pi)**(1/3),internal_reference_bulk=len(bulkinside),internal_alive_bulk=int(np.sum(final['d'][bulkinside]<1)),internal_alive_fraction=float(np.mean(final['d'][bulkinside]<1)) if len(bulkinside) else None,external_reference_bonds=len(outside),external_alive=int(np.sum(final['d'][outside]<1)),external_bond_ids=outside+1))
  case['spatial_analysis']=dict(first_bulk_soft_time_s=float(on) if np.isfinite(on) else None,first_bulk_soft_bond_ids=firstids+1,first_bulk_soft_midpoints_m=mid[firstids],first_50us_softened_bulk=int(earlymask.sum()),first_50us_notch_fraction=float(np.mean(tip[earlymask])) if earlymask.any() else None,first_50us_grip_fraction=float(np.mean(grip[earlymask])) if earlymask.any() else None,final_bulk_failed=int(np.sum(natural&bulk)),final_notch_failed=int(np.sum(natural&bulk&tip)),final_grip_failed=int(np.sum(natural&bulk&grip)),final_new_cross_plane_failed=int(np.sum(natural&remaining)),initial_remaining_cross_plane_bonds=int(remaining.sum()),final_remaining_cross_plane_bonds=int(np.sum(cross&(final['d']<1))),accepted_keys=trace,free_components=groupstats,region_note='Regions notch_delta and top_grip_delta may overlap only if geometry permits; counts explicitly reported. Cross-plane failure advance is a set of original bond intersections, not a reconstructed continuous crack surface.')
  case['quantitative_response']=dict(peak_abs_top_reaction_N=float(np.max(abs(hist['reaction_N']))),peak_reaction_time_s=float(hist['time_s'][np.argmax(abs(hist['reaction_N']))]),final_W_J=float(hist['W_J'][-1]),final_U_J=float(hist['U_J'][-1]),final_K_J=float(hist['K_J'][-1]),final_D_J=float(hist['D_J'][-1]),final_contact_U_J=float(hist['Uc_J'][-1]),final_R_J=float(hist['R_J'][-1]),final_K_over_U=float(hist['K_J'][-1]/max(hist['U_J'][-1],1e-30)),max_contact_pairs=int(hist['contact_pairs'].max()))
  comparisons.append(dict(case=case['case'],first_bulk_soft_time_s=case['spatial_analysis']['first_bulk_soft_time_s'],first_50us_notch_fraction=case['spatial_analysis']['first_50us_notch_fraction'],bulk_failed=case['spatial_analysis']['final_bulk_failed'],interface_failed=case['last_row']['interface_failed'],largest_free=case['maximum_free_nodes'],remaining_cross=case['spatial_analysis']['final_remaining_cross_plane_bonds'],**case['quantitative_response']))
 validation['direction_comparison']=comparisons;validation['comparison_limit']='Only prescribed displacement direction changed; resulting reaction/work/inertia necessarily differ. This is not a matched-energy test or calibrated Mode-I/Mode-II toughness identification. A negative result remains a completed physical experiment.'
 dump(O/'SMALL_SPECIMEN_VALIDATION.json',validation);print(json.dumps(comparisons,default=cv,indent=2))
if __name__=='__main__':main()
