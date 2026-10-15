from common import *
import sys,time,resource
sys.path.insert(0,str(R))
from analysis.dual_track_v1.core_a import delta_y
T=time.perf_counter();a=json.loads((O/'STATIC_EQUILIBRIUM_AUDIT.json').read_text());d=np.load(O/'final_candidate.npz');rt=np.load(O/'reference_topology.npz');tr=json.loads((O/'INDEPENDENT_EQUILIBRIUM_CANDIDATE.json').read_text());X=d['X'];x=d['x'];b=d['bonds'];alive=d['alive'].astype(bool);s=d['stretch'];u=x-X
mid=X[b[:,0]]+.5*delta_y(X,b);mid[:,1]=(mid[:,1]+PERIOD/2)%PERIOD-PERIOD/2;cmid=x[b[:,0]]+.5*delta_y(x,b);cmid[:,1]=(cmid[:,1]+PERIOD/2)%PERIOD-PERIOD/2
labels=np.load(OLD/'topology_audit.npz')['corrected_bond_regions'];assert len(labels)==len(b)
names=['notch_front','notch_mouth','loading_ends','other_X_free_surface','interior'];regions={};qps=[0,.01,.05,.25,.5,.75,.95,.99,1.]
for k,name in enumerate(names+['periodic_Y_seam_overlap']):
 idx=np.flatnonzero(alive&((labels==k) if k<5 else (rt['periodic_image_y']!=0)));j=idx[np.argmax(s[idx])];ja=idx[np.argmax(abs(s[idx]))]
 regions[name]=dict(bonds=len(idx),max_positive_stretch=float(max(0,s[j])),max_absolute_stretch=float(abs(s[ja])),max_s_over_s00=float(s[j]/S00),reference_bond_id=int(j+1),atom_ids=b[j]+1,reference_endpoints_m=X[b[j]],current_endpoints_m=x[b[j]],reference_midpoint_m=mid[j],current_midpoint_m=cmid[j],max_absolute_stretch_bond_id=int(ja+1),quantile_levels=qps,stretch_quantiles=np.quantile(s[idx],qps),distribution_source='spatial_response.npz')
np.savez_compressed(O/'spatial_response.npz',reference_bond_id=np.arange(1,len(b)+1),bond_regions=labels,periodic_seam_tag=rt['periodic_image_y']!=0,alive=alive,stretch=s,reference_midpoint_m=mid,current_midpoint_m=cmid,atom_id=d['atom_id'],displacement_m=u)
ratio=regions['notch_front']['max_positive_stretch']/regions['loading_ends']['max_positive_stretch'];leading=ratio>1
dump(O/'NOTCH_VS_LOADING_END_STRETCH.json',dict(status='ACCEPTED_ELASTIC_RESPONSE' if a['accepted'] else 'UNACCEPTED_CANDIDATE_SPATIAL_DIAGNOSTIC_ONLY',regions=regions,notch_to_loading_end_max_stretch_ratio=ratio,loading_end_to_notch_ratio=1/ratio,notch_leads_loading_end=leading,partition='Five exclusive categories using frozen corrected finite-X physical edges at +/-2.5mm. Periodic Y seam is an overlapping tag, NOT a free surface.',original_corrected_labels_sha256=sha(OLD/'topology_audit.npz'),dynamic_first_failure_location='NOT_MEASURED'))
cs=np.sqrt(G/1000);wave=dict(c_s_m_s=cs,H_crossing_s=.001/cs,W_crossing_s=.005/cs,candidate_24ms_div_W_crossing=.024/(.005/cs),quasistatic_24ms_certified=False,formal_dynamics=0)
if a['accepted']:
 magnitude=np.linalg.norm(u,axis=1);im=int(np.argmax(magnitude));opening=float(u[rt['high'],2].mean()-u[rt['low'],2].mean());face=4e-6;U0=tr['U0_LPS_J'];U=tr['last']['U_LPS_J'];W=tr['last']['Fdotu_J'];du=U-U0
 response=dict(status='ACCEPTED_SINGLE_LOW_LOAD_ENDPOINT',max_displacement_m=float(magnitude[im]),max_displacement_atom_id=im+1,max_displacement_reference_m=X[im],max_displacement_current_m=x[im],mean_displacement_vector_m=u.mean(0),mean_displacement_magnitude_m=float(magnitude.mean()),top_mean_displacement_m=u[rt['high']].mean(0),bottom_mean_displacement_m=u[rt['low']].mean(0),mean_opening_m=opening,F_face_N=face,secant_compliance_m_N=opening/face,U0_LPS_J=U0,U_loaded_LPS_J=U,delta_U_LPS_J=du,endpoint_Fdotu_J=W,contact_U_J=tr['last']['contact_U_J'],Pi_J=tr['last']['Pi_J'],two_delta_U_over_endpoint_Fdotu=2*du/W,relative_deviation_from_half_Fdotu=2*du/W-1,scope='Single zero-to-2Pa secant, not differential compliance or a measured continuous loading curve. Endpoint Fdotu is NOT integrated path work.',GC_EFFECTIVE=None,G_SEPARATION_PROTOCOL=None,signed_topology_release_J=None,new_dynamic_crack_area_m2=None)
 # One static-shape mass/compliance estimate, not an eigenmode solve.
 shape=u/opening;Meff=float(1000*V*np.sum(shape**2));ksec=face/opening;Tproxy=2*np.pi*np.sqrt(Meff/ksec);wave.update(static_shape_effective_mass_kg=Meff,secant_generalized_stiffness_N_m=ksec,static_shape_period_proxy_s=Tproxy,candidate_24ms_div_static_shape_proxy=.024/Tproxy,proxy_scope='Assume small-strain linear elastic energy and use the normalized measured static displacement as ONE trial shape. Not a measured eigenperiod, not a rigorous bound at this nonlinear endpoint, not a damping/adiabatic criterion.')
else:
 response=dict(status='NOT_ADMITTED_RESPONSE_NOT_CERTIFIED',secant_compliance_m_N=None,mean_opening_m=None,GC_EFFECTIVE=None,G_SEPARATION_PROTOCOL=None,signed_topology_release_J=None,new_dynamic_crack_area_m2=None,diagnostic_source='INDEPENDENT_EQUILIBRIUM_CANDIDATE.json; no certified response computed')
dump(O/'LOW_LOAD_COMPLIANCE_AND_ENERGY.json',response);dump(O/'LOADING_TIMESCALE_AUDIT.json',wave)
text=f'''# Mode-I 加载时间尺度：没有动态验收\n\nG=940 Pa、ρ=1000 kg/m³，c_s={cs:.9f} m/s。H=1 mm 的剪切穿越约{1/cs:.6f} ms；W=5 mm 约{5/cs:.6f} ms。24 ms仅为W穿越时间的{wave['candidate_24ms_div_W_crossing']:.3f}倍。\n\n这些是传播尺度，不是最慢弯曲周期。开裂臂可能有更慢、更柔顺的弯曲模式；真实裂纹推进又增加动力时间尺度。不存在仅凭若干波穿越就认证准静态的证据。\n\n'''
if a['accepted']:text+=f"基于唯一已接受静态位移，令形状r_i=u_i/q，q为平均开度。单形状有效质量ρVΣ|r_i|²={Meff:.9g} kg，单点割线刚度F/q={ksec:.9g} N/m，条件性时间尺度2π√(M_eff/k_sec)={Tproxy*1e3:.6f} ms。这是利用小应变近线性假设的静态形状代理，不是本征模态、严格最慢周期或实际振动数据；没有新增载荷点、Hessian谱扫描或动力模拟。24 ms/该代理约{.024/Tproxy:.3f}，不能据此批准24 ms方案。\n\n"
text+='下一轮必须独立冻结近临界加载协议、惯性/平衡门槛、分离面积与能量审计，并在低惯性证据成立后解释断裂阻力。本轮不启动高载或24 ms动态试验。\n'
(O/'LOADING_TIMESCALE_AUDIT.md').write_text(text)
dump(O/'POSTPROCESS_RESOURCES.json',dict(wall_s=time.perf_counter()-T,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,new_solves=0));print('RESPONSE',response,flush=True)
