"""One actual public-API evaluation and explicit missing-bubble-input response."""
from common import *
from streaming_interface import evaluate,FLOW_CLASSES
from acoustic_streaming_shear import MU,FLUID_RHO,GAP,SPEED,WAVELENGTH
request=dict(model_label='PRESCRIBED_STOKES_SHEAR_TEMPLATE',rho_kg_m3=FLUID_RHO,mu_Pa_s=MU,gap_m=GAP,upper_speed_m_s=SPEED,wavelength_m=WAVELENGTH,time_s=.0024,local_surface_coordinates_m=[[0,-175e-6,0],[0,0,0],[0,175e-6,0]],solid_to_fluid_normals=[[0,0,1]]*3)
result=evaluate(request);shear=result['tangential_traction_Pa'];assert shear[0,0]<0 and shear[1,0]==0 and shear[2,0]>0 and not result['normal_traction_Pa'].any()
missing=evaluate(dict(model_label='BUBBLE_DRIVEN_STREAMING_MODEL',frequency_Hz=2.5e6));assert missing['status']=='STREAMING_LOAD_NOT_YET_DETERMINED' and missing['flow'] is None
dump(O/'PUBLIC_INTERFACE_CHECK.json',dict(status='PASS',request=request,computed_shear_Pa=shear,computed_velocity_m_s=result['velocity_m_s'],computed_gradient_s_inv=result['gradient_s_inv'],missing_bubble_response=missing,flow_classes=FLOW_CLASSES,native_steps=0,scope='New unified API only; no extra fracture case'))
print('PUBLIC_INTERFACE_PASS')
