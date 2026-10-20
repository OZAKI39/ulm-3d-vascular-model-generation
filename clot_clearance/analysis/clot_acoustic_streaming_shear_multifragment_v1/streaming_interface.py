"""Unified provenance-aware SI entrypoint for future acoustic-fluid coupling."""
from acoustic_streaming_shear import *
FLOW_CLASSES={
 'FIRST_ORDER_OSCILLATORY_FLOW':'Acoustic-period velocity; pressure phasors alone do not supply the near-wall viscous solution.',
 'SECOND_ORDER_MEAN_STREAMING':'Cycle-averaged second-order flow with consistent oscillatory boundary conditions; not yet solved here.',
 'PRESCRIBED_STOKES_SHEAR_TEMPLATE':'Explicitly imposed non-acoustic viscous boundary forcing, computed in this round.',
}
def evaluate(request):
 label=request['model_label']
 if label=='BUBBLE_DRIVEN_STREAMING_MODEL':return bubble_model_input_gate(request)
 if label!='PRESCRIBED_STOKES_SHEAR_TEMPLATE':raise ValueError('Unsupported physical model; inverse requirements belong to the mechanical load adapter')
 # The frozen validated model has a defined parameter domain. Do not silently
 # infer viscosity, gap, mean boundary velocity or clot registration from sound.
 expected=dict(rho_kg_m3=FLUID_RHO,mu_Pa_s=MU,gap_m=GAP,upper_speed_m_s=SPEED,wavelength_m=WAVELENGTH)
 for key,value in expected.items():
  if request.get(key)!=value:raise ValueError(f'{key} must explicitly match this frozen template ({value})')
 time=float(request['time_s']);assert 0<=time<=END
 coordinates=np.asarray(request['local_surface_coordinates_m']);normals=np.asarray(request['solid_to_fluid_normals']);assert coordinates.shape==normals.shape and coordinates.shape[-1]==3
 assert np.all((coordinates[:,2]>=0)&(coordinates[:,2]<=GAP)), 'Coordinates must be in the explicitly defined local fluid gap'
 assert np.max(abs(np.linalg.norm(normals,axis=-1)-1))<1e-10
 flow=StokesMode();fields=flow.field(coordinates[:,1],coordinates[:,2],time);shear,normal=traction(fields['gradient_s_inv'],fields['pressure_increment_Pa'],normals)
 return dict(status='COMPUTED_TEMPLATE',model_label=label,evidence='PHYSICS_APPROXIMATION',**fields,tangential_traction_Pa=shear,normal_traction_Pa=normal,time_s=time,model_spec=spec(),acoustic_inputs_used=False,bubble_inputs_used=False,surface_scope='Stress projection in the prescribed local gap. Arbitrary normals do not create a new curved no-slip boundary solution.')
