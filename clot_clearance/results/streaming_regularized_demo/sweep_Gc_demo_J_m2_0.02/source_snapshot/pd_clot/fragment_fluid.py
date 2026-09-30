"""Prescribed velocity provider and replaceable verification transport operator."""
import numpy as np
from .streaming import pipe_field, AnalyticTestStreaming


class FragmentFluid:
    def __init__(self, config):
        self.c=config
        self.traction_provider=AnalyticTestStreaming(config['streaming'],config['pipe'],config['simulation']['representative_frequency_Hz'])

    def velocity(self, positions, time_s):
        c=self.c;s=c['streaming']
        velocity,_,_=pipe_field(positions,c['pipe'])
        velocity*=c['transport']['background_velocity_multiplier']
        d=positions-np.asarray(s['bubble_center_m']);L=s['streaming_length_scale_m']
        envelope=np.exp(-.5*np.sum(d*d,axis=1)/L**2)
        # Curl of scalar streamfunction in x-z plane: divergence-free vortex.
        speed=s['streaming_velocity_scale_m_s']
        factor=s['mean_fraction']+s['oscillatory_amplitude']*np.sin(2*np.pi*c['simulation']['representative_frequency_Hz']*time_s+s['phase_rad'])
        vortex=np.column_stack((-d[:,2],np.zeros(len(d)),d[:,0]))/L
        velocity+=speed*factor*envelope[:,None]*vortex
        return velocity

    def surface_force(self, cloud, positions, integrity, exposed, attached, time_s):
        weights=cloud.weight*integrity
        eta=positions[cloud.pairs[:,1]]-positions[cloud.pairs[:,0]]
        length=np.linalg.norm(eta,axis=1)
        direction=eta/np.maximum(length[:,None],np.finfo(float).tiny)
        moment=np.zeros_like(positions)
        np.add.at(moment,cloud.pairs[:,0],weights[:,None]*direction)
        np.add.at(moment,cloud.pairs[:,1],-weights[:,None]*direction)
        norm=np.linalg.norm(moment,axis=1)
        normals=np.zeros_like(positions)
        valid=norm>self.c['surface']['normal_moment_tolerance']
        normals[valid]=-moment[valid]/norm[valid,None]
        # Free fragments use relaxation only, avoiding traction+drag double loading.
        load=exposed & attached & valid & ~cloud.fixed
        traction=np.zeros_like(positions)
        traction[load]=self.traction_provider.traction(positions[load],normals[load],time_s)
        area=cloud.spacing**2*self.c['surface']['area_factor']
        return traction*area,normals,traction


def relaxation_half_step(velocity, fluid_velocity, mobile, dt_half, tau):
    if tau<=0:raise ValueError('tau_h_s must be positive')
    result=velocity.copy()
    result[mobile]=fluid_velocity[mobile]+(velocity[mobile]-fluid_velocity[mobile])*np.exp(-dt_half/tau)
    return result
