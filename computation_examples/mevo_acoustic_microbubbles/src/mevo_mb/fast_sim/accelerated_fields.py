"""Same Warp kernels and force refresh, one packed D2H field read instead of four.

All launches and device copy use Warp's same current stream. The host numpy()
read synchronizes that stream. No external owner or zero-copy lifetime is needed.
Walls retain the unchanged two independent checks and their original transfers.
"""
import numpy as np
import warp as wp
from .accelerated_water import WaterAccelerated

@wp.kernel
def pack_fields(pressure:wp.array(dtype=wp.vec2d),primary:wp.array(dtype=wp.vec3d),total:wp.array(dtype=wp.vec3d),error:wp.array(dtype=wp.int32),out:wp.array2d(dtype=wp.float64)):
    i=wp.tid()
    out[i,0]=pressure[i][0];out[i,1]=pressure[i][1]
    for k in range(3):
        out[i,2+k]=primary[i][k];out[i,5+k]=total[i][k]
    out[i,8]=wp.float64(error[0])

class WaterPacked(WaterAccelerated):
    def __init__(self,*a,**kw):
        super().__init__(*a,**kw)
        self.primary_device=wp.empty_like(self.sim.force)
        self.field_buffer=wp.zeros((self.sim.total,9),dtype=wp.float64,device=self.sim.device)
        wp.load_module(module=pack_fields.module,device=self.sim.device)
    def fields(self,x):
        coord=(x-self.basis.origin)/self.basis.spacing
        if np.any(coord<1) or np.any(coord>=np.array(self.basis.shape)-2):
            raise ValueError('Outside stored cubic field interior; no clipping')
        s=self.sim;s.x.assign(x);s.fields();wp.copy(self.primary_device,s.force);s.pairs()
        wp.launch(pack_fields,s.total,inputs=[s.pressure,self.primary_device,s.force,s.error,self.field_buffer],device=s.device)
        a=self.field_buffer.numpy()
        if np.any(a[:,8]):raise ValueError(f'Original acoustic/APECSS fail-stop bits={int(a[:,8].max())}')
        if s.flow_enabled:raise ValueError('Nonzero-flow configuration')
        primary=a[:,2:5].copy();total=a[:,5:8].copy()
        return dict(pressure_Pa=a[:,0]+1j*a[:,1],primary_force_N=primary,secondary_force_N=total-primary,total_force_N=total)
