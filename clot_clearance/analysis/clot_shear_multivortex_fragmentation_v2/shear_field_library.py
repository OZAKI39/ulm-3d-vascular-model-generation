"""Two bounded velocity templates, not solved acoustic streaming.

The nominal shear is beta * outer-speed and bounds the stationary-surface
driving traction. Actual slip traction also contains the separately audited drag.
"""
from common import *

DEFAULT = dict(nominal_shear_Pa=1000., angle_rad=0., center_x_m=0.,
               center_y_m=0., width_m=125e-6, separation_m=300e-6,
               drive_duration_s=OFF)
BOUNDS = dict(nominal_shear_Pa=[0.,6000.], angle_rad=[-np.pi,np.pi],
              center_x_m=[-100e-6,100e-6], center_y_m=[-50e-6,50e-6],
              width_m=[100e-6,200e-6], separation_m=[200e-6,400e-6],
              drive_duration_s=[.0014,.0028])

def plateau(r, inner, outer):
    """C2 plateau taper, exactly zero beyond the fixed laboratory region."""
    q=np.clip((np.abs(r)-inner)/(outer-inner),0.,1.)
    return 1.-q*q*q*(10.+q*(-15.+6.*q))

class ShearField:
    def __init__(self, mode, theta=None, surface_z=375e-6):
        assert mode in CASES
        self.mode=mode; self.theta=DEFAULT.copy()
        if theta is not None:
            assert set(theta)<=set(DEFAULT)
            self.theta.update(theta)
        for key,value in self.theta.items():
            assert np.isfinite(value) and BOUNDS[key][0]<=value<=BOUNDS[key][1],(key,value)
        self.z0=surface_z; self.beta=BETA
        self.U=self.theta['nominal_shear_Pa']/self.beta
        self.segment=(0.,END,3200.,3200.)
        phi=self.theta['angle_rad']
        self.e1=np.array([np.cos(phi),np.sin(phi),0.])
        self.e2=np.array([-np.sin(phi),np.cos(phi),0.])
        self.center=np.array([self.theta['center_x_m'],self.theta['center_y_m'],surface_z])
        self.stop=self.theta['drive_duration_s']
        self.rise=self.stop*5/7; self.hold=self.stop*6/7

    def envelope(self,t):
        t=np.asarray(t); rise=.5*(1.-np.cos(np.pi*np.clip(t/self.rise,0.,1.)))
        down=.5*(1.+np.cos(np.pi*np.clip((t-self.hold)/(self.stop-self.hold),0.,1.)))
        return np.where(t<0.,0.,np.where(t<self.rise,rise,np.where(t<self.hold,1.,down)))

    def spatial(self,x,n):
        """Current positions and oriented current facet unit normals."""
        q=np.asarray(x)-self.center; a=q@self.e1; b=q@self.e2
        window=plateau(a,400e-6,700e-6)*plateau(b,250e-6,450e-6)*plateau(q[:,2],100e-6,300e-6)
        if self.mode=='uniform':
            raw=np.broadcast_to(self.e1,x.shape).copy()
            flow=raw-n*np.sum(raw*n,axis=1)[:,None]
        elif self.mode=='counter_rotating':
            sigma=self.theta['width_m']; separation=self.theta['separation_m']
            a1=(a+separation/2)/sigma; a2=(a-separation/2)/sigma; bb=b/sigma
            g1=np.exp(-.5*(a1*a1+bb*bb)); g2=np.exp(-.5*(a2*a2+bb*bb))
            # psi = sigma*exp(1/2)/2 * (g1-g2). Each Gaussian gradient
            # contributes at most 1/2, hence |grad psi| <= 1 globally.
            ga=np.exp(.5)/2*(-a1*g1+a2*g2)
            gb=np.exp(.5)/2*(-bb*g1+bb*g2)
            gradient=ga[:,None]*self.e1+gb[:,None]*self.e2
            flow=np.cross(n,gradient)
        else:
            # Actual 650 x 350 um top: nearest 250 um count is 3 x 1.
            # Freeze a square 250 um periodic lattice, centered on the specimen.
            # Three cell centers lie inside the actual top; nine cells intersect it, with six clipped strips and two clipped end cells.
            scale=250e-6;k=np.pi/scale
            ga=np.cos(k*(a+375e-6))*np.sin(k*(b+125e-6))
            gb=np.sin(k*(a+375e-6))*np.cos(k*(b+125e-6))
            gradient=ga[:,None]*self.e1+gb[:,None]*self.e2
            flow=np.cross(n,gradient)
        return flow*window[:,None]

    def __call__(self,x,n,t):
        return self.amplitude(t)/BETA*float(self.envelope(t))*self.spatial(x,n)

    def amplitude(self,t):
        a,b,u,v=self.segment
        return u+(v-u)*np.clip((t-a)/max(b-a,1e-30),0.,1.)

    def applicability(self):
        return dict(label=LABEL,nominal_stationary_shear_cap_Pa=self.theta['nominal_shear_Pa'],
                    outer_speed_bound_m_s=self.U,mu_Pa_s=MU,h_eff_m=H_EFF,
                    Re_h_bound=RHO*self.U*H_EFF/MU,Re_width_bound=RHO*self.U*self.theta['width_m']/MU,
                    h_over_width=H_EFF/self.theta['width_m'],
                    acoustic_realizability='UNESTABLISHED',
                    meaning='Hypothetical thin-layer slip law; no first/second-order bubble flow solution. Re bounds disallow interpreting the peak as a low-Re Stokes prediction.')
