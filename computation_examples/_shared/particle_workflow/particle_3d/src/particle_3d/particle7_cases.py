"""Synthetic infrastructure fixtures, explicitly not production suspension physics."""
from functools import lru_cache
from pathlib import Path
import numpy as np
import pyvista as pv
from .inlet_flux import InletFluxSampler,frozen_boundary_flux
from .injection_population import PopulationSource,InjectionScheduler,LinearProfile
from .injection_admission import FiniteSizeAdmission,PlugPathAdmission
from .particle7_lifecycle import PopulationLifecycle
from .particle7_bridge import DynamicParticleBridge
from .lammps_neighbors import ValidationNeighborPolicy
from .validation_boundary import ValidationBoundaryClassifier
from .wall_geometry import WallGeometry
from .audit import read_frozen
from .field import FrozenFEMField

REPO=Path(__file__).resolve().parents[3]
SONOVUE=Path('/home/lzy/projects/sonovue_size_distribution_v0')


def rectangle(a,b,c,d): return np.array([[a,b,c],[a,c,d]],float)
def surface(triangles):
    return pv.PolyData(np.asarray(triangles).reshape(-1,3),np.column_stack([np.full(len(triangles),3),np.arange(len(triangles)*3).reshape(-1,3)]).ravel())


def channel(width=200e-6,length=1e-8,Q=1e-12):
    w=width/2; velocity=Q/width**2
    inlet=rectangle([-w,-w,0],[w,-w,0],[w,w,0],[-w,w,0])
    outlet=inlet.copy(); outlet[:,:,2]=length
    wall=[]
    for a,b in zip(inlet[[0,0,0,1],[0,1,2,2]],inlet[[0,0,1,1],[1,2,2,0]]):
        # Orient side triangles outward (interior left of boundary edge).
        c=b+[0,0,length]; d=a+[0,0,length]; wall.extend(rectangle(a,b,c,d))
    return InletFluxSampler(inlet,np.full((2,3),velocity)),WallGeometry(wall),ValidationBoundaryClassifier({'OUTLET_01':surface(outlet)}),np.array([0.,0.,velocity])


def synthetic_engine(*,bridge=True,seed=2026092107,width=200e-6,length=1e-8,guard=512):
    sampler,wall,classifier,v=channel(width,length)
    source=PopulationSource(SONOVUE,seed); scheduler=InjectionScheduler(source,LinearProfile([0],[sampler.Q_m3_s]))
    admission=PlugPathAdmission(sampler,source,width=width,velocity=v,guard=guard)
    b=DynamicParticleBridge([],ValidationNeighborPolicy(20e-6,.5e-6,'P7_INFRASTRUCTURE_VALIDATION_ONLY')) if bridge else None
    return PopulationLifecycle(scheduler,admission,classifier,bridge=b)


@lru_cache(maxsize=1)
def real_inlet():
    fem=REPO/'formal_3D_flow_solver/FEM_SimVascular'
    provenance,mesh,flow,boundaries=read_frozen(fem)
    audit,samplers=frozen_boundary_flux(mesh,flow,boundaries)
    return provenance,mesh,FrozenFEMField.from_grids(mesh,flow),boundaries,audit,samplers,WallGeometry.from_frozen(fem)
