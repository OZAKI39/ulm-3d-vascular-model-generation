"""Resource estimates only. No material extrapolation or solver allocation."""
import time
import numpy as np
from scipy.fft import next_fast_len
from .common import *

def main():
 start=time.perf_counter();cfg=read(CONFIG)['D']['pdms_preflight'];f=cfg['frequency_Hz']
 reference=read(ROOT/'results/straight_tube_exp1/f10m/materials/REFERENCE.json')
 # Conservative reference speed is a planning bound, NOT admitted 2.5MHz material data.
 speed=reference['zero_frequency_sound_speed_m_s'];rows=[]
 for ppw in cfg['PPW']:
  h=speed/f/ppw;pad=cfg['PML_cells'];physical=np.array([.0325,.009,.009])+2*cfg['water_envelope_m']
  # Old kWaveArray BLI support allowance retained separately from PML.
  margin=36;required=np.ceil(physical/h).astype(int)+2*margin
  full=np.array([next_fast_len(int(n+2*pad)) for n in required]);shape=full-2*pad
  dt=1/(np.ceil(1500/(cfg['CFL']*h*f))*f);duration=2*np.linalg.norm(shape*h)/speed+16/f;nt=int(np.ceil(duration/dt))
  vox=int(np.prod(full));ratio=vox/(378*192*192)
  memory=15100/1024*ratio*1.2;seconds=373.70579146797536/3*ratio*(nt/3504)*1.5
  rows.append(dict(PPW=ppw,spacing_m=h,shape=shape,expanded_shape=full,PML_cells=pad,BLI_margin_cells=margin,
     axes_bounds_m=np.column_stack([-shape*h/2,shape*h/2]),CFL=1500*dt/h,dt_s=dt,Nt=nt,duration_s=duration,
     reference_longitudinal_wavelength_m=speed/f,water_wavelength_m=1500/f,
     estimated_peak_GiB=memory,estimated_single_source_wall_s=seconds,
     vram_admitted=bool(memory<=cfg['gpu_budget_GiB']),time_admitted=bool(seconds<=120)))
 reasons=['2.5MHz lies outside the saved 3-7MHz measured reference band; the original material module rejects extrapolation',
          'End-seal material/thickness and exterior coupling are not measured',
          'Conditional grid estimates exceed 20GiB safe budget and 120s diagnostic allowance',
          'Equivalent-fluid PDMS omits shear waves; no experimental calibration']
 write(BASE/'diagnostics/PDMS_PREFLIGHT.json',dict(status='PDMS_Y_FIELD_NOT_ADMITTED',PDMS_FIELD_ADMITTED='FAIL',PDMS_FIELD_NOT_RUN=True,
   actual_sources_run=0,material_reference=reference,frequency_Hz=f,planning_speed_m_s=speed,
   speed_scope='Literature-derived zero-frequency limit for resource planning ONLY; not a 2.5MHz material input',
   material_regions='PDMS bulk; water-filled central 1.6mm channel and two 1.5mm vertical ports. Both X ends geometrically sealed. Seal acoustic material unresolved.',
   exterior='Conditional 1mm water envelope; not measured coupling',
   source_geometry=read(OLD/'CONFIG_SNAPSHOT.json')['geometry'],candidates=rows,reasons=reasons,
   literature=dict(title='Acoustic Characterization of Polydimethylsiloxane for Microscale Acoustofluidics',authors='Guangyao Xu et al.',year=2020,
   DOI='10.1103/PhysRevApplied.13.054069',url='https://www.bruus-lab.dk/TMF/publications/Bruus/Bruus_153.pdf',
   frequency_range_Hz=[3e6,7e6],data='Tables II and III, Sylgard 184 10:1, cured 85 C, measured 25 C; archived F10M reference. Published Kramers-Kronig fits may be extrapolated as a new assumption, but no such unvalidated input admitted here.'),
   estimate_method='Measured F10L total 373.705791s /3; 15100MiB peak, expanded grid378x192x192,Nt3504. Voxel linear scaling; time also Nt; memory1.2 and time1.5 heterogeneity factors. Not measured F10S runtime or guaranteed upper bounds.',
   wall_s=time.perf_counter()-start,identity=identity(),training_ready=False))
 print([(r['PPW'],r['estimated_peak_GiB'],r['estimated_single_source_wall_s']) for r in rows])
if __name__=='__main__':main()
