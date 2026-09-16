from pathlib import Path
import tempfile, unittest, sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from prepare_lammps_particles import load_population,write_case

class ConverterTests(unittest.TestCase):
    def test_si_and_sphere_data_roundtrip(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'input.csv';p.write_text('bubble_id,diameter_um\n0,0.8\n1,5.2\n')
            ids,d=load_population(p)
            np.testing.assert_array_equal(ids,[1,2]);np.testing.assert_array_equal(d,np.array([0.8,5.2])*1e-6)
            spec={'N':2,'box_m':[100e-6]*3,'placement':'two','surface_gap_m':10e-6,'velocity_m_s':[1e-3,-2e-3,3e-3],'omega_rad_s':[0,0,0],'technical_density_kg_m3':1000,'contact':False,'mpi_ranks':1,'dt_s':1e-6,'steps':2}
            out=Path(t)/'out';write_case(p,spec,out)
            lines=(out/'particles.data').read_text().split('Atoms # sphere\n\n')[1].split('\nVelocities')[0]
            a=np.loadtxt(lines.splitlines());np.testing.assert_array_equal(a[:,2],d)
            np.testing.assert_array_equal(a[:,3],[1000,1000]);self.assertGreater(a[1,4]-a[0,4],(d[0]+d[1])/2)
            self.assertIn('units si',(out/'in.lammps').read_text())
    def test_reject_duplicate_or_out_of_order_id(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'p.csv';p.write_text('bubble_id,diameter_um\n1,2\n0,1\n')
            with self.assertRaises(ValueError):load_population(p)
    def test_reject_nonfinite_diameter(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'p.csv';p.write_text('bubble_id,diameter_um\n0,nan\n')
            with self.assertRaises(ValueError):load_population(p)

if __name__=='__main__':unittest.main(verbosity=2)
