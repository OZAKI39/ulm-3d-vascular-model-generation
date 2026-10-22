from pathlib import Path
import json
import numpy as np
import pytest
from network_1d0d.balance_design import load_frozen_design
from network_1d0d.parameterized_fem_handoff import export_frozen_design,measure_extension_geometry,transfer_prediction
from network_1d0d.fem_h0_case import serialize_fixed_pressure_xml,audit_xml_change

ROOT=Path(__file__).resolve().parents[2]
FROZEN=ROOT/'reports/parameterized_0d_v1/best_feasible_balance/frozen_balance_design.yaml'
SOURCE=ROOT.parent/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1'
PORTS=ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'


@pytest.fixture(scope='module')
def handoff():
    return export_frozen_design(FROZEN,SOURCE,PORTS)


def test_parameterized_handoff_preserves_pressure_difference(handoff):
    np.testing.assert_allclose(handoff['pairwise_raw_pressure_differences_Pa'],handoff['pairwise_shifted_pressure_differences_Pa'],rtol=1e-12,atol=1e-12)
    assert handoff['source_parameter_status']=='FROZEN_DESIGN'


def test_parameterized_handoff_uses_configured_mu(handoff):
    ports=json.loads(PORTS.read_text())['ports']
    double=measure_extension_geometry(SOURCE,ports,mu_pa_s=2*handoff['mu_pa_s'])
    np.testing.assert_allclose([r['R_extension_Pa_s_m3'] for r in double],
                              [2*r['R_extension_Pa_s_m3'] for r in handoff['extension_measurements']],rtol=1e-13)


def test_analytical_H0_handoff_regression(handoff):
    f=load_frozen_design(FROZEN)
    baseline=transfer_prediction(f['baseline']['prediction'],handoff['extension_measurements'])
    gold=json.loads((ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data/roi_fixed_pressure_bc_H0.json').read_text())
    for field in ('R_extension_Pa_s_m3','pressure_cap_raw_Pa','pressure_cap_shifted_Pa'):
        np.testing.assert_allclose([r[field] for r in baseline['ports']],[r[field] for r in gold['ports']],rtol=2e-12,atol=1e-9)


def test_only_three_outlet_pressure_values_change(handoff):
    old=(SOURCE/'run/solver.xml').read_bytes()
    values=[r['pressure_cap_shifted_Pa'] for r in handoff['ports']]
    new=serialize_fixed_pressure_xml(old,values)
    changes=audit_xml_change(old,new,values)
    assert len(changes)==3 and all('OUTLET_' in r['field'] for r in changes)
    bad=new.replace(b'<Density>1056.0</Density>',b'<Density>1000.0</Density>')
    with pytest.raises(ValueError,match='Unexpected'):audit_xml_change(old,bad,values)


def test_no_cfd_field_or_feedback_read(monkeypatch):
    import pyvista as pv
    reader=pv.read;reads=[]
    def geometry_only(path,*args,**kwargs):
        assert Path(path)==SOURCE/'SV_MESH/mesh-complete.exterior.vtp'
        reads.append(path);return reader(path,*args,**kwargs)
    monkeypatch.setattr(pv,'read',geometry_only)
    h=export_frozen_design(FROZEN,SOURCE,PORTS)
    assert len(reads)==1 and h['CFD_feedback_used'] is False
