import ast,json
from pathlib import Path
import numpy as np
from particle_3d.audit import sha256
from particle_3d.wall_geometry import WallGeometry


def test_particle2_dependency_and_distribution_locked(p3_repo,p3_population):
    p2=json.loads((p3_repo/'particle_3d/reports/particle2/PARTICLE2_VALIDATION.json').read_text())
    assert p2['git_commit']=='c704e08d39f6134da300fc91cc93c8561314a17b'
    assert p2['manual_visual_review']=='PASS'
    for group in ['source_sha256','data_sha256','figure_sha256']:
        for path,digest in p2[group].items():assert sha256(p3_repo/path)==digest,path
    meta=json.loads((p3_repo/'particle_3d/reports/particle3/data/07_selection.json').read_text())
    assert meta['population_array_sha256']==p3_population.metadata['sample_structured_array_sha256']
    assert len(set(meta['indices']))==128


def test_wall_input_integrity_and_only_wall_is_solid(p3_repo,real_wall):
    assert len(real_wall.triangles)==45221
    assert real_wall.provenance['sha256']=='06c5d1d1975cf9470096d9b7e560e2173c9733c847d5bc364aab7355ef12a11b'
    # Loader receives the complete manifest; its actual dataset is solely WALL.
    assert 'WALL.vtp' in real_wall.provenance['path']
    assert real_wall.provenance['solid_boundaries']==['WALL']
    assert not real_wall.triangles.flags.writeable


def test_original_triangle_winding_and_inward_owner(p3_repo,real_wall):
    from particle_3d.field import FrozenFEMField
    from particle_3d.particle3_cases import normal_audit
    stored=json.loads((p3_repo/'particle_3d/reports/particle3/data/01_wall_normals.json').read_text())
    field=FrozenFEMField.from_frozen(p3_repo/'formal_3D_flow_solver/FEM_SimVascular')
    current=normal_audit(real_wall,field,[row['triangle_id'] for row in stored])
    for before,row in zip(stored,current):
        assert row['owner_inward_distance_m']>0
        assert row['owning_tetra']==before['owning_tetra']
        np.testing.assert_array_equal(row['vertices_m'],real_wall.triangles[row['triangle_id']])
        np.testing.assert_allclose(row['normal_in'],-row['normal_out'],atol=0,rtol=0)


def test_particle3_scope_has_no_external_physics_solver(p3_repo):
    names=['particle_shapes','convex_triangle','wall_geometry','wall_gap','wall_contact','physical_time_refinement','rbc_capillary_surrogate','particle3_motion']
    for name in names:
        tree=ast.parse((p3_repo/f'particle_3d/src/particle_3d/{name}.py').read_text())
        imports=[node.module or '' for node in ast.walk(tree) if isinstance(node,ast.ImportFrom)]
        imports += [a.name for node in ast.walk(tree) if isinstance(node,ast.Import) for a in node.names]
        assert not any(any(x in mod.lower() for x in ['lammps','adhesion','lubrication','particle4','penalty']) for mod in imports)
    scope=json.loads((p3_repo/'particle_3d/reports/particle3/data/00_particle3_scope_and_wall_contract.json').read_text())
    assert scope['no_cfd_executed'] and not scope['particle4_started']
    assert not scope['production_particle_timestep_frozen']
