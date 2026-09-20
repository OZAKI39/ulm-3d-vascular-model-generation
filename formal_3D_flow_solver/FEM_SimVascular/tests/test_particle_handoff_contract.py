import re
from fem_sync_support import ROOT, read


def test_new_chat_can_find_all_frozen_inputs_and_units():
    text=(ROOT/'PARTICLE_HANDOFF.md').read_text()
    records=[read('frozen_reference/flow/flow_field_manifest.json'),read('frozen_reference/mesh_manifest.json')]
    records+=list(read('frozen_reference/boundary_manifest.json')['boundaries'].values())
    for r in records:
        assert r['path'] in text and r['sha256'] in text and (ROOT/r['path']).is_file()
    for term in ['FEM DEVELOPMENT','FROZEN','SV1.3Q','ILU_REBUILD_POLICY_WINNER_FOUND','Frozen Field Contract',
                 'Velocity','Pressure','NOT STORED IN FROZEN VTU','m/s','Pa','kg/m³','m³/s','新聊天第一件事','Particle-0']:
        assert term in text


def test_deferred_validation_and_future_contract_are_not_misrepresented():
    text=(ROOT/'PARTICLE_HANDOFF.md').read_text()
    assert 'NOT PERFORMED BY USER DECISION' in text and 'DEFERRED BY USER DECISION' in text
    for name in ['position[3]','velocity[3]','pressure','velocity_gradient[3,3]','vorticity[3]','strain_rate[3,3]','inside_lumen','tetra_id']:
        assert name in text
    assert 'RBC / MB → FEM 不允许' in text and 'red_blood_cell_transport.py' in text
    roadmap=(ROOT/'PARTICLE_RESEARCH_ROADMAP.md').read_text()
    assert all(f'Particle-{i}' in roadmap for i in range(9))


def test_handoff_document_local_links_resolve():
    for name in ['PARTICLE_HANDOFF.md','FROZEN_FEM_BASELINE.md','README.md','SYNC_REPORT_FEM_PARTICLE_HANDOFF.md']:
        text=(ROOT/name).read_text()
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if not target.startswith('https://'):
                assert (ROOT/target.split('#')[0]).exists(),(name,target)
