from sv12_support import accepted,artifact
def test_pressure_uses_area_weighted_ports_and_sections():
    final=accepted();data=artifact('pressure_sections')
    assert data['port_area_averages_pa']==final['area_average_pressure_pa']
    assert data['raw_pressure_offset_applied_pa']==0
    assert all(row['area_m2']>0 and row['triangles']>0 for row in data['sections'])
