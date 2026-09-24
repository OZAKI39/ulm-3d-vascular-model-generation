from particle_3d.injection_population import PopulationSource
def test_separate(source):
 clone=PopulationSource.restore(source.state()); source.orientation(3000); source.rng['ADMISSION_RETRY'].random(2000)
 assert [source.next_mb() for _ in range(100)]==[clone.next_mb() for _ in range(100)]
 assert source.next_rbc()['geometry']==clone.next_rbc()['geometry']
