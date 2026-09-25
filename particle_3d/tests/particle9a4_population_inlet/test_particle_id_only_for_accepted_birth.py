from particle_3d.population_inlet_p9a4 import generate

def test_accepted_ids_and_rejected_null_ids(simple_source):
 ledger=generate(simple_source,128)
 ids=[r['particle_id'] for r in ledger.rows if r['status']=='ACCEPTED_BIRTH']
 assert ids==list(range(1,len(ids)+1))
 assert all(r['particle_id'] is None for r in ledger.rows if r['status']=='REJECTED_SOURCE_EVENT')
