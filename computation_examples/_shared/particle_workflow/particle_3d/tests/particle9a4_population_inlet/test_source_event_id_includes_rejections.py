from particle_3d.population_inlet_p9a4 import generate

def test_every_proposal_has_an_id(simple_source):
 ledger=generate(simple_source,128)
 assert [r['source_event_id'] for r in ledger.rows]==list(range(1,129))
 assert any(r['status']=='REJECTED_SOURCE_EVENT' for r in ledger.rows)
 assert len(ledger.rows)==128
