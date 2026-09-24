def test_ids(scheduler):
 rows=list(scheduler.through(.12)); assert [x['particle_id'] for x in rows]==list(range(1,len(rows)+1))
 assert scheduler.next_particle_id==len(rows)+1
