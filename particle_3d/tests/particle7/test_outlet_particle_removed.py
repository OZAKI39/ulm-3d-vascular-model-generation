def test_no_dead_active(engine):
 engine.step_to(.003); ids={e['particle_id'] for e in engine.exits}
 assert ids and not ids.intersection(engine.active)
 assert not ids.intersection(p.particle_id for p in engine.bridge.read())
 assert not ids.intersection(i for pair in engine.bridge.raw_pairs for i in pair)
