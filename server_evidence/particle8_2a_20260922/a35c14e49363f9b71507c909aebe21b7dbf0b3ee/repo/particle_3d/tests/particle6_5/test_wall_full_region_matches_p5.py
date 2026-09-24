def test_full_matches_historical(scans):
 rows=[r for r in scans[0] if 2e-9<r['h_geom_m'] and r['chi']<=.01]
 assert rows and all(r['coefficient_kg_s']==r['old_p5_default_zeta_kg_s'] for r in rows)
