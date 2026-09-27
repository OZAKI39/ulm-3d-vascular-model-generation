from particle_3d.population_inlet_p9a4 import generate
from particle_3d.continuous_infusion import canonical_bytes,PopulationLedger
import pytest

def test_workers_1_3_6_bitwise_same(real_source):
 ledgers=[generate(real_source,192,workers=w) for w in [1,3,6]]
 assert len({canonical_bytes(x.rows) for x in ledgers})==1
 assert len({canonical_bytes(x.birth_events()) for x in ledgers})==1

def test_reversed_merge_and_duplicate_rejection(simple_source):
 rows=[simple_source.proposal(i) for i in range(1,65)]
 a=PopulationLedger(simple_source.identity);a.merge(rows[::-1]);b=PopulationLedger(simple_source.identity);b.merge(rows)
 assert canonical_bytes(a.rows)==canonical_bytes(b.rows)
 with pytest.raises(ValueError):PopulationLedger(simple_source.identity).merge(rows+[rows[0]])
