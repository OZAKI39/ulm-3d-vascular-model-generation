from copy import deepcopy
import json
import pytest
from particle_3d.formal_cohort_p9a5 import *

def test_formal_extensions_are_contiguous(population,core):
 old=core
 for n in [750,1000,1250,1500,1750,2000,2500,3000,3500,4000,4500,5000]:
  new=cohort(population,n,core['master_seed']);require_prefix(old,new)
  assert [e['particle_id'] for e in new['events'][old['count']:]]==list(range(old['count']+1,n+1))
  old=new
 assert next_size(500,False)==750 and next_size(2000,False)==2500
 assert next_size(500,True)==500 and next_size(5000,False)==5000
