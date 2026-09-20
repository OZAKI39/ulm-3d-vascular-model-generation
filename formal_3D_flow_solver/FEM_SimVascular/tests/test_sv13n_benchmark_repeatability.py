from sv13n_support import *
def test_formal_benchmarks_explicitly_deferred():
 for name in ("cpu_version_benchmark","gpu_benchmark"):
  d=read(name);assert d["status"]=="NOT_RUN" and d["classification"]=="DEFERRED"
  assert "median_s" not in d and "speedup" not in d
