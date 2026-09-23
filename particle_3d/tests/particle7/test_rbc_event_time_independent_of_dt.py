from particle_3d.injection_population import InjectionScheduler
def test_partition(scheduler):
 saved=scheduler.state(); a=list(scheduler.through(.03)); b=InjectionScheduler.restore(saved); rows=[]
 for t in [.001,.002,.013,.023,.03]: rows.extend(b.through(t))
 assert a==rows
