"""Preproduction storage estimate includes measured adaptive substeps and rejected trials."""
import math

def make_plan(benchmark,first_worker_result,free_bytes):
 ages={r['particle_id']:r['trajectory_age_s'] for r in first_worker_result['rows']}
 # bytes / physical age accounts for all accepted substeps and rejected-trial logs.
 # bytes / accepted step times nominal steps would underestimate adaptive output.
 costs=benchmark['benchmark_cost_rows'];dt=benchmark['dt_s']
 worst_rate=max(r['bytes']/max(ages[r['particle_id']],dt) for r in costs)
 worst_bytes=math.ceil(worst_rate*12.*5000*1.25)
 safe_disk=free_bytes*.8
 maximum=5000 if worst_bytes<=safe_disk else min(5000,math.floor(safe_disk/(worst_bytes/5000)/250)*250)
 if maximum<500:raise ValueError('Conservative disk allowance cannot support CORE500; retain data and review resources before production')
 return dict(dt_s=dt,preferred_N_max=5000,resource_limited_N_max=maximum,hard_limit_reason=None if maximum==5000 else 'MEASURED_WORST_BYTES_PER_PHYSICAL_AGE_WITH_ADAPTIVE_REFINEMENT',current_free_disk_bytes=free_bytes,worst_case_5000_storage_bytes=worst_bytes,worst_case_5000_runtime_seconds=benchmark['worst_case_5000_runtime_seconds'],cost_method='Maximum measured total bytes / physical trajectory age x12 seconds x5000 x1.25; includes refined accepted steps and rejected-trial logs. CPU uses unchanged provider-call guard. Empirical conservative envelope, not a strict analytical bound.',maximum_measured_bytes_per_physical_second=worst_rate,declared_before_production=True)
