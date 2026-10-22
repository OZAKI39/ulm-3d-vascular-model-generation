"""Pre-acquisition computational-guard choice from completed diagnostic timing."""

def choose_natural_guard(scaling,extended,budget):
    count=budget['natural_initial_scheduled'];limit=budget['natural_batch_wall_budget_hours']*3600
    base=count/scaling['parallel_tracks_hour']*3600
    # Fixed original 720/2200 stop frequency, not new-sample outcomes. Each
    # continuation factor records total wall cost for those 720 original IDs.
    estimates={1:base}
    for trial in extended['factors']:
        if not trial['all_original_prefixes_exact']:raise ValueError('Cannot select from an invalid continuation')
        if not trial['performance'].get('comparable_full_fresh_batch',True):
            raise ValueError('Resume lookup time cannot predict full natural integration cost')
        estimates[trial['guard_factor']]=base+count/2200*trial['performance']['wall_seconds']
    rows=[dict(guard_factor=k,predicted_seconds=v,predicted_seconds_with_margin=1.15*v,
               fits_predeclared_wall_budget=1.15*v<=limit) for k,v in sorted(estimates.items())]
    fitting=[r['guard_factor'] for r in rows if r['fits_predeclared_wall_budget']]
    if not fitting:raise ValueError('Even baseline forecast exceeds predeclared natural wall budget')
    return dict(chosen_guard_factor=max(fitting),natural_scheduled=count,budget_seconds=limit,estimates=rows,
        timing_safety_factor=1.15,original_stop_fraction=720/2200,
        rule='Largest tested guard predicted to fit fixed 3-hour natural budget with 15% timing margin. Baseline natural cost plus scaled original-stop continuation cost.',
        prediction_is_not_a_runtime_guarantee=True,selection_before_new_birth_ledger=True,
        no_new_sample_outcome_or_outlet_conditioning=True,no_change_to_dt_physics_admission_or_sampling=True)
