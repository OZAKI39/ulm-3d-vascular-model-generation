"""Only bouncer ownership varies; both fluid sides retain bounce-back."""
def bind_bouncers(mir, coordinator, membrane, outer, inner, policy):
    if policy not in ('shared','independent'):
        raise ValueError('UNKNOWN_BOUNCER_OWNERSHIP')
    handles=[]
    for name,pv in (('outer',outer),('inner',inner)):
        if policy=='independent' or not handles:
            bounce=mir.Bouncers.Mesh('membrane_bounce_'+name,'bounce_back')
            coordinator.registerBouncer(bounce);handles.append(bounce)
        coordinator.setBouncer(handles[-1],membrane,pv)
    return handles
