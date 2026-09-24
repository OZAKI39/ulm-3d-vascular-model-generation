 format = 'asciispatial'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'tau1_reference_scaled_smoke'
 basename = 'tracking/state_000001/tau1_reference_scaled_smoke_state_000001'
 glob_rank = 0
 glob_nprocs = 4
 sub_rank = 0
 sub_nprocs = 4
 resultfile = 'tracking/state_000001/tau1_reference_scaled_smoke_state_000001_p*'
 nDofs = 1
 nElems = 182320
 time_control = {
    min = {
        iter = 1 
    },
    max = {
        iter = 1 
    },
    interval = {
        iter = 1 
    },
    check_iter = 1,
    delay_check = false 
}
 shape = {
    {
        kind = 'all' 
    } 
}
 varsys = {
    systemname = 'fluid',
    variable = {
        {
            name = 'density_phy',
            ncomponents = 1 
        },
        {
            name = 'velocity_phy',
            ncomponents = 3 
        } 
    },
    nScalars = 4,
    nStateVars = 2,
    nAuxScalars = 4,
    nAuxVars = 2 
}
