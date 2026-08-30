 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'periodic_pipe_force_axis_n27_tau1'
 basename = 'tracking/periodic_pipe_force_axis_n27_tau1_profile'
 glob_rank = 0
 glob_nprocs = 2
 sub_rank = 0
 sub_nprocs = 2
 resultfile = 'tracking/periodic_pipe_force_axis_n27_tau1_profile_p*'
 nDofs = 1
 nElems = 90720
 time_control = {
    min = {
        iter = 0 
    },
    max = {
        iter = 5000 
    },
    interval = {
        iter = 50 
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
    systemname = 'fluid_incompressible',
    variable = {
        {
            name = 'vel_analy',
            ncomponents = 1 
        },
        {
            name = 'vel_error',
            ncomponents = 1 
        } 
    },
    nScalars = 2,
    nStateVars = 2,
    nAuxScalars = 4,
    nAuxVars = 2 
}
