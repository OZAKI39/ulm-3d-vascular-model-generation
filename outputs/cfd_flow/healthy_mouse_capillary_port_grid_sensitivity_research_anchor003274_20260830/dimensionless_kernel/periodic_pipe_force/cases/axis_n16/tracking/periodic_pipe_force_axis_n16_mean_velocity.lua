 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'periodic_pipe_force_axis_n16'
 basename = 'tracking/periodic_pipe_force_axis_n16_mean_velocity'
 glob_rank = 0
 glob_nprocs = 2
 sub_rank = 0
 sub_nprocs = 2
 resultfile = 'tracking/periodic_pipe_force_axis_n16_mean_velocity_p*'
 nDofs = 1
 nElems = 19968
 time_control = {
    min = {
        iter = 0 
    },
    max = {
        iter = 2000 
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
            name = 'velocity_phy',
            ncomponents = 3 
        } 
    },
    nScalars = 3,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2 
}
