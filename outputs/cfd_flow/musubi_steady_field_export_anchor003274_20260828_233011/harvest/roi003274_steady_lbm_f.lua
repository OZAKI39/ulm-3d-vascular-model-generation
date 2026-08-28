 format = 'asciispatial'
 solver = 'mus_harvesting_v2.0.0-4-g4e8b27'
 simname = 'roi003274_steady_lbm'
 basename = '/tmp/u3d/x/roi003274_steady_lbm_f'
 glob_rank = 0
 glob_nprocs = 1
 sub_rank = 0
 sub_nprocs = 1
 resultfile = '/tmp/u3d/x/roi003274_steady_lbm_f_p*'
 nDofs = 1
 nElems = 221109
 time_control = {
    min = { 
    },
    max = { 
    },
    interval = { 
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
            name = 'pressure_phy',
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
