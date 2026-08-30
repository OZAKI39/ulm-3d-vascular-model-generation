 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'musubi_force_one_step_oracle'
 basename = 'tracking/musubi_force_one_step_oracle_momentum'
 glob_rank = 0
 glob_nprocs = 1
 sub_rank = 0
 sub_nprocs = 1
 resultfile = 'tracking/musubi_force_one_step_oracle_momentum_p*'
 nDofs = 1
 nElems = 64
 time_control = {
    min = {
        iter = 0 
    },
    max = {
        iter = 2 
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
    systemname = 'fluid_incompressible',
    variable = {
        {
            name = 'velocity_phy',
            ncomponents = 3 
        },
        {
            name = 'density_phy',
            ncomponents = 1 
        } 
    },
    nScalars = 4,
    nStateVars = 2,
    nAuxScalars = 4,
    nAuxVars = 2 
}
