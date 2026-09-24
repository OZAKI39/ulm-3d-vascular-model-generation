 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'a3274'
 basename = 'tracking/u/a3274_u'
 glob_rank = 0
 glob_nprocs = 8
 sub_rank = 0
 sub_nprocs = 8
 resultfile = 'tracking/u/a3274_u_p*'
 nDofs = 1
 nElems = 221309
 time_control = {
    min = {
        iter = 0 
    },
    max = {
        iter = 1000000 
    },
    interval = {
        iter = 100 
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
            name = 'vel_mag_phy',
            ncomponents = 1 
        } 
    },
    nScalars = 1,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2 
}
