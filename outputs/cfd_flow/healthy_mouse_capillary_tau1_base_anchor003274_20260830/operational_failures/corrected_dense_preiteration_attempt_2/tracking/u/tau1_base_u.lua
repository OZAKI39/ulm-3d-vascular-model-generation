 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'tau1_base'
 basename = 'tracking/u/tau1_base_u'
 glob_rank = 0
 glob_nprocs = 4
 sub_rank = 0
 sub_nprocs = 4
 resultfile = 'tracking/u/tau1_base_u_p*'
 nDofs = 1
 nElems = 182320
 time_control = {
    min = {
        iter = 0 
    },
    max = {
        iter = 3117935 
    },
    interval = {
        iter = 59875 
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
