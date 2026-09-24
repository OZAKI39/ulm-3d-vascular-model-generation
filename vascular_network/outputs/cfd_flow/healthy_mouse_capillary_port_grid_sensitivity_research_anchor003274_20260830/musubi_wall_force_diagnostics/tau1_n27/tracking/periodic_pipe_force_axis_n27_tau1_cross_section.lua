 format = 'asciispatial'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'periodic_pipe_force_axis_n27_tau1'
 basename = 'tracking/periodic_pipe_force_axis_n27_tau1_cross_section'
 glob_rank = 0
 glob_nprocs = 2
 sub_rank = 0
 sub_nprocs = 2
 resultfile = 'tracking/periodic_pipe_force_axis_n27_tau1_cross_section_p*'
 nDofs = 1
 nElems = 1120
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
        kind = 'canoND',
        object = {
            {
                origin = {   13.147413150334976E-06,   13.147413150334976E-06,   14.958834517714462E-06 },
                vec = {
                    {    3.622842734758971E-06,    0.000000000000000E+00,    0.000000000000000E+00 },
                    {    0.000000000000000E+00,    3.622842734758971E-06,    0.000000000000000E+00 } 
                },
                segments = {
                    -1,
                    -1 
                },
                distribution = 'equal' 
            } 
        } 
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
