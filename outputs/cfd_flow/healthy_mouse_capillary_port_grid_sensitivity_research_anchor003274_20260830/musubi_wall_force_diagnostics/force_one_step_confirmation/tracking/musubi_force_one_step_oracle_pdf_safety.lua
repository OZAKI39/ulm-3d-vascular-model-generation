 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'musubi_force_one_step_oracle'
 basename = 'tracking/musubi_force_one_step_oracle_pdf_safety'
 glob_rank = 0
 glob_nprocs = 1
 sub_rank = 0
 sub_nprocs = 1
 resultfile = 'tracking/musubi_force_one_step_oracle_pdf_safety_p*'
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
            name = 'pdf',
            ncomponents = 19,
            state_varpos = { 1, 2, 3, 4, 5, 6, 7, 8,
                9, 10, 11, 12, 13, 14, 15, 16,
                17, 18, 19 } 
        } 
    },
    nScalars = 19,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2 
}
