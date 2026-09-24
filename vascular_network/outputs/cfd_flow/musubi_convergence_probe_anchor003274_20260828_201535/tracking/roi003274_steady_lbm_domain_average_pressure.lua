 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'roi003274_steady_lbm'
 basename = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/musubi_convergence_probe_anchor003274_20260828_201535/tracking/roi003274_steady_lbm_domain_average_pressure'
 glob_rank = 0
 glob_nprocs = 8
 sub_rank = 0
 sub_nprocs = 8
 resultfile = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/musubi_convergence_probe_anchor003274_20260828_201535/tracking/roi003274_steady_lbm_domain_average_pressure_p*'
 nDofs = 1
 nElems = 221109
 time_control = {
    min = {
        iter = 162464 
    },
    max = {
        iter = 167464 
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
            name = 'pressure_phy',
            ncomponents = 1 
        } 
    },
    nScalars = 1,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2 
}
