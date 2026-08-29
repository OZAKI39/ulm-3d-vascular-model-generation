 format = 'ascii'
 solver = 'Musubi_v2.0.0-4-g4e8b27'
 simname = 'roi003274_adaptive_flux_steady'
 basename = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/mcclure_adaptive_flux_steady_anchor003274_20260829_142633/tracking/p/roi003274_adaptive_flux_steady_p'
 glob_rank = 0
 glob_nprocs = 8
 sub_rank = 0
 sub_nprocs = 8
 resultfile = '/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/outputs/cfd_flow/mcclure_adaptive_flux_steady_anchor003274_20260829_142633/tracking/p/roi003274_adaptive_flux_steady_p_p*'
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
            name = 'pressure_phy',
            ncomponents = 1 
        } 
    },
    nScalars = 1,
    nStateVars = 1,
    nAuxScalars = 4,
    nAuxVars = 2 
}
