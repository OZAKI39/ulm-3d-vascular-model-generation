# Calls the unmodified author package. No copied BAQP solver source.
using Peridynamics, LinearAlgebra, DelimitedFiles
const PD=Peridynamics
const OUT=abspath(@__DIR__)
X=Matrix(readdlm(joinpath(OUT,"input_X.csv"),','))'
V=vec(readdlm(joinpath(OUT,"input_volume.csv"),','))
meta=vec(readdlm(joinpath(OUT,"input_parameters.csv"),','))
delta,mu,bulk=meta
quadratic_kernel(l,d)=(1-l/d)^2
for name in ["affine","nonaffine"]
    body=Body(RKCMaterial(kernel=quadratic_kernel,model=SaintVenantKirchhoff(),monomial=:C1,epsilon=0),Matrix(X),V)
    material!(body;horizon=delta,rho=1056.,G=mu,K=bulk,epsilon_c=100.)
    ts=VelocityVerlet(steps=1)
    dh=PD.threads_data_handler(body,ts,1)
    PD.initialize!(dh,ts)
    (;mat,system,storage,paramsetup)=dh.chunks[1]
    xin=Matrix(readdlm(joinpath(OUT,"input_x_"*name*".csv"),','))'
    ids=[findfirst(k->all(X[:,k].==system.position[:,i]),axes(X,2)) for i in axes(system.position,2)]
    storage.position .= xin[:,ids]
    PD.calc_force_density!(dh,0.,0.)
    psi=PD.export_field(Val(:strain_energy_density),mat,system,storage,paramsetup,0.)
    writedlm(joinpath(OUT,"julia_"*name*"_ids.csv"),ids,',')
    writedlm(joinpath(OUT,"julia_"*name*"_force.csv"),transpose(storage.b_int.*transpose(system.volume)),',')
    writedlm(joinpath(OUT,"julia_"*name*"_Fbar.csv"),transpose(storage.defgrad),',')
    writedlm(joinpath(OUT,"julia_"*name*"_energy.csv"),[sum(psi.*system.volume)],',')
    bondrows=Vector{Vector{Float64}}()
    for i in PD.each_point_idx(system), b in PD.each_bond_idx(system,i)
        j=system.bonds[b].neighbor
        ids[i]<ids[j] || continue
        xi=PD.get_vector_diff(system.position,i,j)
        eta=PD.get_vector_diff(storage.position,i,j)
        Ft=PD.bond_avg(PD.get_tensor(storage.defgrad,i),PD.get_tensor(storage.defgrad,j),xi,eta,system.bonds[b].length)
        push!(bondrows,vcat(Float64(ids[i]),Float64(ids[j]),vec(Matrix(Ft))))
    end
    writedlm(joinpath(OUT,"julia_"*name*"_Ftilde.csv"),transpose(hcat(bondrows...)),',')
    println(name," energy=",sum(psi.*system.volume)," points=",length(V)," bonds=",length(bondrows))
end
println("Julia version: ",VERSION)
