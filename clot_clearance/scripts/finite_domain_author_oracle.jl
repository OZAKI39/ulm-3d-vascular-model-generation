# Diagnostic adapter invoking the unmodified, pinned author package, without integration.
using Peridynamics, LinearAlgebra, DelimitedFiles
const PD=Peridynamics
const OUT=abspath(ARGS[1])
X=Matrix(readdlm(joinpath(OUT,"input_X.csv"),','))'
V=vec(readdlm(joinpath(OUT,"input_volume.csv"),','))
delta,mu,bulk=vec(readdlm(joinpath(OUT,"input_parameters.csv"),','))
quadratic_kernel(delta,length)=(1-length/delta)^2
body=Body(RKCMaterial(kernel=quadratic_kernel,model=SaintVenantKirchhoff(),monomial=:C1,epsilon=0),Matrix(X),V)
material!(body;horizon=delta,rho=1056.,G=mu,K=bulk,epsilon_c=100.)
# VelocityVerlet supplies storage layout only; submit!/time integration are never called.
ts=VelocityVerlet(steps=1)
dh=PD.threads_data_handler(body,ts,1)
PD.initialize!(dh,ts)
(;mat,system,storage,paramsetup)=dh.chunks[1]
xin=Matrix(readdlm(joinpath(OUT,"input_x.csv"),','))'
ids=[findfirst(k->all(X[:,k].==system.position[:,i]),axes(X,2)) for i in axes(system.position,2)]
storage.position .= xin[:,ids]
PD.calc_force_density!(dh,0.,0.)
psi=PD.export_field(Val(:strain_energy_density),mat,system,storage,paramsetup,0.)
writedlm(joinpath(OUT,"julia_ids.csv"),ids,',')
writedlm(joinpath(OUT,"julia_force.csv"),transpose(storage.b_int.*transpose(system.volume)),',')
writedlm(joinpath(OUT,"julia_energy.csv"),[sum(psi.*system.volume)],',')
println("INTACT AFFINE ONLY; Julia ",VERSION,"; particles=",length(V),"; directed bonds=",length(system.bonds))
