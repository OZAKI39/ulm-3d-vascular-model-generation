# Invoke only installed package kinematic routines; no force or integration calls.
using Peridynamics, LinearAlgebra, DelimitedFiles
const PD=Peridynamics
const OUT=abspath(ARGS[1])
X=Matrix(readdlm(joinpath(OUT,"input_X.csv"),','))'
V=vec(readdlm(joinpath(OUT,"input_volume.csv"),','))
ids=Int.(vec(readdlm(joinpath(OUT,"input_representative_ids.csv"),',')))
delta,mu,bulk=vec(readdlm(joinpath(OUT,"input_parameters.csv"),','))
quadratic_kernel(delta,length)=(1-length/delta)^2
# Explicitly disable adaptive floor, Tikhonov, and SVD truncation using the actual API.
mat=RKCMaterial(kernel=quadratic_kernel,model=SaintVenantKirchhoff(),monomial=:PD2,lambda=0.,beta=0.)
@assert mat.epsilon==0 && mat.lambda==0 && mat.beta==0
body=Body(mat,Matrix(X),V)
material!(body;horizon=delta,rho=1056.,G=mu,K=bulk,epsilon_c=100.)
# This solver object selects storage fields only. No initialize!/submit!/step! call.
dh=PD.threads_data_handler(body,VelocityVerlet(steps=1),1)
(;system,storage,paramsetup)=dh.chunks[1]
lookup=Dict(Tuple(X[:,i])=>i for i in axes(X,2))
original=[lookup[Tuple(system.position[:,i])] for i in axes(system.position,2)]
local_ids=[findfirst(==(i),original) for i in ids]
@assert all(storage.bond_active)
weightrows=Vector{Vector{Float64}}()
for (pid,i) in zip(ids,local_ids)
    params=PD.get_params(paramsetup,i)
    PD.rkc_weights!(storage,system,mat,params,0.,0.,i)
    for bondid in PD.each_bond_idx(system,i)
        j=system.bonds[bondid].neighbor
        g=PD.get_vector(storage.gradient_weight,bondid)
        push!(weightrows,[pid,original[j],g[1],g[2],g[3]])
    end
end
writedlm(joinpath(OUT,"julia_gradient_weights.csv"),reduce(vcat,permutedims.(weightrows)),',')
for field in ["axial","shear","cross"]
    input=Matrix(readdlm(joinpath(OUT,"input_"*field*"_x.csv"),','))'
    storage.position .= input[:,original]
    output=zeros(length(ids),9)
    for (k,i) in enumerate(local_ids)
        PD.rkc_defgrad!(storage,system,mat,PD.get_params(paramsetup,i),0.,0.,i)
        F=PD.get_tensor(storage.defgrad,i)
        # Explicit row-major output for NumPy convention.
        output[k,:]=[F[a,b] for a in 1:3 for b in 1:3]
    end
    writedlm(joinpath(OUT,"julia_"*field*"_F.csv"),output,',')
end
println("Julia ",VERSION,"; package ",pkgversion(PD),"; path ",pathof(PD))
println("Full PD2, unchanged quadratic kernel and volumes; epsilon=lambda=beta=0")
println("Representative physical node IDs: ",ids,"; directed weights: ",length(weightrows))
println("KINEMATIC ONLY: no force, stress, Z, equilibrium or time step called")
