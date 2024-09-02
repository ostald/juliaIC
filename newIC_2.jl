using DifferentialEquations
using BenchmarkTools

include("juliaIC.jl");
include("loadElspec.jl");
include("interpolate_temp.jl")
include("ion_prod.jl")
include("ionchem.jl")
using .ionchem

#todo
# - clean up
# - ne not assigned??
# - stepfunctions
# - check form/shape of raw data to be interpolated => how is it supplied, how to standardize?
# - eitenne production

#load reactions, define particles etc.
#const path_reactions_file = "test_data/Reaction rates full set ext.txt"
#const dndt, particles, reactions, ode_raw, dndt_str, reactions_str = juliaIC.initIC(path_reactions_file)
#rrates = [r[4] for r in reactions]

#load ELSPEC output to define time, height, ion densities, temperatures, production rates.
con = loadmat("test_data/ElSpec-iqt_IC_0.mat")
ts, te, h, nion, T, e_prod = getparams(con)
nh = length(h)
np = length(ionchem.particles)
n0 = zeros(np, nh)

#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
nion_mapping = ["ne", "N2", "O2", "O", "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to particle names
for (ind1, ion) in enumerate(nion_mapping)
    ind2 = findfirst(==(ion), [p[2] for p in ionchem.particles])
    if nothing == ind2 println(ion * " not found")
    else n0[ind2, :] = nion[ind1, :, 1]
    end
end

#interpolate temperatures
temp_itp = interpolate_temp(ts, T)

#interpolate production
e_prod_itp = interpolate_q(ts, e_prod)

nprod = assign_prod(e_prod_f, ionchem.particles, n0)

tspan = (ts[1], ts[end])
tspan = (0, 10)

sol2 = ionchem.ic(tspan, n0, nprod, temp_itp, nh)
@time sol2 = ionchem.ic(tspan, n0, nprod, temp_itp, nh)

ni = stack(sol2.u, dims =1)

#be carefule; plots can be generated without transposing, but will look wierd
using Plots
heatmap(sol2.t, h, ni[:, 2, :]')
heatmap(ts, h, e_prod')



#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


"""
When iterating over all the indices for an array, it is better to iterate over eachindex(A) 
    instead of 1:length(A). Not only will this be faster in cases where A is IndexCartesian, 
    but it will also support arrays with custom indexing, such as OffsetArrays. If only the 
    values are needed, then is better to just iterate the array directly, i.e. for a in A.
"""
