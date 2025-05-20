using BenchmarkTools

include("loadElspec.jl");
include("interpolate_temp.jl")
include("ion_prod.jl")
include("ionchem.jl")
using .ionchem

#todo
# - clean up & simplify (X, rr, temp_2 in ionchem.ic => really necessary for allocation?)
# - ne not assigned?? => done, check!!!
# - stepfunctions
# - check form/shape of raw data to be interpolated => how is it supplied, how to standardize?
# - eitenne production
# - temperature correction of raw electron density ne = P/(1 + Te/Ti) => in ElSpec??

#load ELSPEC output to define time, height, ion densities, temperatures, production rates.
con = loadmat("test_data/ElSpec-iqt_IC_0.mat")
ts, te, h, nion, T, e_prod = getparams(con)
nh = length(h)
np = length(ionchem.particles)
n0 = zeros(np, nh)

#interpolation timesteps:
t_itp = [ts[1]; (ts[2:end-1] + te[2:end-1])./2; te[end]]

#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
nion_mapping = ["e-", "N2", "O2", "O", "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to ionchem.particle names
for (ind1, ion) in enumerate(nion_mapping)
    ind2 = findfirst(==(ion), [p[2] for p in ionchem.particles])
    if nothing == ind2 println(ion * " not found")
    #assign initial densities
    else n0[ind2, :] = nion[ind1, :, 1]
    end
end

#interpolate temperatures
temp_itp = interpolate_temp(t_itp, T)

#interpolate production and assign production rates to species
e_prod_itp = interpolate_q(t_itp, e_prod)
ni_prod = assign_prod(e_prod_f, e_prod_itp, ionchem.particles, n0)

#integration period
tspan = (ts[1], te[end])
#tspan = (0, 10)

@time sol = ionchem.ic(tspan, n0, ni_prod, temp_itp, nh, (ts + te)./2)

ni = stack(sol.u, dims =1)

#be carefule; plots can be generated without transposing, but will look wierd
# using CairoMakie
# fig, ax, hm = heatmap(sol.t, h, ni[:, 2, :])
# Colorbar(fig)
# heatmap(ts, h, e_prod)



#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


"""
When iterating over all the indices for an array, it is better to iterate over eachindex(A)
    instead of 1:length(A). Not only will this be faster in cases where A is IndexCartesian,
    but it will also support arrays with custom indexing, such as OffsetArrays. If only the
    values are needed, then is better to just iterate the array directly, i.e. for a in A.
"""
