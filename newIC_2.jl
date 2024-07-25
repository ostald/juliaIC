using DifferentialEquations
using BenchmarkTools

include("juliaIC.jl");
include("loadElspec.jl");
include("interpolate_temp.jl")
include("ion_prod.jl")

#todo
#nprod ?
#const ts, te, n0, nprod, e_prod, temp, temp_notInterp, nh, ordering = a; ?
#stepfunctions
#e_prod_itp = [Interpolator(ts, e_prod[:, ih]) for ih in 1:size(e_prod, 2)]
#assign all other deinsities too
#eitenne production

#load reactions, define particles etc.
const path_reactions_file = "test_data/Reaction rates full set ext.txt"
const dndt, particles, reactions, ode_raw, dndt_str, reactions_str = juliaIC.initIC(path_reactions_file)

#load ELSPEC output to define time, height, ion densities, temperatures, production rates.
con = loadmat("test_data/ElSpec-iqt_IC_0.mat")
ts, te, h, nion, T, e_prod = getparams(con)
nh = length(h)
np = length(particles)
n0 = zeros(np, nh)

#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
nion_mapping = ["ne", "N2", "O2", "O", "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to particle names
for (ind1, ion) in enumerate(nion_mapping)
    ind2 = findfirst(==(ion), [p[2] for p in particles])
    if nothing == ind2 println(ion * " not found")
    else n0[ind2, :] = nion[ind1, :, 1]
    end
end

#interpolate temperatures
temp_itp = interpolate_temp(ts, T)

#interpolate production
e_prod_itp = interpolate_q(ts, e_prod)

nprod = assign_prod(e_prod_f, particles, n0)


function myODEf(dn, n, p, t)
    rrates, nprod, dndt, temp_itp = p
    temp_2 = temp(temp_itp, t)
    for j in 1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, rrates, temp_2, n, t)
    end
    nothing
end

function myODEf_speed(dn, n, p, t)
    rrates, nprod, dndt, temp_itp = p
    temp_2 = temp(temp_itp, t)
    rr = [r(temp_2) for r in rrates]
    for j in 1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, rr, 0, n, t)
    end
    nothing
end


rrates = [r[4] for r in reactions]

tspan = (ts[1], ts[end])
#tspan = (0, 1)
prob = ODEProblem(myODEf_speed, n0, tspan, (rrates, nprod, dndt, temp_itp))

@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);
@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);
@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);


@profview sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3)


#sol = solve(prob, SciPyDiffEq.BDF(), reltol = 1e-7, abstol = 1e-3)

ni = stack(sol.u, dims =1)


#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


function solveIC_allAtOnce(n0, ts, te, nprod_julia, temp, reactions, dndt)
    tspan = (ts, te)
    prob = ODEProblem(dummyf, n0, tspan, (reactions, nprod_julia, dndt, temp))
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save)
    return sol    
end

#be carefule; plots can be generated without transposing, but will look wierd
using Plots
heatmap(sol.t, h, ni[:, 2, :]')
heatmap(ts, h, e_prod')



"""
When iterating over all the indices for an array, it is better to iterate over eachindex(A) 
    instead of 1:length(A). Not only will this be faster in cases where A is IndexCartesian, 
    but it will also support arrays with custom indexing, such as OffsetArrays. If only the 
    values are needed, then is better to just iterate the array directly, i.e. for a in A.
"""