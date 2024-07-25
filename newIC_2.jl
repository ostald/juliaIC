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
    reactions, nprod, dndt, temp_itp = p
    temp_2 = temp(temp_itp, t)
    for j in 1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, reactions, temp_2, n, t)
    end
    nothing
end


tspan = (ts[1], ts[end])
#tspan = (0, 1)
prob = ODEProblem(myODEf, n0, tspan, (reactions, nprod, dndt, temp_itp))
@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3)



#sol = solve(prob, SciPyDiffEq.BDF(), reltol = 1e-7, abstol = 1e-3)

dn = zeros(size(n0))

dummyf(dn, n0, (reactions, nprod, dndt, temp_itp), 0.1)


@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);
@profview sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3)



#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


function solveIC_allAtOnce(n0, ts, te, nprod_julia, temp, reactions, dndt)
    tspan = (ts, te)
    prob = ODEProblem(dummyf, n0, tspan, (reactions, nprod_julia, dndt, temp))
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save)
    return sol    
end

"""
p = plot(sol[1],
    xaxis = "Time (t)",
    yaxis = ("u(t) (in μm)"),
    label = reshape([p[2] for p in particles], (1, :)),
    )#ylimits = (1, 1e12))
"""


