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
const dndt, particles, reactions, ode_raw, dndt_str = juliaIC.initIC(path_reactions_file)

#load ELSPEC output to define time, height, ion densities, temperatures, production rates.
con = loadmat("test_data/ElSpec-iqt_IC_0.mat")
ts, te, h, nion, T, e_prod = getparams(con)
nh = length(h)
np = length(particles)
n0 = zeros(np, nh)

#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
for i in 1:np
    if particles[i][2] == "N2"
        n0[i, :] = nion[2, :, 1]
        println(i)
    elseif particles[i][2] == "O2"
        n0[i, :] = nion[3, :, 1]
    elseif particles[i][2] == "O(1D)"
        n0[i, :] = nion[4, :, 1]
    elseif particles[i][2] == "NO+"
        n0[i, :] = nion[6, :, 1]
    elseif particles[i][2] == "O2+"
        n0[i, :] = nion[7, :, 1]
    elseif particles[i][2] == "O+(4S)"
        n0[i, :] = nion[8, :, 1]
    end
end

#revise!
const t_save = range(ts, te, 1001)
const t_step = range(ts[1], te[end], step=1e-2)

#interpolate temperatures
temp_itp = interpolate_temp(ts, T)
@benchmark temp(temp_itp, 1)

#interpolate production
e_prod_itp = interpolate_q(ts, e_prod)
#println(size(e_prod_f(e_prod_itp, 0.1)))
@benchmark e_prod_f(e_prod_itp, 0.1)

#production for all species:
#1. no production
function zerof(t)
    return zeros(62) #idea: instead of hardcoding, use variable nheihgts to automatically generate a hardcoded version
end
#make array of no production
nprod_julia = Array{Any}(undef, length(particles))
nprod_julia.= zerof
#assign production where there is
i_O  = findall(p -> p[2] == "O", particles)[1]
i_O2 = findall(p -> p[2] == "O2", particles)[1]
i_N2 = findall(p -> p[2] == "N2", particles)[1]
for i in 1:length(particles)
    if particles[i][2] == "e-"
        nprod_julia[i] = (t) -> e_prod_f(e_prod_itp, t)
    end
    if particles[i][2] == "O+(4S)"
        nprod_julia[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_O , :]*0.56 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
        println(size(nprod_julia[i](0)))
    end
    if particles[i][2] == "O2+"
        nprod_julia[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_N2, :]*0.92 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
    end
    if particles[i][2] == "N2+"
        nprod_julia[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_O2, :]      ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
    end
end

const nprod_julia_c = nprod_julia
nprod = nprod_julia

# for i in 1:length(particles)
#     print(i, " ")
#     #println(dndt[i])
#     println(size(dndt[i](nprod, reactions, temp(temp_itp, 0.1), n0, 0.1)))
#     println(dndt[i](nprod, reactions, temp(temp_itp, 0.1), n0, 0.1))
# end
# println("chkpnt2")

function dummyf(n, p, t)
    dn = zeros(size(n))
    reactions, nprod, dndt, temp_itp = p
    temp_2 = temp(temp_itp, t)
    #rr = [r[4](temp_2) for r in reactions]
    for j in 1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, reactions, temp_2, n, t)
    end 
    return dn
end

tspan = (ts[1], ts[end])
tspan = (0, 1)
prob = ODEProblem(dummyf, n0, tspan, (reactions, nprod, dndt, temp_itp))
#sol = solve(prob, SciPyDiffEq.BDF(), reltol = 1e-7, abstol = 1e-3)

#dummyf(u0, (reactions, nprod, dndt, temp), 0)


@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3)
@profview sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3)



#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


function solveIC_allAtOnce(n0, ts, te, nprod_julia, temp, reactions, dndt)
    tspan = (ts, te)
    prob = ODEProblem(dummyf, n0, tspan, (reactions, nprod_julia, dndt, temp))
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save)
    return sol    
end


p = plot(sol[1],
    xaxis = "Time (t)",
    yaxis = ("u(t) (in μm)"),
    label = reshape([p[2] for p in particles], (1, :)),
    )#ylimits = (1, 1e12))



