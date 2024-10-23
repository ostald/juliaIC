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

#modify ts
ts[2:end] = ts[2:end] .+ 5*60 #5min
te[2:end] = ts[2:end] .+ 5*60



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
    rrates, nprod, dndt, temp_itp, temp_2, rr, X = p
    temp_2 = temp(temp_itp, t)
    rr = [r(temp_2, X) for r in rrates]
    for j in axes(n, 1)#1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, rr, 0, n, t)
    end
    nothing
end



rrates = [r[4] for r in reactions]

tspan = (ts[1], ts[end])
#tspan = (0, 1)

#n0t = copy(n0')

X = ones(size(T[1, :, 1]))
temp_2 = temp(temp_itp, 0.1)
rr = [r(temp_2, X) for r in rrates]




prob = ODEProblem(myODEf, n0, tspan, (rrates, nprod, dndt, temp_itp, temp_2, rr, X))

@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);
@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);
@time sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);


n0_ = copy(n0)
@profview myODEf(n0_, n0, (rrates, nprod, dndt, temp_itp, temp_2, rr, X), 0.1)

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
heatmap(sol.t, h, ni[:, 9, :]')
heatmap(ts, h, e_prod')



"""
When iterating over all the indices for an array, it is better to iterate over eachindex(A) 
    instead of 1:length(A). Not only will this be faster in cases where A is IndexCartesian, 
    but it will also support arrays with custom indexing, such as OffsetArrays. If only the 
    values are needed, then is better to just iterate the array directly, i.e. for a in A.
"""

function conductivity(ni, temperatures)
    # Resonant ion neutral Momentum Transfer collision frequencies for conductivities
    # from Schunk and Nagy, p. 107:
    # densities in cm^-3 !!!!
    # H+, H     > 50 K      2.65 × 10e−10 n(H)  Tr^(1/2) (1 − 0.083 log10(Tr))^2
    # He+, He   > 50 K      8.73 × 10e−11 n(He) Tr^(1/2) (1 − 0.093 log10(Tr))^2
    # N+, N     > 275 K     3.83 × 10e−11 n(N)  Tr^(1/2) (1 − 0.063 log10(Tr))^2
    # O+, O     > 235 K     3.67 × 10e−11 n(O)  Tr^(1/2) (1 − 0.064 log10(Tr))^2
    # N2+, N2   > 170 K     5.14 × 10e−11 n(N2) Tr^(1/2) (1 − 0.069 log10(Tr))^2
    # O2+, O2   > 800 K     2.59 × 10e−11 n(O2) Tr^(1/2) (1 − 0.073 log10(Tr))^2
    # H+, O     > 300 K     6.61 × 10e−11 n(O)  Ti^(1/2) (1 − 0.047 log10(Ti))^2
    # O+, H     > 300 K     4.63 × 10e−12 n(H)  (Tn + Ti / 16)^(1/2)
    # CO+, CO   > 525 K     3.42 × 10e−11 n(CO) Tr^(1/2) (1 − 0.085 log10(Tr))^2
    # CO2+, CO2 > 850 K     2.85 × 10e−11 n(CO2)Tr^(1/2) (1 − 0.083 log10(Tr))^2

    # ki = eB/(ν_in m_i)
    # ΣP = n0 e /  B  * ki / (1 + ki^2)

    amu = 1.66054e-27 #kg

    Tr = temperatures
    Ti = Tr
    Tn = Tr

    [nH, nHe, nN, nO, nN2, nO2, nCO, nCO2] = [ni[in, :, :] for in in axes(ni, 1)]
    
    nu_Hp_H    = 2.65e-10 .* nH   .* Tr.^(1/2) .* (1 .- 0.083 .* log(10, Tr)).^2
    nu_Hep_He  = 8.73e-11 .* nHe  .* Tr.^(1/2) .* (1 .- 0.093 .* log(10, Tr)).^2
    nu_Np_N    = 3.83e-11 .* nN   .* Tr.^(1/2) .* (1 .- 0.063 .* log(10, Tr)).^2
    nu_Op_O    = 3.67e-11 .* nO   .* Tr.^(1/2) .* (1 .- 0.064 .* log(10, Tr)).^2
    nu_N2p_N2  = 5.14e-11 .* nN2  .* Tr.^(1/2) .* (1 .- 0.069 .* log(10, Tr)).^2
    nu_O2p_O2  = 2.59e-11 .* nO2  .* Tr.^(1/2) .* (1 .- 0.073 .* log(10, Tr)).^2
    nu_Hp_O    = 6.61e-11 .* nO   .* Ti.^(1/2) .* (1 .- 0.047 .* log(10, Ti)).^2
    nu_Op_H    = 4.63e-12 .* nH   .* (Tn + Ti  / 16)^(1/2)                        
    nu_COp_CO  = 3.42e-11 .* nCO  .* Tr.^(1/2) .* (1 .- 0.085 .* log(10, Tr)).^2
    nu_CO2p_CO2= 2.85e-11 .* nCO2 .* Tr.^(1/2) .* (1 .- 0.083 .* log(10, Tr)).^2

    ki_Hp_H    = e.*B ./ (nu_Hp_H     .*  1 .* amu)
    ki_Hep_He  = e.*B ./ (nu_Hep_He   .*  4 .* amu)
    ki_Np_N    = e.*B ./ (nu_Np_N     .*  7 .* amu)
    ki_Op_O    = e.*B ./ (nu_Op_O     .*  8 .* amu)
    ki_N2p_N2  = e.*B ./ (nu_N2p_N2   .* 14 .* amu)
    ki_O2p_O2  = e.*B ./ (nu_O2p_O2   .* 16 .* amu)
    ki_Hp_O    = e.*B ./ (nu_Hp_O     .*  1 .* amu)
    ki_Op_H    = e.*B ./ (nu_Op_H     .*  8 .* amu)
    ki_COp_CO  = e.*B ./ (nu_COp_CO   .* 12 .* amu)
    ki_CO2p_CO2= e.*B ./ (nu_CO2p_CO2 .* 22 .* amu)

    condP_Hp_H    = n0 .* e .* ki_Hp_H        ./ (B .* (1 .+ ki_Hp_H    .^2))
    condP_Hep_He  = n0 .* e .* ki_Hep_He      ./ (B .* (1 .+ ki_Hep_He  .^2))
    condP_Np_N    = n0 .* e .* ki_Np_N        ./ (B .* (1 .+ ki_Np_N    .^2))
    condP_Op_O    = n0 .* e .* ki_Op_O        ./ (B .* (1 .+ ki_Op_O    .^2))
    condP_N2p_N2  = n0 .* e .* ki_N2p_N2      ./ (B .* (1 .+ ki_N2p_N2  .^2))
    condP_O2p_O2  = n0 .* e .* ki_O2p_O2      ./ (B .* (1 .+ ki_O2p_O2  .^2))
    condP_Hp_O    = n0 .* e .* ki_Hp_O        ./ (B .* (1 .+ ki_Hp_O    .^2))
    condP_Op_H    = n0 .* e .* ki_Op_H        ./ (B .* (1 .+ ki_Op_H    .^2))
    condP_COp_CO  = n0 .* e .* ki_COp_CO      ./ (B .* (1 .+ ki_COp_CO  .^2))
    condP_CO2p_CO2= n0 .* e .* ki_CO2p_CO2    ./ (B .* (1 .+ ki_CO2p_CO2.^2))

    condH_Hp_H    = n0 .* e .* ki_Hp_H    .^2 ./ (B .* (1 .+ ki_Hp_H    .^2))
    condH_Hep_He  = n0 .* e .* ki_Hep_He  .^2 ./ (B .* (1 .+ ki_Hep_He  .^2))
    condH_Np_N    = n0 .* e .* ki_Np_N    .^2 ./ (B .* (1 .+ ki_Np_N    .^2))
    condH_Op_O    = n0 .* e .* ki_Op_O    .^2 ./ (B .* (1 .+ ki_Op_O    .^2))
    condH_N2p_N2  = n0 .* e .* ki_N2p_N2  .^2 ./ (B .* (1 .+ ki_N2p_N2  .^2))
    condH_O2p_O2  = n0 .* e .* ki_O2p_O2  .^2 ./ (B .* (1 .+ ki_O2p_O2  .^2))
    condH_Hp_O    = n0 .* e .* ki_Hp_O    .^2 ./ (B .* (1 .+ ki_Hp_O    .^2))
    condH_Op_H    = n0 .* e .* ki_Op_H    .^2 ./ (B .* (1 .+ ki_Op_H    .^2))
    condH_COp_CO  = n0 .* e .* ki_COp_CO  .^2 ./ (B .* (1 .+ ki_COp_CO  .^2))
    condH_CO2p_CO2= n0 .* e .* ki_CO2p_CO2.^2 ./ (B .* (1 .+ ki_CO2p_CO2.^2))

    condH = sum(condH_Hp_H    ,
                condH_Hep_He  ,
                condH_Np_N    ,
                condH_Op_O    ,
                condH_N2p_N2  ,
                condH_O2p_O2  ,
                condH_Hp_O    ,
                condH_Op_H    ,
                condH_COp_CO  ,
                condH_CO2p_CO2) 
    condH = 1

    return nothing

end

function lookup(particles, p)
    return findall(isequal(p), stack(particles))[1][2]
end

[print(particles[lookup(particles, p)]) for p in ["H", "N(4S)", "O", "N2", "O2"]]

densities = [ni[:, 18, :], ni[:, 8, :], ni[:, 5, :], ni[:, 16, :], ni[:, 14, :], 0, 0]
nuO2 = conductivity(densities, 300)

temperatures = stack([temp(temp_itp, t) for t in sol.t])
Te, Ti, Tn = temperatures
Tr = (Ti + Tn) / 2
    

B = 50000e-9 #T
e = 1.6e-19 #C
mi = 2.66e-26 #kg
n0 = ni[:,2, :]

ki = e.*B ./ (nuO2 .* mi)
condP = n0 .* e .* ki    ./ (B .* (1 .+ ki.^2))
condH = n0 .* e .* ki.^2 ./ (B .* (1 .+ ki.^2))

heatmap(sol.t, h, condP')

