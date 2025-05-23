using MAT
#using Plots
include("interpolate_temp.jl")
include("ion_prod.jl")
include("ionchem.jl")
include("get_msis.jl")
using .ionchem
using Dates

include("interpolate_temp.jl")
include("chemistry.jl");

using DifferentialEquations
using CairoMakie
cm = CairoMakie


#load reactions, define particles etc.
path_reactions_file = "test_data/Reaction rates full set ext.txt"
dndt, particles, reactions, ode_raw, dndt_str, reactions_str = chemistry.initIC(path_reactions_file)
rrates = [r[4] for r in reactions]

path_to_ion_rates = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_474s/Qzt_all_L.mat"
path_to_neutral_atm = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_474s/neutral_atm.mat"

data = matread(path_to_ion_rates)
total_ionization = data["QO2i"] + data["QOi"] + data["QN2i"]
total_ionization = total_ionization[1:end-3, :] * 1e4 # to get units of (/m³/s), because I messed up things a bit
h = data["h_atm"][1:end-3] # in (m)
ts = data["t"] # in (s)

#ts[1] = ts[1] - 1800

heatmap(ts, h./1e3, max.(6, log10.(total_ionization)))#, xlimits=(0,0.1))


#Te, Tn from neutral_atm.mat
neutral_atm =  matread(path_to_neutral_atm)
Te = neutral_atm["Te"][1:end-3]
Tn = neutral_atm["Tn"][1:end-3]
Ti = Tn
T = [Te, Ti, Tn]


ne  = neutral_atm["ne"][1:end-3]
nO  = neutral_atm["nO"][1:end-3]
nO2 = neutral_atm["nO2"][1:end-3]
nN2 = neutral_atm["nN2"][1:end-3]
nion = stack([ne, nN2, nO2, nO])'#, nAr, nNOp, nO2p, nOp]

nh = length(h)
nt = length(ts)
np = length(ionchem.particles)
n0 = zeros(np, nh)

#dt = 0.0001

#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
nion_mapping = ["e-", "N2", "O2", "O"]#, "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to ionchem.particle names
for (ind1, ion) in enumerate(nion_mapping)
    ind2 = findfirst(==(ion), [p[2] for p in particles])
    if nothing == ind2 println(ion * " not found")
    #assign initial densities
    else n0[ind2, :] = nion[ind1, :]
    end
end

n0[findall(p -> p[2] == "NO+", particles)[1], :] = 1/3 .* ne
n0[findall(p -> p[2] == "O2+", particles)[1], :] = 1/3 .* ne
n0[findall(p -> p[2] == "O+(4S)", particles)[1], :] = 1/3 .* ne

loc = [76, 5] #hardcoded :(


#fix prodiction!!
e_prod = total_ionization';
e_prod[1, :] = e_prod[60, :]
#e_prod_itp = [Interpolator(ts, e_prod[:, ih]) for ih in axes(e_prod, 2)]

ni_prod = Array{Any}(undef, length(particles))

zerof = gen_zerof(nh)
ni_prod.= zerof

i_O  = findall(p -> p[2] == "O", particles)[1]
i_O2 = findall(p -> p[2] == "O2", particles)[1]
i_N2 = findall(p -> p[2] == "N2", particles)[1]


ni_prod[findall(p -> p[2] == "e-", particles)[1]]     = (t) -> e_prod[1, :]
ni_prod[findall(p -> p[2] == "O+(4S)", particles)[1]] = (t) -> e_prod[1, :] .* n0[i_O , :]*0.56 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
ni_prod[findall(p -> p[2] == "O2+", particles)[1]]    = (t) -> e_prod[1, :] .* n0[i_O2, :]      ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
ni_prod[findall(p -> p[2] == "N2+", particles)[1]]    = (t) -> e_prod[1, :] .* n0[i_N2, :]*0.92 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])


X = ones(nh)
rr = [r(T, X) for r in rrates]

tspan = (ts[1], ts[end])

#cb_f = []
function cb_f(integrator)
    ind = findfirst(ts .>= integrator.t)
    ni_prod[findall(p -> p[2] == "e-", particles)[1]]     = (t) -> e_prod[ind, :]
    ni_prod[findall(p -> p[2] == "O+(4S)", particles)[1]] = (t) -> e_prod[ind, :] .* n0[i_O , :]*0.56 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
    ni_prod[findall(p -> p[2] == "O2+", particles)[1]]    = (t) -> e_prod[ind, :] .* n0[i_O2, :]      ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
    ni_prod[findall(p -> p[2] == "N2+", particles)[1]]    = (t) -> e_prod[ind, :] .* n0[i_N2, :]*0.92 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
    integrator.p = ni_prod, dndt, rr
    return integrator.p
end

t_cb = ts
cb = PresetTimeCallback(t_cb, cb_f)

function myODEf(dn, n, p, t)
    ni_prod, dndt_, rr = p
    for j in axes(n, 1)
        dn[j, :] .= dndt_[j](ni_prod, rr, 0, n, t)
        #        return (nprodd, rr, tem, nn, tt) -> Base.invokelatest(f, nprodd, rr, tem, nn, tt)
    end
    nothing
end


prob = ODEProblem(myODEf, n0, tspan, (ni_prod, dndt, rr))
sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, callback = cb);
#return sol

ni = stack(sol.u, dims =1)

using JLD2
resdir = "/mnt/data/oliver/alfventrain474s/"
jldsave(joinpath(resdir, "ic_activeIonosphere_moret.jld2"); ts, h, T, e_prod, ni, particles)



cm.set_theme!(Theme(colormap = :hawaii))


# Electron desnity in time and height
fig, ax, hm = cm.heatmap(sol.t,
                h/1e3,
                ni[:, 2, :], 
                axis=(xlabel="Time [s]", 
                        ylabel="Height [km]",
                    )
                #colorscale = log10, 
                #colorrange = (1e5, 1e11)
                )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Electron Density [m⁻³]")
cm.display(fig)

##

# Ionization rate in time and height
fig, ax, hm = cm.heatmap(ts,
                        h/1e3, 
                        e_prod,
                        colorscale = log10, 
                        colorrange = (10^(6), 1e11),
                        axis=(xlabel="Time [s]", 
                            ylabel="Height [km]",
                            )
                        )
cb = cm.Colorbar(fig[1, 2], 
                    hm,
                    label="Ionization rate [m⁻³ s⁻¹]"
                )
cm.display(fig)

##


heatmap(sol.t, h./1e3, log10.(ni[:, 2, :]'))
