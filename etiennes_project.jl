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

include("ic_io.jl")

initialize = "loop"
res_filename = "ic_loopIonosphere.jld2"
resdir = "/mnt/data/oliver/alfventrain474s/"

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

e_prod = total_ionization';

neutral_atm =  matread(path_to_neutral_atm)

nh = length(h)
nt = length(ts)
np = length(ionchem.particles)
n0 = zeros(np, nh)

if initialize == "preheat"
    ts[1] = ts[1] - 1800
    e_prod[1, :] = e_prod[101, :]
    ne  = neutral_atm["ne"][1:end-3] .* 0 .+ 1
    n0[findall(p -> p[2] == "NO+", particles)[1], :] = 1/3 .* ne
    n0[findall(p -> p[2] == "O2+", particles)[1], :] = 1/3 .* ne
    n0[findall(p -> p[2] == "O+(4S)", particles)[1], :] = 1/3 .* ne

elseif initialize == "loop"
    _, ni, _, _, _, _, _ = load_ic("/mnt/data/oliver/alfventrain474s/ic_coldIonosphere.jld2")
    n0 = ni[end, :, :]
    #ne = ne[:, end]
    #n0[findall(p -> p[2] == "NO+", particles)[1], :] = nNOp[:, end]
    #n0[findall(p -> p[2] == "O2+", particles)[1], :] = nO2p[:, end]
    #n0[findall(p -> p[2] == "O+(4S)", particles)[1], :] = nOp_4S[:, end]
elseif initialize == "cold"
    ne  = neutral_atm["ne"][1:end-3] .* 0 .+ 1
    n0[findall(p -> p[2] == "NO+", particles)[1], :] = 1/3 .* ne
    n0[findall(p -> p[2] == "O2+", particles)[1], :] = 1/3 .* ne
    n0[findall(p -> p[2] == "O+(4S)", particles)[1], :] = 1/3 .* ne

end

#heatmap(ts, h./1e3, max.(6, log10.(total_ionization)))#, xlimits=(0,0.1))


#Te, Tn from neutral_atm.mat
Te = neutral_atm["Te"][1:end-3]
Tn = neutral_atm["Tn"][1:end-3]
Ti = Tn
T = [Te, Ti, Tn]


if initialize != "loop"
    nO  = neutral_atm["nO"][1:end-3]
    nO2 = neutral_atm["nO2"][1:end-3]
    nN2 = neutral_atm["nN2"][1:end-3]
    nion = stack([ne, nN2, nO2, nO])'#, nAr, nNOp, nO2p, nOp]



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
end

loc = [76, 5] #hardcoded :(

"""
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = "/Library/Frameworks/Python.framework/Versions/3.10/bin/python3"
using PythonCall
iri2016 = pyimport("iri2016 => IRI")
date = sum([2018, 12, 07] .* [1, 1/100, 1/10000])
iri2016.IRI(date.astype(datetime), h, loc[0], loc[1])
"""

max_ionization = stack([findmax(total_ionization[:, idt]) for idt in axes(total_ionization, 2)])

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

filter = unique(i -> round(sol.t[i], digits = 4), eachindex(sol.t))

tsol = sol.t[filter]
ni = stack(sol.u[filter], dims =1)

save_ic(joinpath(resdir, res_filename), tsol, ni, h, T, e_prod, particles, ts)
#jldsave(joinpath(resdir, "ic_coldIonosphere.jld2"); tsol, ni, h, T, e_prod, particles, ts)

assign_densities(ni, particles)
ni_ion = nNOp .+ nO2p .+ nOp_4S .+ nOp_2P .+ nOp_2D .+ nN2p .+ nNp .+ nHp .+ nO2p_a4P



using CairoMakie
cm = CairoMakie
cm.set_theme!(Theme(colormap = :hawaii))

##
# Electron desnity in time and height
fig, ax, hm = cm.heatmap(tsol,
                h/1e3,
                max.(1e5, ni[:, 2, :]), 
                axis=(xlabel="Time [s]", 
                        ylabel="Height [km]",
                        limits=((-10, 3), nothing)
                    ),
                colorscale = log10 ,
                #colorrange = (1e8, 1e12)
                )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Electron Density [m⁻³]")
cm.display(fig)

##

# Plot all desnities in time and height

for pix in particles
fig, ax, hm = cm.heatmap(tsol,
                h/1e3,
                max.(1e5, ni[:, pix[1], :]), 
                axis=(xlabel="Time [s]", 
                        ylabel="Height [km]",
                        #limits=((-1800, 0), nothing)
                    ),
                colorscale = log10 ,
                #colorrange = (1e8, 1e12)
                )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = pix[2] *" Density [m⁻³]")
cm.display(fig)
end

##

# Ionization rate in time and height
fig, ax, hm = cm.heatmap(ts,
                        h/1e3, 
                        e_prod,
                        colorscale = log10, 
                        colorrange = (10^(6), 1e11),
                        axis=(xlabel="Time [s]", 
                            ylabel="Height [km]",
                            limits=((0, 0.5), nothing),
                            )
                        )
cb = cm.Colorbar(fig[1, 2], 
                    hm,
                    label="Ionization rate [m⁻³ s⁻¹]"
                )
cm.display(fig)

##

#Charge balance
fig, ax, hm = cm.heatmap(tsol,
                h/1e3,
                transpose((ni_ion .- ne ) ./ ne) ,
                #max.(1e5, ni[:, pix[1], :]), 
                axis=(xlabel="Time [s]", 
                        ylabel="Height [km]",
                        #limits=((-1800, 0), nothing),
                    ),
                #colorscale = log10 ,
                #colorrange = (1e8, 1e12)
                )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Charge imbalance")
cm.display(fig)

#tsol, ni, h, T, e_prod, particles, ts = load_ic(joinpath(resdir, res_filename))

fig, ax, lin = scatter(0, 0, axis=(title="Charge conservation",),)
for idx in 1:nh
    lines!(tsol, ((ni_ion .- ne ) ./ ne)[idx, :] )
end
display(fig)