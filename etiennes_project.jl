using MAT
using Plots
include("interpolate_temp.jl")
include("ion_prod.jl")
include("ionchem.jl")
include("get_msis.jl")
using .ionchem
using Dates


path_to_ion_rates = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_474s/Qzt_all_L.mat"
path_to_neutral_atm = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/Alfven_train_474s/neutral_atm.mat"

data = matread(path_to_ion_rates)
total_ionization = data["QO2i"] + data["QOi"] + data["QN2i"]
total_ionization = total_ionization * 1e4 # to get units of (/m³/s), because I messed up things a bit
z = data["h_atm"] # in (m)
t = data["t"] # in (s)

heatmap(t, z./1e3, max.(6, log10.(total_ionization)))

ts = t
te = t .+ 0.001
ts[1] = ts[1] - 60*30 #set first timestep 30 min in past
dt = 0.0001

h = z

e_prod = total_ionization';

#Te, Tn from neutral_atm.mat
neutral_atm =  matread(path_to_neutral_atm)
Te = neutral_atm["Te"]
Tn = neutral_atm["Tn"]
Ti = Tn


ne  = neutral_atm["ne"]
nO  = neutral_atm["nO"]
nO2 = neutral_atm["nO2"]
nN2 = neutral_atm["nN2"]
nion = stack([ne, nN2, nO2, nO])'#, nAr, nNOp, nO2p, nOp]

loc = [76, 5] #hardcoded :(

t_itp = sort([ts .+ dt; te .- dt])
e_prod = repeat(e_prod, inner = (2, 1))
T = permutedims(repeat(stack([Te, Ti, Tn]), outer = (1, 1, length(t_itp))), (3, 1 ,2))

function temp_f(dummy, t)
    return [Te, Ti, Tn]
end

nh = length(h)
nt = length(ts)
np = length(ionchem.particles)
n0 = zeros(np, nh)

#interpolation timesteps:
#t_itp = [ts[1]; (ts[2:end-1] + te[2:end-1])./2; te[end]]

particles = ionchem.particles
#assign densities
#nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
nion_mapping = ["e-", "N2", "O2", "O"]#, "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to ionchem.particle names
for (ind1, ion) in enumerate(nion_mapping)
    ind2 = findfirst(==(ion), [p[2] for p in ionchem.particles])
    if nothing == ind2 println(ion * " not found")
    #assign initial densities
    else n0[ind2, :] = nion[ind1, :]
    end
end

n0[findall(p -> p[2] == "NO+", particles)[1], :] = 1/3 .* ne
n0[findall(p -> p[2] == "O2+", particles)[1], :] = 1/3 .* ne
n0[findall(p -> p[2] == "O+(4S)", particles)[1], :] = 1/3 .* ne


#interpolate temperatures
temp_itp = interpolate_temp(t_itp, T)

#interpolate production and assign production rates to species
e_prod_itp = interpolate_q(t_itp, e_prod)
ni_prod = assign_prod(e_prod_f, e_prod_itp, ionchem.particles, n0)
#integration period

tspan = (t_itp[1], t_itp[end])
#ts, te, h, nion, T, e_prod, loc = getparams(con)

@time sol = ionchem.ic(tspan, n0, ni_prod, temp_itp, nh, (ts + te)./2)

 #shift first timestep back half an hour

 #no interpolation as of now

