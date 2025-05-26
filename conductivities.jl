using SatelliteToolboxGeomagneticField
using LinearAlgebra
using JLD2

#also do non-resonant collisions!!!
# to do:
# - magnetic field model (B not constant, but check on variation!!) => done
# - manage constants
# - load and save densities
# - load remperatures

qe = 1.60217663e-19 #Coulomb

@load "/mnt/data/oliver/alfventrain474s/ic_activeIonosphere.jld2"

"""
data = load("/mnt/data/oliver/alfventrain474s/ic_coldIonosphere.jld2")
particles = data["particles"]
ts = data["ts"]
#te = data["te"]
ni = data["ni"]
"""


ne = transpose(ni[:, 2, :])
nH = transpose(ni[:, 18, :])
nHe = 0
nN = transpose(ni[:, 7, :] .+ ni[:, 8, :])
nO = transpose(ni[:, 3, :] .+ ni[:, 4, :] .+ ni[:, 5, :])
nCO = 0
nN2 = transpose(ni[:, 16, :])
nO2 = transpose(ni[:, 14, :])
nCO2 = 0


#--- Magnetic field for conductivities ------------
"""
include("loadElspec.jl")

con = loadmat("/Users/ost051/Documents/PhD/IonChem/juliaIC/test_data/ElSpec-iqt_IC_0.mat")

iri = con["iri"]
Tn = iri[:, 1, :]
Ti = iri[:, 2, :]
Te = iri[:, 3, :]
Tr = (Tn .+ Ti)./2

Tn = ones(size(nH)) .* 300
Ti = ones(size(nH)) .* 300
Te = ones(size(nH)) .* 300
Tr = ones(size(nH)) .* 300
"""


#T = [Te, Ti, Tn]
Tn = repeat(T[3], inner=(1, size(ne)[2]))
Ti = repeat(T[2], inner=(1, size(ne)[2]))
Te = repeat(T[1], inner=(1, size(ne)[2]))
Tr = (Ti .+ Tn) ./2

"""
date = sum(con["btime"][1:3] .* [1, 1/100, 1/10000])
loc = con["loc"]
h = con["h"] .* 1e3
h = h[1:end-1]
"""


loc = [76, 5] #hardcoded :(
date = sum([2018, 12, 07] .* [1, 1/100, 1/10000])
R = Val(:geodetic)

B = norm.(igrfd.(date, h, loc[1], loc[2], R))

"""
using Plots
plot(B, h/1e3)
"""
#---------------------------------------------------------

include("collisionFrequencies.jl")

nu_Hp   = sum_nu(  make_nu_Hp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr, Ti))
nu_Hep  = sum_nu( make_nu_Hep(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_Cp   = sum_nu(  make_nu_Cp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2))
nu_Np   = sum_nu(  make_nu_Np(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_Op   = sum_nu(  make_nu_Op(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_COp  = sum_nu( make_nu_COp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_N2p  = sum_nu( make_nu_N2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_NOp  = sum_nu( make_nu_NOp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2))
nu_O2p  = sum_nu( make_nu_O2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_CO2p = sum_nu(make_nu_CO2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))


amu = 1.66054e-27 #kg
mH  =  1*amu
mHe =  2*amu
mC  = 12*amu
mN  = 14*amu
mO  = 16*amu
mCO = 28*amu
mN2 = 28*amu
mNO = 30*amu
mO2 = 32*amu
mCO2= 44*amu


function mobility_coeff(B, nu_i, m_i)
    # ki = eB/(ν_in m_i)
    # i for ion, nu_i is already summed for all neutral species
    k_i = qe .* B ./ (nu_i .* m_i)
    return k_i
end

function Pederson_cond(ne, k_i, B)
    condP = ne .* qe .* k_i ./ (B .* (1 .+ k_i.^2))
end

function Hall_cond(ne, k_i, B)
    condH = ne .* qe .* k_i .^2 ./ (B .* (1 .+ k_i.^2))
end


k_Hp  = mobility_coeff(B, nu_Hp, mH)
pc_Hp = Pederson_cond(ne, k_Hp, B)
hc_Hp = Hall_cond(ne, k_Hp, B)

k_Hep  = mobility_coeff(B, nu_Hep, mHe)
pc_Hep = Pederson_cond(ne, k_Hep, B)
hc_Hep = Hall_cond(ne, k_Hep, B)

k_Cp  = mobility_coeff(B, nu_Cp, mC)
pc_Cp = Pederson_cond(ne, k_Cp, B)
hc_Cp = Hall_cond(ne, k_Cp, B)

k_Np  = mobility_coeff(B, nu_Np, mN)
pc_Np = Pederson_cond(ne, k_Np, B)
hc_Np = Hall_cond(ne, k_Np, B)

k_Op  = mobility_coeff(B, nu_Op, mO)
pc_Op = Pederson_cond(ne, k_Op, B)
hc_Op = Hall_cond(ne, k_Op, B)

k_COp  = mobility_coeff(B, nu_COp, mCO)
pc_COp = Pederson_cond(ne, k_COp, B)
hc_COp = Hall_cond(ne, k_COp, B)

k_N2p  = mobility_coeff(B, nu_N2p, mN2)
pc_N2p = Pederson_cond(ne, k_N2p, B)
hc_N2p = Hall_cond(ne, k_N2p, B)

k_NOp  = mobility_coeff(B, nu_NOp, mNO)
pc_NOp = Pederson_cond(ne, k_NOp, B)
hc_NOp = Hall_cond(ne, k_NOp, B)

k_O2p  = mobility_coeff(B, nu_O2p, mO2)
pc_O2p = Pederson_cond(ne, k_O2p, B)
hc_O2p = Hall_cond(ne, k_O2p, B)

k_CO2p  = mobility_coeff(B, nu_CO2p, mCO2)
pc_CO2p = Pederson_cond(ne, k_CO2p, B)
hc_CO2p = Hall_cond(ne, k_CO2p, B)


pc_total = pc_Hp + pc_Hep + pc_Cp + pc_Np + pc_Op + pc_COp + pc_N2p + pc_NOp + pc_O2p + pc_CO2p
hc_total = hc_Hp + hc_Hep + hc_Cp + hc_Np + hc_Op + hc_COp + hc_N2p + hc_NOp + hc_O2p + hc_CO2p

using Plots
p = Plots
fig = heatmap(tsol, 
        h./1e3, 
        max.(-23, log10.(pc_total)),
        clims=(-23, -20),
        ylabel = "Height [km]",
        xlabel = "Time [s]",
        colorbar_title = "log10 Pederson Conductivity []",
        xlim = (0, 3),
        ylim = (100, 200) 
        )

fig = heatmap(tsol, 
        h./1e3, 
        max.(-23, log10.(hc_total)),
        #clims=(-23, -20),
        ylabel = "Height [km]",
        xlabel = "Time [s]",
        colorbar_title = "log10 Hall Conductivity []",
        xlim = (0, 3),
        ylim = (100, 400) 
        )


plot( ne[:, 1]   , h./1e3, xscale=:log10, xlabel = "Electron Density [m-3]", ylabel="Height [km]", label="0 s")
plot!(ne[:, 1000], h./1e3, label="{1 s")
plot!(ne[:, 3000], h./1e3, label="{3 s")



#electron density
heatmap(tsol,
        h./1e3, 
        log10.(ne),
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow,
        #clims=(-12, -7),
        colorbar_title = "log10 Electron Density [m⁻³]",
        #colorscale=log10
        xlims=(0, 0.5)
        )
