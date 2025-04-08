using SatelliteToolboxGeomagneticField
using LinearAlgebra


#also do non-resonant collisions!!!
# to do:
# - magnetic field model (B not constant, but check on variation!!) => done
# - manage constants
# - load and save densities
# - load remperatures




#--- Magnetic field for conductivities ------------
include("loadElspec.jl")

con = loadmat("/Users/ost051/Documents/PhD/IonChem/juliaIC/test_data/ElSpec-iqt_IC_0.mat")

date = sum(con["btime"][1:3] .* [1, 1/100, 1/10000])
loc = con["loc"]
h = con["h"] .* 1e3
R = Val(:geodetic)

B = norm.(igrfd.(date, h, loc[1], loc[2], R))

using Plots
plot(B, h/1e3)
#---------------------------------------------------------

#--- collision rates -------------------------------------
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

#---------------------------------------------------------


amu = 1.66054e-27 #kg
mH  =  1*amu
mHe =  2*amu
mN  = 14*amu
mO  = 16*amu
mN2 = 28*amu
mO2 = 32*amu
mCO = 28*amu
mCO2= 44*amu


function mobility_coeff(B, nu, m)
    # ki = eB/(ν_in m_i)
    k = e .* B ./ (nu .* m)
    return k
end

function Pederson_cond(ne, k, B)
    condP = ne .* e .* k ./ (B .* (1 .+ k.^2))
end

function Hall_cond(ne, k, B)
    condH = ne .* e .* k .^2 ./ (B .* (1 .+ k.^2))
end

