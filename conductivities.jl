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
#where n_s is given in cm^-3
nu_Hp_H    = 2.65e-10 .* nH   .* Tr.^(1/2) .* (1 .- 0.083 .* log(10, Tr)).^2
nu_Hep_He  = 8.73e-11 .* nHe  .* Tr.^(1/2) .* (1 .- 0.093 .* log(10, Tr)).^2
nu_Np_N    = 3.83e-11 .* nN   .* Tr.^(1/2) .* (1 .- 0.063 .* log(10, Tr)).^2
nu_Op_O    = 3.67e-11 .* nO   .* Tr.^(1/2) .* (1 .- 0.064 .* log(10, Tr)).^2
nu_N2p_N2  = 5.14e-11 .* nN2  .* Tr.^(1/2) .* (1 .- 0.069 .* log(10, Tr)).^2
nu_O2p_O2  = 2.59e-11 .* nO2  .* Tr.^(1/2) .* (1 .- 0.073 .* log(10, Tr)).^2
nu_Hp_O    = 6.61e-11 .* nO   .* Ti.^(1/2) .* (1 .- 0.047 .* log(10, Ti)).^2
nu_Op_H    = 4.63e-12 .* nH   .* (Tr ./8).^(1/2)                        
nu_COp_CO  = 3.42e-11 .* nCO  .* Tr.^(1/2) .* (1 .- 0.085 .* log(10, Tr)).^2
nu_CO2p_CO2= 2.85e-11 .* nCO2 .* Tr.^(1/2) .* (1 .- 0.083 .* log(10, Tr)).^2

#The collision frequency coefficients C_in × 10^10 for nonresonant ion–neutral interactions.
# collision frequency nu_in = C_in * n_n with neutral density n_n in cm^-3
Ion     H       He      N       O       CO      N2      O2      CO2  
H+      Ra      10.6    26.1    R       35.6    33.6    32.0    41.4 
He+     4.71    R       11.9    10.1    16.9    16.0    15.3    20.0 
C+      1.69    1.71    5.73    4.94    8.74    8.26    8.01    10.7 
N+      1.45    1.49    R       4.42    7.90    7.47    7.25    9.73 
O+      R       1.32    4.62    R       7.22    6.82    6.64    8.95 
CO+     0.74    0.79    2.95    2.58    R       4.24    4.49    6.18 
N2+     0.74    0.79    2.95    2.58    4.84    R       4.49    6.18 
NO+     0.69    0.74    2.79    2.44    4.59    4.34    4.27    5.89 
O2+     0.65    0.70    2.64    2.31    4.37    4.13    R       5.63  
CO2+    0.47    0.51    2.00    1.76    3.40    3.22    3.18    R

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


function mobility_coeff(B, nu_in, m_i)
    # ki = eB/(ν_in m_i)
    # i for ion, n for neutral species
    k_i = e .* B ./ (nu_in .* m_i)
    return k_i
end

function Pederson_cond(ne, k_i, B)
    condP = ne .* e .* k_i ./ (B .* (1 .+ k_i.^2))
end

function Hall_cond(ne, k_i, B)
    condH = ne .* e .* k_i .^2 ./ (B .* (1 .+ k_i.^2))
end

