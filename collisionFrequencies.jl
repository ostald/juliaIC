
# Momentum transfer collision frequencies for resonant ion–neutral interactions.5,10 Densities are in cm−3
# H+, H         > 50    2.65 × 10−10 n(H)   Tr^ 1/2  * (1 − 0.083 log10 Tr)^2  
# He+, He       > 50    8.73 × 10−11 n(He)  Tr^ 1/2  * (1 − 0.093 log10 Tr)^2  
# N+, N         > 275   3.83 × 10−11 n(N)   Tr^ 1/2  * (1 − 0.063 log10 Tr)^2  
# O+, O         > 235   3.67 × 10−11 n(O)   Tr^ 1/2  * (1 − 0.064 log10 Tr)^2  
# N2+ , N2      > 170   5.14 × 10−11 n(N2)  Tr^ 1/2  * (1 − 0.069 log10 Tr)^2  
# O2+ , O2      > 800   2.59 × 10−11 n(O2)  Tr^ 1/2  * (1 − 0.073 log10 Tr)^2  
# H+, O         > 300   6.61 × 10−11 n(O)   Ti^ 1/2  * (1 − 0.047 log10 Ti)^2  
# O+, H         > 300   4.63 × 10−12 n(H)   (Tn + Ti/16)^1/2  
# CO+, CO       > 525   3.42 × 10−11 n(CO)  Tr^ 1/2  * (1 − 0.085 log10 Tr)^2  
# CO2+, CO2     > 850   2.85 × 10−11 n(CO2) Tr^ 1/2  * (1 − 0.083 log10 Tr)^2


#where n_s is given in m^-3
nu_Hp_H(nH, Tr)         = 2.65e-16 .* nH   .* Tr.^(1/2) .* (1 .- 0.083 .* log10.(Tr)).^2
nu_Hep_He(nHe, Tr)      = 8.73e-17 .* nHe  .* Tr.^(1/2) .* (1 .- 0.093 .* log10.(Tr)).^2
nu_Np_N(nN, Tr)         = 3.83e-17 .* nN   .* Tr.^(1/2) .* (1 .- 0.063 .* log10.(Tr)).^2
nu_Op_O(nO, Tr)         = 3.67e-17 .* nO   .* Tr.^(1/2) .* (1 .- 0.064 .* log10.(Tr)).^2
nu_N2p_N2(nN2, Tr)      = 5.14e-17 .* nN2  .* Tr.^(1/2) .* (1 .- 0.069 .* log10.(Tr)).^2
nu_O2p_O2(nO2, Tr)      = 2.59e-17 .* nO2  .* Tr.^(1/2) .* (1 .- 0.073 .* log10.(Tr)).^2
nu_Hp_O(nO, Ti)         = 6.61e-17 .* nO   .* Ti.^(1/2) .* (1 .- 0.047 .* log10.(Ti)).^2    #be aware that here it is Ti! (Schunk and Nagy, P.107)
nu_Op_H(nH, Tr)         = 4.63e-18 .* nH   .* (Tr ./8).^(1/2)                        
nu_COp_CO(nCO, Tr)      = 3.42e-17 .* nCO  .* Tr.^(1/2) .* (1 .- 0.085 .* log10.(Tr)).^2
nu_CO2p_CO2(nCO2, Tr)   = 2.85e-17 .* nCO2 .* Tr.^(1/2) .* (1 .- 0.083 .* log10.(Tr)).^2


#The collision frequency coefficients C_in × 10^10 for nonresonant ion–neutral interactions.
# collision frequency nu_in = C_in * n_n with neutral density n_n in cm^-3
#Ion     H       He      N       O       CO      N2      O2      CO2  
#H+      Ra      10.6    26.1    R       35.6    33.6    32.0    41.4 
#He+     4.71    R       11.9    10.1    16.9    16.0    15.3    20.0 
#C+      1.69    1.71    5.73    4.94    8.74    8.26    8.01    10.7 
#N+      1.45    1.49    R       4.42    7.90    7.47    7.25    9.73 
#O+      R       1.32    4.62    R       7.22    6.82    6.64    8.95 
#CO+     0.74    0.79    2.95    2.58    R       4.24    4.49    6.18 
#N2+     0.74    0.79    2.95    2.58    4.84    R       4.49    6.18 
#NO+     0.69    0.74    2.79    2.44    4.59    4.34    4.27    5.89 
#O2+     0.65    0.70    2.64    2.31    4.37    4.13    R       5.63  
#CO2+    0.47    0.51    2.00    1.76    3.40    3.22    3.18    R

#densities here in m^-3 (thus 1e-16, compare to above 1e10, where the minus is missing, I think)
#nu_Hp_H(nH)       = R    
nu_Hep_H(nH)      = 4.71 .* 1e-16 .* nH
nu_Cp_H(nH)       = 1.69 .* 1e-16 .* nH
nu_Np_H(nH)       = 1.45 .* 1e-16 .* nH
#nu_Op_H(nH)       = R 
nu_COp_H(nH)      = 0.74 .* 1e-16 .* nH
nu_N2p_H(nH)      = 0.74 .* 1e-16 .* nH
nu_NOp_H(nH)      = 0.69 .* 1e-16 .* nH
nu_O2p_H(nH)      = 0.65 .* 1e-16 .* nH
nu_CO2p_H(nH)     = 0.47 .* 1e-16 .* nH
 
nu_Hp_He(nHe)     = 10.6 .* 1e-16 .* nHe
#nu_Hep_He(nHe)    = R 
nu_Cp_He(nHe)     = 1.71 .* 1e-16 .* nHe
nu_Np_He(nHe)     = 1.49 .* 1e-16 .* nHe
nu_Op_He(nHe)     = 1.32 .* 1e-16 .* nHe
nu_COp_He(nHe)    = 0.79 .* 1e-16 .* nHe
nu_N2p_He(nHe)    = 0.79 .* 1e-16 .* nHe
nu_NOp_He(nHe)    = 0.74 .* 1e-16 .* nHe
nu_O2p_He(nHe)    = 0.70 .* 1e-16 .* nHe
nu_CO2p_He(nHe)   = 0.51 .* 1e-16 .* nHe
 
nu_Hp_N(nN)       = 26.1 .* 1e-16 .* nN
nu_Hep_N(nN)      = 11.9 .* 1e-16 .* nN
nu_Cp_N(nN)       = 5.73 .* 1e-16 .* nN
#nu_Np_N(nN)       = R
nu_Op_N(nN)       = 4.62 .* 1e-16 .* nN
nu_COp_N(nN)      = 2.95 .* 1e-16 .* nN
nu_N2p_N(nN)      = 2.95 .* 1e-16 .* nN
nu_NOp_N(nN)      = 2.79 .* 1e-16 .* nN
nu_O2p_N(nN)      = 2.64 .* 1e-16 .* nN
nu_CO2p_N(nN)     = 2.00 .* 1e-16 .* nN
 
#nu_Hp_O(nO)       = R 
nu_Hep_O(nO)      = 10.1 .* 1e-16 .* nO
nu_Cp_O(nO)       = 4.94 .* 1e-16 .* nO
nu_Np_O(nO)       = 4.42 .* 1e-16 .* nO
#nu_Op_O(nO)       = R 
nu_COp_O(nO)      = 2.58 .* 1e-16 .* nO
nu_N2p_O(nO)      = 2.58 .* 1e-16 .* nO
nu_NOp_O(nO)      = 2.44 .* 1e-16 .* nO
nu_O2p_O(nO)      = 2.31 .* 1e-16 .* nO
nu_CO2p_O(nO)     = 1.76 .* 1e-16 .* nO

nu_Hp_CO(nCO)     = 35.6 .* 1e-16 .* nCO 
nu_Hep_CO(nCO)    = 16.9 .* 1e-16 .* nCO 
nu_Cp_CO(nCO)     = 8.74 .* 1e-16 .* nCO 
nu_Np_CO(nCO)     = 7.90 .* 1e-16 .* nCO 
nu_Op_CO(nCO)     = 7.22 .* 1e-16 .* nCO 
#nu_COp_CO(nCO)    = R   
nu_N2p_CO(nCO)    = 4.84 .* 1e-16 .* nCO 
nu_NOp_CO(nCO)    = 4.59 .* 1e-16 .* nCO 
nu_O2p_CO(nCO)    = 4.37 .* 1e-16 .* nCO 
nu_CO2p_CO(nCO)   = 3.40 .* 1e-16 .* nCO 

nu_Hp_N2(nN2)     = 33.6 .* 1e-16 .* nN2
nu_Hep_N2(nN2)    = 16.0 .* 1e-16 .* nN2
nu_Cp_N2(nN2)     = 8.26 .* 1e-16 .* nN2
nu_Np_N2(nN2)     = 7.47 .* 1e-16 .* nN2
nu_Op_N2(nN2)     = 6.82 .* 1e-16 .* nN2
nu_COp_N2(nN2)    = 4.24 .* 1e-16 .* nN2
#nu_N2p_N2(nN2)    = R
nu_NOp_N2(nN2)    = 4.34 .* 1e-16 .* nN2
nu_O2p_N2(nN2)    = 4.13 .* 1e-16 .* nN2
nu_CO2p_N2(nN2)   = 3.22 .* 1e-16 .* nN2

nu_Hp_O2(nO2)     = 32.0 .* 1e-16 .* nO2
nu_Hep_O2(nO2)    = 15.3 .* 1e-16 .* nO2
nu_Cp_O2(nO2)     = 8.01 .* 1e-16 .* nO2
nu_Np_O2(nO2)     = 7.25 .* 1e-16 .* nO2
nu_Op_O2(nO2)     = 6.64 .* 1e-16 .* nO2
nu_COp_O2(nO2)    = 4.49 .* 1e-16 .* nO2
nu_N2p_O2(nO2)    = 4.49 .* 1e-16 .* nO2
nu_NOp_O2(nO2)    = 4.27 .* 1e-16 .* nO2
#nu_O2p_O2(nO2)    = R
nu_CO2p_O2(nO2)   = 3.18 .* 1e-16 .* nO2

nu_Hp_CO2(nCO2)   = 41.4 .* 1e-16 .* nCO2
nu_Hep_CO2(nCO2)  = 20.0 .* 1e-16 .* nCO2
nu_Cp_CO2(nCO2)   = 10.7 .* 1e-16 .* nCO2
nu_Np_CO2(nCO2)   = 9.73 .* 1e-16 .* nCO2
nu_Op_CO2(nCO2)   = 8.95 .* 1e-16 .* nCO2
nu_COp_CO2(nCO2)  = 6.18 .* 1e-16 .* nCO2
nu_N2p_CO2(nCO2)  = 6.18 .* 1e-16 .* nCO2
nu_NOp_CO2(nCO2)  = 5.89 .* 1e-16 .* nCO2
nu_O2p_CO2(nCO2)  = 5.63 .* 1e-16 .* nCO2
#nu_CO2p_CO2(nCO2) = R

make_nu_Hp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr, Ti) = [
      nu_Hp_H(nH, Tr)       
    , nu_Hp_He(nHe) 
    , nu_Hp_N(nN)       
    , nu_Hp_O(nO, Ti)       
    , nu_Hp_CO(nCO) 
    , nu_Hp_N2(nN2)     
    , nu_Hp_O2(nO2)     
    , nu_Hp_CO2(nCO2) 
    ]

make_nu_Hep(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_Hep_H(nH)       
  , nu_Hep_He(nHe, Tr) 
  , nu_Hep_N(nN)       
  , nu_Hep_O(nO)       
  , nu_Hep_CO(nCO) 
  , nu_Hep_N2(nN2)     
  , nu_Hep_O2(nO2)     
  , nu_Hep_CO2(nCO2) 
  ]

make_nu_Cp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2) = [
    nu_Cp_H(nH)       
  , nu_Cp_He(nHe) 
  , nu_Cp_N(nN)       
  , nu_Cp_O(nO)       
  , nu_Cp_CO(nCO) 
  , nu_Cp_N2(nN2)     
  , nu_Cp_O2(nO2)     
  , nu_Cp_CO2(nCO2) 
  ]

make_nu_Np(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_Np_H(nH)       
  , nu_Np_He(nHe) 
  , nu_Np_N(nN, Tr)       
  , nu_Np_O(nO)       
  , nu_Np_CO(nCO) 
  , nu_Np_N2(nN2)     
  , nu_Np_O2(nO2)     
  , nu_Np_CO2(nCO2) 
  ]

make_nu_Op(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_Op_H(nH, Tr)       
  , nu_Op_He(nHe) 
  , nu_Op_N(nN)       
  , nu_Op_O(nO, Tr)       
  , nu_Op_CO(nCO) 
  , nu_Op_N2(nN2)     
  , nu_Op_O2(nO2)     
  , nu_Op_CO2(nCO2) 
  ]

make_nu_COp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_COp_H(nH)       
  , nu_COp_He(nHe) 
  , nu_COp_N(nN)       
  , nu_COp_O(nO)       
  , nu_COp_CO(nCO, Tr) 
  , nu_COp_N2(nN2)     
  , nu_COp_O2(nO2)     
  , nu_COp_CO2(nCO2) 
  ]

make_nu_N2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_N2p_H(nH)       
  , nu_N2p_He(nHe) 
  , nu_N2p_N(nN)       
  , nu_N2p_O(nO)       
  , nu_N2p_CO(nCO) 
  , nu_N2p_N2(nN2, Tr)     
  , nu_N2p_O2(nO2)     
  , nu_N2p_CO2(nCO2) 
  ]

make_nu_NOp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2) = [
    nu_NOp_H(nH)       
  , nu_NOp_He(nHe) 
  , nu_NOp_N(nN)       
  , nu_NOp_O(nO)       
  , nu_NOp_CO(nCO) 
  , nu_NOp_N2(nN2)     
  , nu_NOp_O2(nO2)     
  , nu_NOp_CO2(nCO2) 
  ]


make_nu_O2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_O2p_H(nH)       
  , nu_O2p_He(nHe) 
  , nu_O2p_N(nN)       
  , nu_O2p_O(nO)       
  , nu_O2p_CO(nCO) 
  , nu_O2p_N2(nN2)     
  , nu_O2p_O2(nO2, Tr)     
  , nu_O2p_CO2(nCO2) 
  ]

make_nu_CO2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr) = [
    nu_CO2p_H(nH)       
  , nu_CO2p_He(nHe) 
  , nu_CO2p_N(nN)       
  , nu_CO2p_O(nO)       
  , nu_CO2p_CO(nCO) 
  , nu_CO2p_N2(nN2)     
  , nu_CO2p_O2(nO2)     
  , nu_CO2p_CO2(nCO2, Tr)   
  ]

function sum_nu(nu_is) #sum all collision frequncies, can handle float + matrix, for cases when nCO2 = 0 (i.e. float instead of matrix)
    sum = 0
    for nu in nu_is
        sum = sum .+ nu
    end
    return sum
  end