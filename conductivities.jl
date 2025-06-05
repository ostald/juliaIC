using SatelliteToolboxGeomagneticField
using LinearAlgebra
include("ic_io.jl")

using CairoMakie
cm = CairoMakie
cm.set_theme!(Theme(colormap = :hawaii))

# to do:
# - manage constants

qe = 1.60217663e-19 #Coulomb

tsol, ni, h, T, e_prod, particles, ts = load_ic("/mnt/data/oliver/alfventrain474s/ic_coldIonosphere.jld2")

"""
data = load("/mnt/data/oliver/alfventrain474s/ic_coldIonosphere.jld2")
particles = data["particles"]
ts = data["ts"]
#te = data["te"]
ni = data["ni"]
"""

# Densities might need transposing

nN    = nN_4S  .+ nN_2D 
nO    = nO     .+ nO_1D  .+ nO_1S 
nOp   = nOp_4S .+ nOp_2D .+ nOp_2P
nO2p  = nO2p   .+ nO2p_a4P

nHe   = 0
nCO   = 0
nCO2  = 0
nHep  = 0
nCp   = 0
nCOp  = 0
nCO2p = 0

#implement iri densities

##
# check ion densities
ixt = 3000

fig, ax, lin = lines(ne[:, ixt],
                        h./1e3,    
                        axis=(xlabel = "Density [m-3]",
                                ylabel = "Height [km]", 
                                title = "$(tsol[ixt]) s", 
                                xscale = log10,
                                limits = ((1, 1e13), nothing)
                                ),
                        label="e",
                    )
#lines!(nHp[:, ixt], h./1e3, label="H⁺")
#lines!(nHep[:, ixt], h./1e3, label="He⁺")
#lines!(nCp[:, ixt], h./1e3, label="C⁺")
lines!(nNp[:, ixt], h./1e3, label="N⁺")
lines!(nOp[:, ixt], h./1e3, label="O⁺")
#lines!(nCOp[:, ixt], h./1e3, label="CO⁺")
lines!(nN2p[:, ixt], h./1e3, label="N₂⁺")
lines!(nNOp[:, ixt], h./1e3, label="NO⁺")
lines!(nO2p[:, ixt], h./1e3, label="O₂⁺")
#lines!(nCO2p[:, ixt], h./1e3, label="CO₂⁺")
lines!((nHp .+ nHep .+ nCp .+ nNp .+ nOp .+ nCOp .+ nN2p .+ nNOp .+ nO2p .+ nCO2p)[:, ixt], h./1e3, label="i⁺", linestyle=:dot )

axislegend()
display(fig)
##

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

B = norm.(igrfd.(date, h, loc[1], loc[2], R)) .*1e-9

lines(B, h/1e3)

#---------------------------------------------------------

include("collisionFrequencies.jl")


nu_e    = sum_nu(   make_nu_e(nH, nHe, nO, nCO, nN2, nO2, nCO2, Te))
nu_Hp   = sum_nu(  make_nu_Hp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr, Ti))
nu_Hep  = sum_nu( make_nu_Hep(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_Cp   = sum_nu(  make_nu_Cp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2))
nu_Np   = sum_nu(  make_nu_Np(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_Op   = sum_nu(  make_nu_Op(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr, Ti, Tn))
nu_COp  = sum_nu( make_nu_COp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_N2p  = sum_nu( make_nu_N2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_NOp  = sum_nu( make_nu_NOp(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2))
nu_O2p  = sum_nu( make_nu_O2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))
nu_CO2p = sum_nu(make_nu_CO2p(nH, nHe, nN, nO, nCO, nN2, nO2, nCO2, Tr))

#nu_Op = 0 .* nu_Op

##
f = Figure()
ax = f[1, 1] = Axis(f, xscale = log10, xlabel="Collision Frequency [s⁻¹]", ylabel="Height [km]", limits=((1e-6, 1e5), nothing))
coll_freqs = [nu_e, nu_Hp, nu_Hep, nu_Cp, nu_Np, nu_Op, nu_COp, nu_N2p, nu_NOp, nu_O2p, nu_CO2p]
labels = ["e", "Hp  ","Hep ","Cp  ","Np  ","Op  ","COp ","N2p ","NOp ","O2p ","CO2p"]
for idx in axes(coll_freqs, 1)
    lines!(coll_freqs[idx][:, 1], h/1e3, label=labels[idx])
end
axislegend(ax)
display(f)


##
f = Figure()
ax = f[1, 1] = Axis(f, xscale = log10, xlabel="Collision Frequency [s⁻¹]", ylabel="Height [km]", limits=((1e-6, 1e5), nothing))
lines!(nu_Op_O(nO, Tr)[:, 1], h/1e3, label = "Op_O")
lines!(nu_Op_O2(nO2)[:, 1], h/1e3, label = "Op_O2")
lines!(nu_O2p_O2(nO2, Tr)[:, 1], h/1e3, label = "O2p_O2")
lines!(nu_N2p_N2(nN2, Tr)[:, 1], h/1e3, label = "N2p_N2")
lines!(nu_O2p_N2(nN2)[:, 1], h/1e3, label = "O2p_N2")
axislegend(ax)
display(f)



amu = 1.66054e-27 #kg
me = 9.1093837e-31 #kg
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


##
f = Figure()
ax = f[1, 1] = Axis(f, xscale = log10, xlabel="Collision Frequency [s⁻¹]", ylabel="Height [km]", limits=((1, 1e8), (100, 200)))
lines!(nu_e[:, 1], h/1e3, label = "e")
lines!(nu_O2p[:, 1], h/1e3, label = "O2p")
lines!(nu_N2p[:, 1], h/1e3, label = "N2p")
lines!(nu_NOp[:, 1], h/1e3, label = "NOp")
lines!(nu_Op[:, 1], h/1e3, label = "Op")
lines!(qe .* B ./ mO2, h/1e3, label="O2+ cy")
lines!(qe .* B ./ me, h/1e3, label="e- cy")
axislegend(ax)
display(f)

##
function mobility_coeff(B, nu_i, m_i)
    # ki = eB/(ν_in m_i)
    # i for ion, nu_i is already summed for all neutral species
    k_i = qe .* B ./ (nu_i .* m_i)
    return k_i
end

function Pederson_cond(n_i, k_i, B)
    condP = n_i .* qe .* k_i ./ (B .* (1 .+ k_i.^2))
end

function Hall_cond(n_i, k_i, B)
    condH = n_i .* qe .* k_i .^2 ./ (B .* (1 .+ k_i.^2))
end

k_e = mobility_coeff(B, nu_e, me) #double check!!!
pc_e = Pederson_cond(ne, k_e, B)
hc_e =  Hall_cond(ne, k_e, B)

k_Hp  = mobility_coeff(B, nu_Hp, mH)
pc_Hp = Pederson_cond(nHp, k_Hp, B)
hc_Hp = Hall_cond(nHp, k_Hp, B)

k_Hep  = mobility_coeff(B, nu_Hep, mHe)
pc_Hep = Pederson_cond(nHep, k_Hep, B)
hc_Hep = Hall_cond(nHep, k_Hep, B)

k_Cp  = mobility_coeff(B, nu_Cp, mC)
pc_Cp = Pederson_cond(nCp, k_Cp, B)
hc_Cp = Hall_cond(nCp, k_Cp, B)

k_Np  = mobility_coeff(B, nu_Np, mN)
pc_Np = Pederson_cond(nNp, k_Np, B)
hc_Np = Hall_cond(nNp, k_Np, B)

k_Op  = mobility_coeff(B, nu_Op, mO)
pc_Op = Pederson_cond(nOp, k_Op, B)
hc_Op = Hall_cond(nOp, k_Op, B)

k_COp  = mobility_coeff(B, nu_COp, mCO)
pc_COp = Pederson_cond(nCOp, k_COp, B)
hc_COp = Hall_cond(nCOp, k_COp, B)

k_N2p  = mobility_coeff(B, nu_N2p, mN2)
pc_N2p = Pederson_cond(nN2p, k_N2p, B)
hc_N2p = Hall_cond(nN2p, k_N2p, B)

k_NOp  = mobility_coeff(B, nu_NOp, mNO)
pc_NOp = Pederson_cond(nNOp, k_NOp, B)
hc_NOp = Hall_cond(nNOp, k_NOp, B)

k_O2p  = mobility_coeff(B, nu_O2p, mO2)
pc_O2p = Pederson_cond(nO2p, k_O2p, B)
hc_O2p = Hall_cond(nO2p, k_O2p, B)

k_CO2p  = mobility_coeff(B, nu_CO2p, mCO2)
pc_CO2p = Pederson_cond(nCO2p, k_CO2p, B)
hc_CO2p = Hall_cond(nCO2p, k_CO2p, B)




f = Figure()
ax = f[1, 1] = Axis(f, xscale = log10, xlabel="Mobility Coefficient [1]", ylabel="Height [km]")
lines!(k_O2p[:, 1], h/1e3, label = "O₂⁺")
lines!(k_N2p[:, 1], h/1e3, label = "N₂⁺")
lines!(k_NOp[:, 1], h/1e3, label = "NO⁺")
lines!(k_Op[:, 1], h/1e3, label =  "O⁺")
axislegend(ax)
display(f)

##
f = Figure()
ax = f[1, 1] = Axis(f, xscale = log10, xlabel="O₂⁺ Conductivity dependence on Mobility [1]", ylabel="Height [km]")
lines!(k_O2p[:, 1]./(1 .+ k_O2p[:, 1].^2) , h/1e3, label = "Pederson")
lines!(k_O2p[:, 1].^2 ./(1 .+ k_O2p[:, 1].^2) , h/1e3, label = "Hall")
axislegend(ax)
display(f)

##


pc_total = pc_e +  pc_Hp + pc_Hep + pc_Cp + pc_Np + pc_Op + pc_COp + pc_N2p + pc_NOp + pc_O2p + pc_CO2p
hc_total = hc_e - (hc_Hp + hc_Hep + hc_Cp + hc_Np + hc_Op + hc_COp + hc_N2p + hc_NOp + hc_O2p + hc_CO2p)



##
# check neutral densities
ixt = 1000

fig, ax, lin = lines(nH[:, ixt],
                        h./1e3,    
                        axis=(xlabel = "Density [m-3]",
                                ylabel = "Height [km]", 
                                xscale=log10,
                                title = "$(tsol[ixt]) s"),
                        label="H",
                    )
#lines!(nHe[:, ixt], h./1e3, label="He")
#lines!(nC[:, ixt], h./1e3, label="C")
lines!(nN[:, ixt], h./1e3, label="N")
lines!(nO[:, ixt], h./1e3, label="O")
#lines!(nCO[:, ixt], h./1e3, label="CO")
lines!(nN2[:, ixt], h./1e3, label="N₂")
lines!(nNO[:, ixt], h./1e3, label="NO")
lines!(nO2[:, ixt], h./1e3, label="O₂")
#lines!(nCO2[:, ixt], h./1e3, label="CO2")

axislegend()
display(fig)



##
# check ion densities
ixt = 1000

fig, ax, lin = lines(ne[:, ixt],
                        h./1e3,    
                        axis=(xlabel = "Density [m⁻³]",
                                ylabel = "Height [km]", 
                                title = "$(round(tsol[ixt], digits=3)) s",
                                xscale=log10,
                                limits=((1, 1e12), nothing)
                                ),
                        label="e",
                    )
#lines!(nHp[:, ixt], h./1e3, label="H⁺")
#lines!(nHep[:, ixt], h./1e3, label="He⁺")
#lines!(nCp[:, ixt], h./1e3, label="C⁺")
lines!(nNp[:, ixt], h./1e3, label="N⁺")
lines!(nOp[:, ixt], h./1e3, label="O⁺")
#lines!(nCOp[:, ixt], h./1e3, label="CO⁺")
lines!(nN2p[:, ixt], h./1e3, label="N₂⁺")
lines!(nNOp[:, ixt], h./1e3, label="NO⁺")
lines!(nO2p[:, ixt], h./1e3, label="O₂⁺")
#lines!(nCO2p[:, ixt], h./1e3, label="CO₂⁺")
lines!((nHp .+ nHep .+ nCp .+ nNp .+ nOp .+ nCOp .+ nN2p .+ nNOp .+ nO2p .+ nCO2p)[:, ixt], h./1e3, label="i+", linestyle=:dot )

axislegend()
display(fig)


##
ixt = 1000

fig, ax, lin = scatter(pc_e[:, ixt],
                        h./1e3,    
                        axis=(xlabel = "Pederson Conductivity [S/m]",
                                ylabel = "Height [km]", 
                                title = "$(tsol[ixt]) s",
                                xscale = log10,
                                limits=((1e-15, 1e-5), nothing)
                                ),
                        label="e",
                        marker=:xcross
                    )
#lines!(pc_Hp[:, ixt], h./1e3, label="H+")
#lines!(pc_Hep[:, ixt], h./1e3, label="He+")
#lines!(pc_Cp[:, ixt], h./1e3, label="C+")
lines!(pc_Np[:, ixt], h./1e3, label="N+")
lines!(pc_Op[:, ixt], h./1e3, label="O+")
#lines!(pc_COp[:, ixt], h./1e3, label="CO+")
lines!(pc_N2p[:, ixt], h./1e3, label="N2+")
lines!(pc_NOp[:, ixt], h./1e3, label="NO+")
lines!(pc_O2p[:, ixt], h./1e3, label="O2+")
#lines!(pc_CO2p[:, ixt], h./1e3, label="CO2+")
lines!((pc_Hp .+ pc_Hep .+ pc_Cp .+ pc_Np .+ pc_Op .+ pc_COp .+ pc_N2p .+ pc_NOp .+ pc_O2p .+ pc_CO2p)[:, ixt], h./1e3, label="i+", linestyle=:dot )

axislegend()
display(fig)

##

ixt = 1000

fig, ax, lin = lines(hc_e[:, ixt],
                        h./1e3,    
                        axis=(xlabel = "Hall Conductivity [S/m]",
                                ylabel = "Height [km]", 
                                title = "$(tsol[ixt]) s",
                                xscale = log10,
                                limits=((1e-15, 1e-4), nothing)
                                ),
                        label="e",
                    )
#lines!(hc_Hp[:, ixt], h./1e3, label="H+")
#lines!(hc_Hep[:, ixt], h./1e3, label="He+")
#lines!(hc_Cp[:, ixt], h./1e3, label="C+")
lines!(hc_Np[:, ixt], h./1e3, label="N+")
lines!(hc_Op[:, ixt], h./1e3, label="O+")
#lines!(hc_COp[:, ixt], h./1e3, label="CO+")
lines!(hc_N2p[:, ixt], h./1e3, label="N2+")
lines!(hc_NOp[:, ixt], h./1e3, label="NO+")
lines!(hc_O2p[:, ixt], h./1e3, label="O2+")
#lines!(hc_CO2p[:, ixt], h./1e3, label="CO2+")
lines!((hc_Hp .+ hc_Hep .+ hc_Cp .+ hc_Np .+ hc_Op .+ hc_COp .+ hc_N2p .+ hc_NOp .+ hc_O2p .+ hc_CO2p)[:, ixt], h./1e3, label="i+", linestyle=:dot )
lines!(hc_total[:, ixt], h./1e3, label="e-i", linestyle=:dot )

axislegend()
display(fig)
##
ixt = 1000

fig, ax, lin = lines(pc_total[:, ixt],
                        h./1e3,    
                        axis=(  xlabel = "Conductivity [S/m]",
                                ylabel = "Height [km]", 
                                xscale=log10,
                                title = "$(tsol[ixt]) s",
                                #limits= ((1e-8, 1e-5), (115, 135))
                                ),
                        label="Pederson"
                    )
#lines!(pc_e[:, ixt],
#        h./1e3,    
#        label="P e"
#        )
lines!(hc_total[:, ixt],
        h./1e3,    
        label="Hall"
        )
#lines!(hc_e[:, ixt],
#        h./1e3,    
#        label="H e"
#        )
ax2 = Axis(fig[1, 1], yticklabelcolor = :red, yaxisposition = :right, xaxisposition= :top,  xlabel = "Density [m-3]",
                ylabel = "Height [km]", 
                xscale=log10,
                #limits= (nothing, (115, 135)),
        )

lines!(ne[:, ixt],
        h./1e3,    
        label="ne",
        )
axislegend()
display(fig)
##

fig, ax, hm = cm.heatmap(tsol,
                            h/1e3,
                            pc_total',
                            axis=(xlabel="Time [s]", 
                                    ylabel="Height [km]",
                                    #limits=(nothing, nothing)
                            ),
                            colorscale = log10 ,
                            colorrange = (1e-8, 1e-5)
                            )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Pederson Conductivity [S/m]")
cm.display(fig)
##

fig, ax, hm = cm.heatmap(tsol,
                            h/1e3,
                            hc_total',
                            axis=(xlabel="Time [s]", 
                                    ylabel="Height [km]",
                                    limits=(nothing, nothing)
                            ),
                            colorscale = log10 ,
                            colorrange = (1e-8, 1e-5)
                            )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Hall Conductivity [S/m]")
cm.display(fig)

##

fig, ax, hm = cm.heatmap(tsol,
                            h./1e3,
                            pc_total'./hc_total',
                            axis=(xlabel="Time [s]", 
                                    ylabel="Height [km]",
                                    limits=(nothing, (100, 150))
                            ),
                            #colorscale = log10 ,
                            colorrange = (0.2, 5),
                            )
cb = cm.Colorbar(fig[1, 2], 
                hm, 
                label = "Ratio Conductivity [1]")
cm.display(fig)


##
for idxh in [86, 106, 127, 147, 166]
#idxh = 147
fig, ax, hm = cm.lines(tsol,
                        pc_total'[:, idxh],
                        label="Pederson",
                        axis=(xlabel="Time [s]", 
                                ylabel="Conductivity [S/m]",
                                title = "$(round(h[idxh]/1e3, digits=0)) km",
                                limits=(nothing, (0, 12e-6)),
                            ),
                            #colorscale = log10 ,
                            #colorrange = (0.2, 5)
                            )

lines!(tsol, hc_total'[:, idxh], label= "Hall")
axislegend(position = :lt)

ax2 = Axis(fig[1, 1], 
        yticklabelcolor = :red, 
        yaxisposition = :right,  
        ylabel = "Ratio [1]",
        limits=(nothing, (0, 6))
        )
lines!(tsol, 
        (pc_total'./hc_total')[:, idxh],
        color=:red, 
        label="Ratio",
        )
axislegend()

cm.display(fig)
end
##


fig, ax, li = cm.lines(
        ne[:, 1]   ,
        h./1e3, 
        axis=(xscale=log10, 
            xlabel = "Electron Density [m-3]", 
            ylabel="Height [km]",
        ), 
        label="$(tsol[1]) s")
cm.lines!(ne[:, 1000], h./1e3, label="$(tsol[1000]) s")
cm.lines!(ne[:, 3000], h./1e3, label="$(tsol[3000]) s")
axislegend()
cm.display(fig)


##

fig, ax, lin = lines(pc_NOp[:, ixt], h/1e3, axis=(xscale = log10,), label="P. cond.")
lines!(k_NOp[:, ixt], h/1e3, label="Mob.")
lines!(nNOp[:, ixt], h/1e3, label="Dens.")
lines!(k_NOp[:, ixt]./(1 .+ k_NOp[:, ixt].^2), h/1e3, label="P. Mob")
lines!(k_NOp[:, ixt]./(1 .+ k_NOp[:, ixt].^2) .* nNOp[:, ixt], h/1e3, label="P. Mob * Dens.")
axislegend()
display(fig)