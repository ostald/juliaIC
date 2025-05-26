## In this file: Plots for single iteration
#       - Charge conservation
#       - Densities etc.


using Plots
using .ionchem
particles = ionchem.particles
include("loadElspec.jl")

#load result from IonChem:
resdir = "/home/oliver/Documents/Results/2006-12-12_wtf7";
iter = 20
file = joinpath(resdir, "ic_densities_"*string(iter)*".jld2")
f = load(file)
ni = f["ni"];
particles = f["particles"]
ts = f["ts"]

#automate this
ne       = ni[:, findall(p -> p[2] == "e-"       , particles)[1], :]';
nN2      = ni[:, findall(p -> p[2] == "N2"       , particles)[1], :]';
nO2      = ni[:, findall(p -> p[2] == "O2"       , particles)[1], :]';
nO       = ni[:, findall(p -> p[2] == "O"        , particles)[1], :]';
nNOp     = ni[:, findall(p -> p[2] == "NO+"      , particles)[1], :]';
nO2p     = ni[:, findall(p -> p[2] == "O2+"      , particles)[1], :]';
nOp_4S   = ni[:, findall(p -> p[2] == "O+(4S)"   , particles)[1], :]';
nOp_2D   = ni[:, findall(p -> p[2] == "O+(2D)"   , particles)[1], :]';
nOp_2P   = ni[:, findall(p -> p[2] == "O+(2P)"   , particles)[1], :]';
nN2p     = ni[:, findall(p -> p[2] == "N2+"      , particles)[1], :]';
nNp      = ni[:, findall(p -> p[2] == "N+"       , particles)[1], :]';
nHp      = ni[:, findall(p -> p[2] == "H+"       , particles)[1], :]';
nO2p_a4P = ni[:, findall(p -> p[2] == "O2+(a4P)" , particles)[1], :]';
nionsp = nNOp .+ nO2p .+ nOp_4S .+ nOp_2P .+ nOp_2D .+ nN2p .+ nNp .+ nHp .+ nO2p_a4P



#plot neutral density profiles:
tix = 200 #time index
plt = plot( max.(1e10, nN2[:, tix]), h/1e3, xscale=:log10, label="N2", xlabel = "Density [m⁻³]", ylabel="Height [km]")
plot!(max.(1e10, nO2[:, tix]), h/1e3, label = "O2")
plot!(max.(1e12,  nO[:, tix]), h/1e3, label = "O")
display(plt) #important in loops!!

#electron density
heatmap(ts.-ts[1], 
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



#plot charge conservation, relative
heatmap(ts.-ts[1], 
        h./1e3, 
        log10.(max.(1e-12, abs.((ne .- nionsp) ./ ne))),
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow,
        clims=(-12, -7),
        colorbar_title = "log10 Relative Charge Imbalance [1]",
        #colorscale=log10
        )


#Charge conservation for each height separately
plot(ts.-ts[1], ((ne.-nionsp)./ne)', legend=false)

#N+ Density
heatmap(ts.-ts[1], 
        h./1e3, 
        log10.(nNp),
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow,
        clims=(0, 8),
        colorbar_title = "log10 N+ Density [m-3]",
        )


#Relaltive N+ Density
heatmap(ts.-ts[1], 
        h./1e3, 
        log10.(nNp./ne),
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow,
        #clims=(0, 2e-5),
        clims=(-10, -4),
        colorbar_title = "log10 relative N2 Density [1]",
        )


#variation in O2+:
idx_O2p = 1
for iter in 1:4
        if 1== iter
                file = joinpath(resdir, "ic_densities_"*string(iter-1)*".jld2")
                f = load(file)
                global ni_old = f["ni"]
        end 
        file = joinpath(resdir, "ic_densities_"*string(iter)*".jld2")
        f = load(file)
        ni = f["ni"]

        var_ni = (ni_old .- ni) ./ ni

        plt = heatmap(var_ni[:, idx_O2p, :]')
        display(plt)

        ni_old = ni
end


heatmap(log10.(ni[:, idx_O2p, :])', 
        clims=(10, 13))
        


function plot_density(density, ts, h, clims=[])
        plt = heatmap(ts.-ts[1], 
                h./1e3, 
                log10.(density),
                xlabel="Time [s]", 
                ylabel="Height [km]",
                c=:batlow,
                #clims=(0, 2e-5),
                clims=clims,
        )
        return plt
end

plt = plot_density(nO2p, ts, h, (10, 13))
colorbar_title = "log10 O2+ Density [m-3]"
display(plt)

plot_density(nNOp, ts, h, (9, 11))
plot_density(max.(nN2p, 1), ts, h, (5, 10))

heatmap(nO2p./nNOp,
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow)



con = loadmat(joinpath(resdir, "ElSpec-iqt_IC_"*string(iter)*".mat"))
heatmap(con["iri"][:, 9, :] ./ con["iri"][:, 8, :],
        xlabel="Time [s]", 
        ylabel="Height [km]",
        c=:batlow)




#variation in O2+:
effrr_old = 0

for iter in 1:14
        if 1== iter
                matfile = joinpath(resdir, "ElSpec-iqt_IC_"*string(iter-1)*".mat")
                con = loadmat(matfile)
                effrr_old = con["alpha"]
        end
        matfile = joinpath(resdir, "ElSpec-iqt_IC_"*string(iter)*".mat")
        con = loadmat(matfile)
        effrr = con["alpha"] 
        
        var_effrr = (effrr_old .- effrr) ./ effrr

        plt = heatmap(var_effrr)
        display(plt)

        effrr_old = effrr
end



