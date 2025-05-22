using CairoMakie
using Statistics: mean
include("loadElspec.jl")


#convergence in alpha:
resdir = "/home/oliver/Documents/Results/2006-12-12_wtf7";
files = filter(contains(r"ElSpec-iqt_IC_.*.mat"), readdir(resdir, join=true))

effrr_old = 0

fig, ax, sc = scatter(1, 1, axis=(yscale=log10, limits=(nothing, (1e-7, 1e1)),), color=(:red, 0))

for iter in 1:length(files)-1
        if 1== iter
                matfile = joinpath(resdir, "ElSpec-iqt_IC_"*string(iter-1)*".mat")
                con = loadmat(matfile)
                global effrr_old = con["alpha"]
        end
        matfile = joinpath(resdir, "ElSpec-iqt_IC_"*string(iter)*".mat")
        con = loadmat(matfile)
        effrr = con["alpha"] 
        
        var_effrr = abs.((effrr_old .- effrr)) ./ effrr

        scatter!(iter, maximum(var_effrr), color=:blue)
        scatter!(iter, mean(var_effrr), color=:black)

        effrr_old = effrr
end

current_figure()

##
#____________________________________

