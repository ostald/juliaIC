using MATLAB
include("ic_iter.jl")
include("loadElspec.jl")

#Define some parameters
ppdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf-pp";
fitdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf";
resdir = "/Users/ost051/Documents/PhD/Results/damped_osc_fixed";
elspecdir = "/Users/ost051/Documents/PhD/ELSPEC";

if !isdir(resdir)
    mkdir(resdir);
end

if isfile(joinpath(resdir, "ElSpec-iqt_IC_1.mat"))
    println("resdir not empty; define new resdir?")
    #wait_for_key("press any key to continue")
    input = read(stdin, 1);
    println(" ")
    println("overwriting dir")
end

#call Elspec
mat"addpath($elspecdir)"

for i in 1:44
    iter = i-1
    #call Elspec
    mat"ElSpec_IC_iter($iter, $resdir, $ppdir, $fitdir)"
    #call IC
    ic_iter(iter, resdir)
end




println("finished")
