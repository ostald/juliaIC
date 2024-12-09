using MATLAB
include("ic_iter.jl")
include("loadElspec.jl")



#Define some parameters

"""
ppdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf-pp";
fitdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf";
resdir = "/Users/ost051/Documents/PhD/Results/damped_osc_fixed7";
elspecdir = "/Users/ost051/Documents/PhD/ELSPEC";

btime = [2006, 12, 12, 19, 30, 0.0]; #must be float array!
etime = [2006, 12, 12, 19, 35, 0.0];

experiment = "arc1"

"""
ppdir = "/Users/ost051/Documents/PhD/Data/2022-11-02/02112022/2022-11-02_beata_4@uhfb_pp";
fitdir = "/Users/ost051/Documents/PhD/Data/2022-11-02/02112022/2022-11-02_beata_5@uhfa";
resdir = "/Users/ost051/Documents/PhD/Results/2022-11-02_andres_iqtcl2";
elspecdir = "/Users/ost051/Documents/PhD/ELSPEC";

btime = [2022, 11, 02, 17, 05, 0.0]; #must be float array!
etime = [2022, 11, 02, 17, 10, 0.0];

experiment = "beata"


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

for iter in 0:44
    #call Elspec
    mat"ElSpec_IC_iter($iter, $resdir, $ppdir, $fitdir, $btime, $etime, $experiment)"
    #call IC
    ic_iter(iter, resdir)
end


println("finished")
