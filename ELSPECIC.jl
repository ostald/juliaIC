using MATLAB
using Dates
include("ic_iter.jl")
include("loadElspec.jl")

#Define some parameters

to_evaluate = [ [[2005, 09, 30, 20, 10, 0.0], [2005, 09, 30, 20, 10, 0.0]],
                [[2007, 01, 18, 19, 15, 0.0], [2007, 01, 18, 20, 45, 0.0]],
                [[2007, 12, 12, 00, 00, 0.0], [2007, 12, 12, 04, 30, 0.0]],
                [[2007, 02, 06, 20, 00, 0.0], [2007, 02, 06, 21, 40, 0.0]],
                [[2016, 03, 09, 21, 45, 0.0], [2016, 03, 09, 22, 30, 0.0]]
                ]

for tt in to_evaluate
    date_float = tt[1]
    date_str = string(Date(date_float[1:3]...))
    dir_fit = "/Users/ost051/Documents/PhD/Data" #"/mnt/data/bjorn/EISCAT/Analysed/"
    dir_pp = "" #"/mnt/data/oliver/guisdap_res/pp/"
    #println(joinpath(dir_fit, date_str))
    #println(filter(x->occursin(date_str[1:4], x), readdir(dir_fit)))
    println(joinpath(dir_fit, filter(x->occursin(date_str[1:4], x), readdir(dir_fit))...))
    fitdir = (joinpath(dir_fit, filter(x->occursin(date_str[1:4], x), readdir(dir_fit))...))
    ppdir = (joinpath(dir_pp, filter(x->occursin(date_str[1:4], x), readdir(dir_pp))...))
    ppdir = "/mnt/data/oliver/guisdap_res/pp/04312199.mat"
    resdir = joinpath("/mnt/data/oliver/ic/", date_str)
    btime = tt[1]
    etime = tt[2]
end


ppdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf-pp";

fitdir = "/Users/ost051/Documents/PhD/Data/2006-12-12_arc1_4@uhf";
resdir = "/Users/ost051/Documents/PhD/Results/2006-12-12_newLimitDiv";
elspecdir = "/Users/ost051/Documents/PhD/ELSPEC";

btime = [2006, 12, 12, 19, 30, 0.0]; #must be float array!
etime = [2006, 12, 12, 19, 35, 0.0];

experiment = "arc1"

"""
ppdir = "/Users/ost051/Documents/PhD/Data/2022-11-02/02112022/2022-11-02_beata_4@uhfb_pp";
fitdir = "/Users/ost051/Documents/PhD/Data/2022-11-02/02112022/2022-11-02_beata_5@uhfa";
resdir = "/Users/ost051/Documents/PhD/Results/2022-11-02_andres_newLimitDiv";
elspecdir = "/Users/ost051/Documents/PhD/ELSPEC";

btime = [2022, 11, 02, 17, 05, 0.0]; #must be float array!
etime = [2022, 11, 02, 17, 10, 0.0];

experiment = "beata"
"""

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
