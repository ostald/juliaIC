using BenchmarkTools

include("loadElspec.jl");
include("interpolate_temp.jl")
include("ion_prod.jl")
include("ionchem.jl")
using .ionchem
using MAT


#todo
# - clean up & simplify (X, rr, temp_2 in ionchem.ic => really necessary for allocation?)
# - ne not assigned?? => done, check!!!
# - stepfunctions
# - check form/shape of raw data to be interpolated => how is it supplied, how to standardize?
# - eitenne production
# - temperature correction of raw electron density ne = P/(1 + Te/Ti) => in ElSpec??


function ic_iter(iter, resdir)
    #load ELSPEC output to define time, height, ion densities, temperatures, production rates.
    con = loadmat(joinpath(resdir, "ElSpec-iqt_IC_" * string(iter) * ".mat"))
    ts, te, h, nion, T, e_prod = getparams(con)
    ts[1] = ts[1] - 60*60*30
    nh = length(h)
    nt = length(ts)
    np = length(ionchem.particles)
    n0 = zeros(np, nh)

    #interpolation timesteps:
    t_itp = [ts[1]; (ts[2:end-1] + te[2:end-1])./2; te[end]]

    #assign densities
    #nion = [ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp]
    nion_mapping = ["e-", "N2", "O2", "O", "Ar", "NO+", "O2+", "O+(4S)"] #mapping must correspond to ionchem.particle names
    for (ind1, ion) in enumerate(nion_mapping)
        ind2 = findfirst(==(ion), [p[2] for p in ionchem.particles])
        if nothing == ind2 println(ion * " not found")
        #assign initial densities
        else n0[ind2, :] = nion[ind1, :, 1]
        end
    end

    #interpolate temperatures
    temp_itp = interpolate_temp(t_itp, T)

    #interpolate production and assign production rates to species
    e_prod_itp = interpolate_q(t_itp, e_prod)
    ni_prod = assign_prod(e_prod_f, e_prod_itp, ionchem.particles, n0)

    #integration period
    tspan = (ts[1], te[end])
    #tspan = (0, 10)

    @time sol = ionchem.ic(tspan, n0, ni_prod, temp_itp, nh, (ts + te)./2)

    ni = stack(sol.u, dims =1)

    particles = ionchem.particles
    nN2 = ni[:, findall(p -> p[2] == "N2", particles)[1], :]';
    nO2 = ni[:, findall(p -> p[2] == "O2", particles)[1], :]';
    nO  = ni[:, findall(p -> p[2] == "O" , particles)[1], :]';
    nAr = nion[7, :, :];
    nNOp = ni[:, findall(p -> p[2] == "NO+", particles)[1], :]';
    nO2p = ni[:, findall(p -> p[2] == "O2+", particles)[1], :]';
    nOp_4S = ni[:, findall(p -> p[2] == "O+(4S)", particles)[1], :]';


    #save output, yet to be ordered (code untested)
    elspec_iri_sorted = permutedims(stack([con["iri"][:, 1, :], con["iri"][:, 2, :], con["iri"][:, 3, :], nN2, nO2, nO, nAr, nNOp, nO2p, nOp_4S]), (1, 3, 2))
    eff_rr =   (ionchem.rrates[1](T, ones(nt, nh))' .* ni[:, findall(p -> p[2] == "O2+", particles)[1], :]' .+
                ionchem.rrates[2](T, ones(nt, nh))' .* ni[:, findall(p -> p[2] == "N2+", particles)[1], :]' .+
                ionchem.rrates[3](T, ones(nt, nh))' .* ni[:, findall(p -> p[2] == "NO+", particles)[1], :]')./ni[:, findall(p -> p[2] == "e-", particles)[1], :]'
    #mdict = {"elspec_iri_sorted": elspec_iri_sorted, "eff_rr": eff_rr, "ne_init": ne_init}
    #spio.savemat(direc + 'IC_' + str(iteration) + '.mat', mdict)

    file = joinpath(resdir, "IC_" * string(iter) * ".mat")
    matwrite(file, Dict("elspec_iri_sorted" => elspec_iri_sorted, "eff_rr" => eff_rr);)

end

#be careful; plots can be generated without transposing, but will look wierd
#using Plots
#heatmap(sol.t, h, ni[:, 2, :]')
#heatmap(ts, h, e_prod')



#plot!(range(ts[50],te[100],step=1e-2), [stepf(e_prod[1, :], t, ts, te) for t in range(ts[50],te[100],step=1e-2)])
#plot!(ts[50:100], e_prod[1, 50:100])


"""
When iterating over all the indices for an array, it is better to iterate over eachindex(A) 
    instead of 1:length(A). Not only will this be faster in cases where A is IndexCartesian, 
    but it will also support arrays with custom indexing, such as OffsetArrays. If only the 
    values are needed, then is better to just iterate the array directly, i.e. for a in A.
"""
