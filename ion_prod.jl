using PCHIPInterpolation

function interpolate_q(t, e_prod)
    e_prod_itp = [Interpolator(t, e_prod[:, ih]) for ih in axes(e_prod, 2)]
    return e_prod_itp
end

function e_prod_f(e_prod_itp, t)
    return [i.(t) for i in e_prod_itp]
end

"""
function zerof(t)
    return zeros(62) #idea: instead of hardcoding, use variable nheihgts to automatically generate a hardcoded version
end
"""

function gen_zerof(nh)
    function zerof(t)
        return zeros(nh)
    end
    return zerof
end

#production for all species:
function assign_prod(e_prod_f, e_prod_itp, particles, n0) #n0 sufficient bc major species are constant
    #make array of no production
    ni_prod = Array{Any}(undef, length(particles))
    nh = size(e_prod_itp, 1)
    zerof = gen_zerof(nh)
    ni_prod.= zerof
    #assign production where there is
    i_O  = findall(p -> p[2] == "O", particles)[1]
    i_O2 = findall(p -> p[2] == "O2", particles)[1]
    i_N2 = findall(p -> p[2] == "N2", particles)[1]
    for i in 1:length(particles)
        if particles[i][2] == "e-"
            ni_prod[i] = (t) -> e_prod_f(e_prod_itp, t)
        end
        if particles[i][2] == "O+(4S)"
            ni_prod[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_O , :]*0.56 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
        end
        if particles[i][2] == "O2+"
            ni_prod[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_O2, :]      ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
        end
        if particles[i][2] == "N2+"
            ni_prod[i] = (t) -> e_prod_f(e_prod_itp, t) .* n0[i_N2, :]*0.92 ./(n0[i_O, :]*0.56 .+ n0[i_N2, :]*0.92 .+ n0[i_O2, :])
        end
    end
    return ni_prod
end
