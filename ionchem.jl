module ionchem

include("interpolate_temp.jl")
include("juliaIC.jl");

using DifferentialEquations
export ic

#load reactions, define particles etc.
const path_reactions_file = "test_data/Reaction rates full set ext.txt"
const dndt, particles, reactions, ode_raw, dndt_str, reactions_str = juliaIC.initIC(path_reactions_file)
rrates = [r[4] for r in reactions]

function myODEf(dn, n, p, t)
    rrates, nprod, dndt, temp_itp, temp_2, rr, X = p
    temp_2 = temp(temp_itp, t)
    rr = [r(temp_2, X) for r in rrates]
    for j in axes(n, 1)#1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, rr, 0, n, t)
    end
    nothing
end

function ic(tspan, n0, nprod, temp_itp, nh)
    temp_2 = temp(temp_itp, 0.1)
    X = ones(nh)
    rr = [r(temp_2, X) for r in rrates]

    prob = ODEProblem(myODEf, n0, tspan, (rrates, nprod, dndt, temp_itp, temp_2, rr, X))
    
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3);

    return sol
end

end