module ionchem

include("interpolate_temp.jl")
include("chemistry.jl");

using DifferentialEquations

#load reactions, define particles etc.
const path_reactions_file = "test_data/Reaction rates full set ext.txt"
const dndt, particles, reactions, ode_raw, dndt_str, reactions_str = chemistry.initIC(path_reactions_file)
const rrates = [r[4] for r in reactions]

#evaluate ODEs at a timestep
function myODEf(dn, n, p, t)
    rrates, ni_prod, dndt, temp_itp, T, rr, X = p
    T = temp_f(temp_itp, t)
    rr = [r(T, X) for r in rrates]
    for j in axes(n, 1)
        dn[j, :] .= dndt[j](ni_prod, rr, 0, n, t)
        #        return (nprodd, rr, tem, nn, tt) -> Base.invokelatest(f, nprodd, rr, tem, nn, tt)
    end
    nothing
end

function ic(tspan, n0, ni_prod, temp_itp, nh, t_save = [], t_cb = [], cb_f = [])
    T = temp_f(temp_itp, 0.1)
    X = ones(nh)
    rr = [r(T, X) for r in rrates]
    
    cb = PresetTimeCallback(t_cb, cb_f)

    prob = ODEProblem(myODEf, n0, tspan, (rrates, ni_prod, dndt, temp_itp, T, rr, X))
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save, callback = cb);
    return sol
end

end