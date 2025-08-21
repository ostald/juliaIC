module ionchem

include("loadElspec.jl")
include("interpolate_temp.jl")
include("ion_prod.jl")
include("chemistry.jl")
include("main.jl")

export loadmat, getparams
export interpolate_temp
export interpolate_q, assign_prod, e_prod_f
export solveIC, solveIC_allAtOnce, particles, reactions, initIC
export myODEf, ic

end # module ionchem
