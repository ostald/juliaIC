module ionchem

include("interpolate_temp.jl")
include("ion_prod.jl")
include("chemistry.jl")
include("main.jl")
include("ic_io.jl")
include("get_msis.jl")

export msis
export save_ic, load_ic, assign_densities
export interpolate_temp
export interpolate_q, assign_prod, e_prod_f
export solveIC, solveIC_allAtOnce, particles, reactions, initIC
export myODEf, ic

end # module ionchem
