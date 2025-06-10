using JLD2

function save_ic(path, tsol, ni, h, T, e_prod, particles, ts)
    jldsave(path; 
            tsol, 
            ni, 
            h, 
            T, 
            e_prod, 
            particles, 
            ts)
end

function load_ic(file)
    #@load file
    data = load(file)

    tsol        = data["tsol"]
    ni          = data["ni"]
    h           = data["h"]
    T           = data["T"]
    e_prod      = data["e_prod"]
    particles   = data["particles"]
    ts          = data["ts"]

    #assign_densities(ni, particles)
    return tsol, ni, h, T, e_prod, particles, ts
end

function assign_densities(nspecies, particles, prefix="")
    np = length(particles)
    for idx in 1:np
        name = particles[idx][2]
        name = replace(name,
                       "+" => "p", 
                       "(" => "_", 
                       ")" => "",
                       "-" => ""
                       )
        name = "n"*name
        name = prefix*name
        println("assigning "*name)
        
        global ni
        ni = nspecies
        eval(Meta.parse("global $name = transpose(ni[:, $idx, :])"))
    end
end