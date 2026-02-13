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

function assign_densities(ni, particles)
    #produces a named touple ni_ntup from 3d matrix ni
    #so that it can be called ni_ntup.particle_name
    #organise data in vectors:
    pn = [p[2] for p in particles]
    pnames  = replace.(pn, 
                       "+" => "p", 
                       "(" => "_", 
                       ")" => "",
                       "-" => ""
                       )
    ni_vec = [ni[:, i, :] for i in axes(ni, 2)]
    #put in to dict
    d = Dict(pnames .=> ni_vec)
    #convert into named tuple
    ni_ntup = NamedTuple((Symbol(key),value) for (key,value) in d)
    return ni_ntup
end
