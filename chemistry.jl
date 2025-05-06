module chemistry

export solveIC, solveIC_allAtOnce, particles, reactions, initIC
#other variables or methods can be accessed by typing modulName.method

using Profile
using SciPyDiffEq
#using Plots


function readreactionfile(rfilepath)
    #reads reactions from a file, returns strings
    #removes tabs and splits lines by semicolons
    open(rfilepath, "r") do rfile
        recordLine = false
        reactions = ""
        while !eof(rfile)
            s = readline(rfile)
            if first(s, 12) == "--Reactions:" 
                recordLine = true
                continue
            end
            if s == ""
                recordLine = false
            end
            if recordLine == true
                reactions = reactions *"\n"* s
            end
        end
        reactions1 = chop(reactions, head = 1)
        reactions2 = replace(reactions1, "\t" => "")
        reactions_str = [split(line, ";") for line in split(reactions2, "\n")]
        return reactions_str
    end
end

function getReactionsParticles(reactions_str)
    # parses single reaction info from string into list of particles, and reactions with info in vecotr
    reactions_ = []
    particles_ = []
    for (i, r) in enumerate(reactions_str)
        edus, pros = split(r[2], "=>")
        educts = [replace(edu, ' ' => "") for edu in split(edus, " + ")]
        nproducts = [replace(pro, ' ' => "") for pro in split(pros, " + ")]
        for edu in educts
            if !(edu in particles_)
                append!(particles_, [edu])
            end
        end
        for pro in nproducts
            if !(pro in particles_)
                append!(particles_, [pro])
            end
        end
        reactionrate = strip(r[3], ' ')
        branchingratio = [replace(ratio, ' ' => "") for ratio in split(r[4], ',')]
        reaction_info = [i, r[1], educts, nproducts, reactionrate, branchingratio]            
        append!(reactions_, [reaction_info])
    end
    return particles_, reactions_
end

function orderParticles(particles_, ordering = [])
    if ordering != []
        #stupid ordering algorithm, no idea why it needs to be so complicated
        particles_2 = []
        for o in ordering
            for p in particles_
                if o == replace(p, "-" => "")
                    #println(p)
                    append!(particles_2, [p])
                end
            end
        end
    else 
        particles_2 = particles_
    end
    particles = [[i, p] for (i, p) in enumerate(particles_2)]
    return particles
end


function parseReactionRate(s)
    f = eval(Meta.parse("(Te, Tn, Ti, Tr, X) -> X .* " * s))
    return (T, X) -> Base.invokelatest(f, T[1], T[2], T[3], (T[2]+T[3])/2, X)
end


function setReactions(reactions_, particles)
    #put reactions in a computer friendly format
    #all string to numbers where sensible
    #replace particle names with index number
    #replace reaction string with actual function
    reactions = []
    reactions_str = []
    for r in reactions_
        reaction_number = r[1]
        educts_indices = [findfirst(x -> x[2] == ed, particles) for ed in r[3]]
        #println(ed_ind)
        nproducts_indices = [findfirst(x -> x[2] == pr, particles) for pr in r[4]]
        #println(pr_ind)
        reaction_rate_str = replace(
                replace(
                    replace(
                        replace(
                            replace(
                                replace(
                                    replace(
                                        replace(r[5], "<=" => ".<=")
                                        , "/" => " ./")
                                    , "**" => ".^")
                                , "*" => ".*")
                            , "np." => "")
                        , "where" => "ifelse.")
                    , "exp" => "exp.")
                , "sqrt" => "sqrt.")
        
        #println(reaction_rate_str)
        reaction_rate_function = parseReactionRate(reaction_rate_str)
        branchingratio_float = [if br == "" 1 else parse(Float64, br) end for br in r[6]]
        reaction_info = [reaction_number, educts_indices, nproducts_indices, reaction_rate_function, branchingratio_float]
        append!(reactions, [reaction_info])
        append!(reactions_str, [reaction_rate_str])
    end
    return reactions, reactions_str
end

function get_ode_mat(reactions, particles)
    ode_mat = Array{Union{Float64, Vector}}(undef, length(reactions), length(particles))
    fill!(ode_mat, [0])

    for re in reactions
        ir = re[1] #index
        ed = re[2]
        pr = re[3]
        br = re[5]
        for ed_i in ed
            ode_mat[ir, ed_i] = [ir, -1.0, ed]
        end
        if br == [1]
            br = ones(length(pr))
        end
        for (pr_i, bra) in zip(pr, br)
            ode_mat[ir, pr_i] = [ir, 1.0.*bra, ed]
        end
        for ed_i in ed
            for pr_i in pr
                if ed_i == pr_i
                    ode_mat[ir, ed_i] = [0] #this only works if branching ratio = 1
                end
            end
        end
    end
    return ode_mat
end

function get_ode_raw(ode_mat, particles)
    ode_raw = Array{Vector}(undef, length(particles))
    for i in 1:length(particles)
        ode_raw[i] = [info for info in ode_mat[:, i]]
        if particles[i][2] in ["H", "O", "O2", "N2"]
            ode_raw[i] = [0]
        end
    end
    return ode_raw
end



function initIC(path_reactions_file, ordering = [])
    
    reactions_str = readreactionfile(path_reactions_file)
    particles_, reactions_ = getReactionsParticles(reactions_str)
    particles = orderParticles(particles_, ordering)
    reactions, reactions_str = setReactions(reactions_, particles)
    ode_mat = get_ode_mat(reactions, particles)
    ode_raw = get_ode_raw(ode_mat, particles)

    #start compiling ODE:
    function dndtFromString(s)
        f = eval(Meta.parse("(nprodd, rr, tem, nn, tt) -> " * s))
        return (nprodd, rr, tem, nn, tt) -> Base.invokelatest(f, nprodd, rr, tem, nn, tt)
    end

    #function dndtFromString(s)
    #    f = eval(Meta.parse("(t) -> " * s))
    #    return (t) -> Base.invokelatest(f,t)
    #end

    dndt = Array{Any}(undef, length(particles))
    dndt_str = Array{Any}(undef, length(particles))

    for i in 1:length(particles)
        #dndt_str = ""
        if ode_raw[i] == [0]
            #dndt_str[i] = "zeros(size(tem[1]))"
            dndt_str[i] = "zeros(length(tem[1]))"
            #print(particles[i][2], "\n") #all species that are assumed of constant density (N2, O2, O, H)
        else
            #print(ode_raw[i], "\n")
            dndt_str[i] = "nprodd[$i]"
            for o in ode_raw[i]
                if o == [0] #if not participating in reaction, no terms are added
                else
                    #print(o, "\n")
                    dndt_str[i] = dndt_str[i] * " .+ $(o[2]) .* rr[$(o[1])]"
                    #dndt_str[i] = dndt_str[i] * " \n .+ $(o[2]).* $(reactions_[o[1]][2])" 
                    for ed_i in o[3]
                        #dndt_str[i] = dndt_str[i] * " .*  @view(nn[:, $ed_i])" #for transposed version
                        dndt_str[i] = dndt_str[i] * " .*  @view(nn[$ed_i, :])"
                        #dndt_str[i] = dndt_str[i] * " $(particles[ed_i][2])"
                    end
                end
            end
        end
        #println(dndt_str[i], "\n")
    end

    for i in 1:length(particles)
        dndt[i] = dndtFromString(dndt_str[i])
    end

    return dndt, particles, reactions, ode_raw, dndt_str, reactions_str
end



end


"""


function dummyf2(n, p, t)
    dn = zeros(size(n))
    reactions, nprod, dndt, temp, ih = p
    temp_ = temp(t)
    temp_2 = [temp_[i, :] for i in 1:size(temp_,1)]
    
    for j in 1:size(n)[1]
        #dn[j, :] .= nprod[j](t)
        #println(size(dndt[j](nprod, reactions, temp_2, n, t)))
        dn[j] = dndt[j](nprod, reactions, temp_2, n, t)[ih]
    end 
    return dn
end

function dummyf(n, p, t)
    dn = zeros(size(n))
    reactions, nprod, dndt, temp = p
    temp_ = temp(t)
    temp_2 = [temp_[i, :] for i in 1:size(temp_,1)]
            
    for j in 1:size(n)[1]
        dn[j, :] .= dndt[j](nprod, reactions, temp_2, n, t)
    end 
    return dn
end


function solveIC(n0, ts, te, nprod, temp, t_save, reactions, dndt);
    
    sol = Array{Any}(undef, 62)

    for i in 1:62
        #println(i)
        
        u0 = n0[:, i]
        tspan = (ts, te)
        prob = ODEProblem(dummyf2, u0, tspan, (reactions, nprod, dndt, temp, i))
        sol[i] = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save)
    end
    return sol
    
end



function solveIC_allAtOnce(n0, ts, te, nprod, temp, t_save, reactions, dndt)
    
    tspan = (ts, te)
    prob = ODEProblem(dummyf, n0, tspan, (reactions, nprod, dndt, temp))
    sol = solve(prob, TRBDF2(autodiff=false), reltol = 1e-7, abstol = 1e-3, saveat = t_save)

    return sol
    
end

end
"""