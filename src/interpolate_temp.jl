using PCHIPInterpolation

function interpolate_temp(t, temp)
    #interpolates a set of temperature measurements at times t
    # and return a function continious in time
    # temp must have time as first axis
    temp_itp = Array{Any}(undef, size(temp)[2:3])
    for i in axes(temp, 2)
        for j in axes(temp, 3)
            temp_itp[i, j] = Interpolator(t, temp[:, i, j])
        end
    end
    return temp_itp
end

#old, can be deprecated? i guess so.
function temp_matrix(temp_itp, t)
    #evaluates interpolated temperature temp_int at time t
    # returns matrix of size [3, 62]
    return [temp_itp[i, j](t) for i in axes(temp_itp, 1), j in axes(temp_itp, 2)]
end

function temp_f(temp_itp, t)
    #evaluates interpolated temperature temp_int at time t
    # returns vector of size 3 with vecotrs of size 62 
    # [3[62]]
    return [[temp_itp[i, j](t) for i in axes(temp_itp, 1)] for j in axes(temp_itp, 2)]
end