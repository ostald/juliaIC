using PCHIPInterpolation

function interpolate_q(t, e_prod)
    e_prod_itp = [Interpolator(t, e_prod[:, ih]) for ih in axes(e_prod, 2)]
    return e_prod_itp
end

function e_prod_f(e_prod_itp, t)
    return (i.(t) for i in e_prod_itp)
end