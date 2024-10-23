
using MAT

#load and return contents of ElSpec output
function loadmat(matfile)
    file = matopen(matfile)
    con = read(file, "ElSpecOut")
    close(file)
    return con
end

#get the relevant parameters for IC, some preprocessing applied
#rescale ion densities to restore plasma neutrality
#set start time to 0
#retrieve temperatures
function getparams(con)
    ne = con["ne"] #size: (62, 674)
    
    iri = permutedims(con["iri"], (2, 1, 3))
    Tn_, Ti_, Te_, nN2, nO2, nO, nAr, nNOp, nO2p, nOp = [iri[i, :, :] for i in 1:size(iri)[1]]
    Tn = Tn_
    nNOp, nO2p, nOp = [i .* ne ./ sum([nNOp, nO2p, nOp]) for i in [nNOp, nO2p, nOp]]
    nion = cat(ne, nN2, nO2, nO, nAr, nNOp, nO2p, nOp; dims = 3)
    nion = permutedims(nion, (3, 1, 2))

    par = con["par"]
    Ti = par[:, 2, :]
    Te = par[:, 3, :]
    
    temp = permutedims(cat(Te, Ti, Tn; dims = 3), (2, 1, 3)) #be aware of ordering!!!
    #size: (674, 62, 3)
    
    ts_ = dropdims(con["ts"]; dims = 2)
    te_ = dropdims(con["te"]; dims = 2)
    ts = ts_ .- ts_[1]
    te = te_ .- ts_[1]
    
    e_prod = permutedims(con["q"], (2, 1))
    #size 674x62
    
    h = dropdims(con["h"]; dims = 1)
    
    return ts, te, h, nion, temp, e_prod
end