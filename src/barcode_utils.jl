#!/usr/bin/env julia

function nnzbars(grain, plo, phi, tid; dim::Int=1)
    sd = dim+2
    @assert 0 < sd <= length(plo) "requested dimension is outside the range of calculated bars"
    plo = plo[sd]
    phi = phi[sd]
    lowfilt = grain[sd-1]
    higfilt = grain[sd]
    counter = 0
    for i = 1:length(plo)
        if higfilt[phi[i]] != lowfilt[plo[i]]
            counter += 1
        end
    end
    counter + length(tid[sd])-length(plo)
end

barname2cyclename(grain, plo, phi, tid; dim::Int=1) = barname2cyclename(grain, plo, phi, tid, [1]; dim=dim)
function barname2cyclename(grain, plo, phi, tid, barnumber; dim::Int=1)
    sd = dim+2
    tid = tid[sd]
    plo = plo[sd]
    phi = phi[sd]
    nummortals = length(plo)
    nzcycles = findall(grain[sd][phi] .!= grain[sd-1][plo])
    append!(nzcycles, nummortals+1:length(tid))
    return nzcycles[barnumber]
end
function barname2cyclename(grain, plo, phi, tid, barnumber::Int; dim::Int=1)
    sd = dim+2
    tid = tid[sd]
    plo = plo[sd]
    phi = phi[sd]
    numclasses = length(tid)
    nummortals = length(plo)
    counter = 0
    cyclename = 0
    for i = 1:nummortals
        if grain[sd][phi[i]] != grain[sd-1][plo[i]]
            counter+=1
        end
        if counter == barnumber
            cyclename = i
            break
        end
    end
    if cyclename == 0
        for i = (nummortals+1):length(tid)
            counter+=1
            if counter == barnumber
                cyclename = i
                break
            end
        end
    end
    return cyclename
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    if sd == 2
        # this case must be treated separately b/c julia indexing starts at 1
        numlowlows = 0
        numnlpll = 0
    else
        numlowlows = length(cp[sd-2])-1
        numnlpll = numlowlows-length(plo[sd-2])
    end
    # the -1 accounts for the fact that length(cp[sd-1]) = (# cells of dimension secard-1) - 1
    numnlpl = length(cp[sd-1])-1-length(plo[sd-1])

    summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
    append!(summands,[tid[sd][cyclenumber]])

    sd == 2 && return summands

    supp = falses(numlowlows)
    sc = sd-1
    # this is a bit ridiculus; it's here bc i don't have a verified form of spmmF2 on hand
    for j in summands
        for k in cran(cp[sc],j)
            i = rv[sc][k]
            supp[i] = !supp[i]
        end
    end
    brv = findall(supp[tid[sd-1]])
    bcp = [1,length(brv)+1]
    brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpll)
    brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpll)
    plow2phigtranslator = zeros(Int,numlowlows)
    plow2phigtranslator[plo[sd-1]]=phi[sd-1]
    brv = plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
    append!(brv,summands)
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Array{Int,1})
    numclasses = length(cyclenumber)
    rep  = Array{Array{Int,1},1}(undef,numclasses)
    for p = 1:numclasses
        rep[p] = getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber[p])
    end
    rep
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    numlowlows = sd == 2 ? 0 : length(farfaces[sd-2])
    numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

    summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
    append!(summands,[tid[sd][cyclenumber]])
    brv = ff2aflight(farfaces,firstv,sd-1,summands)
    supp = falses(length(farfaces[sd-2]))
    supp[brv] .= .!supp[brv]

    brv = findall(supp[tid[sd-1]])
    bcp = [1,length(brv)+1]
    brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
    brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

    plow2phigtranslator = zeros(Int,numlowlows)
    plow2phigtranslator[plo[sd-1]]=phi[sd-1]
    # some of the nonzero entries might lie in non-basis rows
    brv = plow2phigtranslator[tid[sd-1][brv]]
    brv = brv[findall(!iszero,brv)]
    append!(brv, summands)
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    numlowlows = sd == 2 ? 0 : length(farfaces[sd-2])
    numlows    = length(farfaces[sd-1])
    numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

    numclasses = length(cyclenumber)
    summands = Array{Array{Int,1},1}(undef,numclasses)
    rep  = Array{Array{Int,1},1}(undef,numclasses)
    summandsupp = falses(numlows)
    for i = 1:numclasses
        summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
        append!(summands[i],[tid[sd][cyclenumber[i]]])
        summandsupp[summands[i]] .= true
    end

    lowgenerators = findall(summandsupp)
    numlowgenerators = length(lowgenerators)
    translator = zeros(Int,numlows)
    translator[lowgenerators] = 1:length(lowgenerators)

    lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

    supp = falses(numlowlows)
    m = size(lowfacemat,1)
    plow2phigtranslator = Array{Int}(undef,numlowlows)
    plow2phigtranslator[plo[sd-1]] = phi[sd-1]
    for i = 1:numclasses

        fill!(supp, false)
        for j = 1:length(summands[i])
            for k = 1:m
                kk = lowfacemat[k,translator[summands[i][j]]]
                supp[kk] = !supp[kk]
            end
        end

        brv = findall(supp[tid[sd-1]])
        bcp = [1,length(brv)+1]
        brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numlowlows)
        brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numlowlows)

        rep[i] = append!(plow2phigtranslator[tid[sd-1][brv]],summands[i])
    end
    rep
end





