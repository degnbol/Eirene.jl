#!/usr/bin/env julia

function nnzbars(D::Dict;dim = 1)
    sd = dim+2
    if !(0 < sd <= length(D["plo"]))
        print("error: requested dimension is outside the range of calculated bars")
        return
    end
    plo = D["plo"][sd]
    phi = D["phi"][sd]
    lowfilt = D["grain"][sd-1]
    higfilt = D["grain"][sd]
    counter::Int = 0
    for i = 1:length(plo)
        if higfilt[phi[i]] != lowfilt[plo[i]]
            counter+=1
        end
    end
    counter += length(D["tid"][sd])-length(plo)
    return counter
end

function barname2cyclename(D::Dict,barnumber = [1];dim = 1)
    if typeof(barnumber) <: Array
        sd = dim+2
        tid = D["tid"][sd]
        plo = D["plo"][sd]
        phi = D["phi"][sd]
        nummortals = length(plo)
        nzcycles = findall(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
        append!(nzcycles,Array((nummortals+1):length(tid)))
        return nzcycles[barnumber]
    elseif typeof(barnumber) <: UnitRange
        sd = dim+2
        tid = D["tid"][sd]
        plo = D["plo"][sd]
        phi = D["phi"][sd]
        nummortals = length(plo)
        nzcycles = findall(D["grain"][sd][phi].!= D["grain"][sd-1][plo])
        append!(nzcycles,nummortals+1:length(tid))
        return nzcycles[barnumber]
    elseif typeof(barnumber)<:Number
        sd = dim+2
        tid = D["tid"][sd]
        plo = D["plo"][sd]
        phi = D["phi"][sd]
        numclasses = length(tid)
        nummortals = length(plo)
        counter = 0
        cyclename = 0
        for i = 1:nummortals
            if D["grain"][sd][phi[i]] != D["grain"][sd-1][plo[i]]
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
end


function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
    if sd == 2
        numlowlows = 0
    else
        numlowlows = length(farfaces[sd-2])
    end
    numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

    summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
    append!(summands,[tid[sd][cyclenumber]])
    brv = ff2aflight(farfaces,firstv,sd-1,summands)
    supp = falses(length(farfaces[sd-2]))
    for k in brv
        supp[k] = !supp[k]
    end

    brv = findall(supp[tid[sd-1]])
    bcp = [1,length(brv)+1]
    brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
    brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

    plow2phigtranslator = zeros(Int,numlowlows)
    plow2phigtranslator[plo[sd-1]]=phi[sd-1]
    brv		= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
    brv = brv[findall(!iszero,brv)]

    brv = append!(brv,summands)

    return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Number)
    if sd == 2
        # this case must be treated separately b/c julia indexing starts at 1
        numlowlows = 0
        numnlpll = 0
    else
        numlowlows = length(cp[sd-2])-1
        numnlpll = numlowlows-length(plo[sd-2])
    end
    numnlpl = length(cp[sd-1])-1-length(plo[sd-1])	# the -1 accounts for the fact that length(cp[sd-1]) = (# cells of dimension secard-1) - 1

    summands = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber)]
    append!(summands,[tid[sd][cyclenumber]])

    if sd == 2
        return summands
    end

    supp = falses(numlowlows)
    sc = sd-1
    for j in summands	# this is a bit ridiculus; it's here bc i don't have a verified form of spmmF2 on hand
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
    brv					= plow2phigtranslator[tid[sd-1][brv]] # some of the nonzero entries might lie in non-basis rows
    brv = append!(brv,summands)
    return brv
end

function getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber::Array{Int,1})
    numclasses = length(cyclenumber)
    rep  = Array{Array{Int,1},1}(undef,numclasses)
    for p = 1:numclasses
        rep[p] = getcycle_cell(rv,cp,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber[p])
    end
    return rep
end

function getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    if sd == 2
        numlowlows = 0
    else
        numlowlows = length(farfaces[sd-2])
    end
    numlows    = length(farfaces[sd-1])
    numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

    numclasses = length(cyclenumber)
    summands = Array{Array{Int,1},1}(undef,numclasses)
    rep  = Array{Array{Int,1},1}(undef,numclasses)
    summandsupp = falses(numlows)
    for i = 1:numclasses
        summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
        append!(summands[i],[tid[sd][cyclenumber[i]]])
        summandsupp[summands[i]].=true
    end

    lowgenerators = findall(summandsupp)
    numlowgenerators = length(lowgenerators)
    translator = zeros(Int,numlows)
    translator[lowgenerators] = 1:length(lowgenerators)

    lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

    supp = falses(numlowlows)
    m = size(lowfacemat,1)
    plow2phigtranslator = Array{Int}(undef,numlowlows)
    plow2phigtranslator[plo[sd-1]]=phi[sd-1]
    for i = 1:numclasses

        supp[:] .= false
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
    return rep
end

function getcycle(D::Dict,sd,cyclenumber)
    if !haskey(D,"Lirv")
        if !haskey(D,"cyclerep")
            println("This object does not store a complete cycle basis.")
        else
            println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
        end
        return
    end
    farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
    Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
    rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    return rrv
end

function getcycle(D::Dict,cyclenumber;dim = 1)
    if !haskey(D,"Lirv")
        println("Error: the function <getcycle(D::Dict,cyclenumber;dim = 1)> assumes a key value for \"Lirv\" in the input object D.  This key value is absent.")
        return
    end
    sd = dim+2
    Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
    Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];
    plo=D["plo"];phi=D["phi"];tid=D["tid"]
    if haskey(D,"farfaces")
        farfaces = D["farfaces"];firstv = D["firstv"];
        rrv = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    else
        rrv = getcycle_cell(D["rv"],D["cp"],Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    end
    return rrv
end

function getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    if sd == 2
        numlowlows = 0
    else
        numlowlows = length(farfaces[sd-2])
    end
    numlows    = length(farfaces[sd-1])
    numnlpl = length(farfaces[sd-1])-length(plo[sd-1])

    numclasses = length(cyclenumber)
    summands = Array{Array{Int,1},1}(undef,numclasses)
    rep  = Array{Int}(undef,numclasses)
    summandsupp = falses(numlows)
    for i = 1:numclasses
        summands[i] = tid[sd][crows(Licp[sd],Lirv[sd],cyclenumber[i])]
        append!(summands[i],[tid[sd][cyclenumber[i]]])
        summandsupp[summands[i]]=true
    end

    lowgenerators = findall(summandsupp)
    numlowgenerators = length(lowgenerators)
    translator = zeros(Int,numlows)
    translator[lowgenerators] = 1:length(lowgenerators)

    lowfacemat = ff2aflight(farfaces,firstv,sd-1,lowgenerators)

    supp = falses(numlowlows)
    m = size(lowfacemat,1)
    plow2phigtranslator = Array{Int}(undef,numlowlows)
    plow2phigtranslator[plo[sd-1]]=phi[sd-1]
    for i = 1:numclasses

        supp[:].= false
        for j = 1:length(summands[i])
            for k = 1:m
                kk = lowfacemat[k,translator[summands[i][j]]]
                supp[kk] = !supp[kk]
            end
        end

        brv = findall(supp[tid[sd-1]])
        bcp = [1,length(brv)+1]
        brv,bcp = spmmF2silentLeft(Lrv[sd-1],Lcp[sd-1],brv,bcp,numnlpl)
        brv,bcp = spmmF2silentLeft(Rrv[sd-1],Rcp[sd-1],brv,bcp,numnlpl)

        rep[i] = length(brv)+length(summands[i])
    end
    return rep
end

function getcyclesize(D::Dict, cyclenumber; dim=1)
    sd = dim+2
    if !haskey(D,"Lirv")
        if !haskey(D,"cyclerep")
            println("This object does not store a complete cycle basis.")
        else
            println("This object does not store a complete cycle basis, only those cycles that represent persistent homology classes.")
        end
        return
    end
    farfaces = D["farfaces"];firstv = D["firstv"];Lirv = D["Lirv"];Licp=D["Licp"];Lrv=D["Lrv"]
    Lcp=D["Lcp"];Rrv=D["Rrv"];Rcp=D["Rcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]

    rrv = getcyclesize(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenumber)
    return rrv
end


function complexrank(C; dim=1)
    sd = dim+1
    if dim > C["input"]["maxdim"]+1 || dim < 0
        return 0
    elseif C["input"]["model"] == "complex"
        return length(C["cp"][sd])-1
    else
        return length(C["farfaces"][sd])
    end
end

function boundaryrank(C; dim=1)
    sd = dim+1;
    complexrank(C, dim=dim) == 0 ? 0 : length(C["plo"][sd])
end

function boundarycorank(C; dim=1)
    sd = dim+1;
    complexrank(C, dim=dim) == 0 ? 0 : complexrank(C, dim=dim)-boundaryrank(C, dim=dim)
end

empteval(f, a, c) = isempty(a) ? c : f(a)

