#!/usr/bin/env julia

function buildcomplex3(symmat::Matrix{Int}, maxsd::Int)
    grain = Array{Array{Int,1}}(undef,maxsd+1)
    farfaces = Array{Array{Int,1}}(undef,maxsd+1)
    prepairs = Array{Array{Int,1}}(undef,maxsd+1)
    firstv = Array{Array{Int,1}}(undef,maxsd+1)
    
    grain[maxsd+1] = Array{Int}(undef,0)
    farfaces[maxsd+1] = Array{Int}(undef,0)
    prepairs[maxsd+1] = Array{Int}(undef,0)
    firstv[maxsd+1] = [1]

    m = size(symmat,1)
    w = offdiagmean(symmat) |> vec

    vperm = sortperm(w, rev=true, alg=MergeSort)
    symmat = symmat[vperm,vperm]

    farfaces[1] = 1:m
    firstv[1] = 1:m+1
    grain[1] = diag(symmat)
    prepairs[1] = Array{Int}(undef,0)

    farfaces[2], firstv[2], grain[2] = generate2faces(symmat)
    prepairs[2] = Array{Int}(undef,0)

    if maxsd == 3
        generate3faces!(farfaces, firstv, grain, prepairs, m, symmat)
        return (farfaces=farfaces, firstv=firstv, grain=grain, prepairs=prepairs, symmat=symmat, nvl2ovl=vperm)
    end

    fpi = Int[]
    ff2pv = Int[]
    pmhist = zeros(Int, m, m)

    for sd = 3:maxsd
        nl = length(farfaces[sd-1])
        nll = length(farfaces[sd-2])
        
        npsupp = trues(nl)
        pflist = Array{Int}(undef,nl)
        jrv = farfaces[sd-1]
        jcp = firstv[sd-1]
        jz = grain[sd-1]
        zll = grain[sd-2]
        izfull = Array{Int}(undef,nll)
        r = Int[]
        z = Int[]
        c = Array{Int}(undef,m+1)
        c[1] = 1
        
        if sd == maxsd-1
            ff2pv = fill(m+1, nl)
        elseif sd == maxsd
            #### sort j-matrix by grain
            alterweight = 1 + maximum(zll) .- zll
            
            lowfilt = alterweight[jrv]
            invertiblevec = integersinsameorderbycolumn2(lowfilt, jcp)
            inversevec0 = Array{Int}(undef,nl)
            inversevec0[invertiblevec] = 1:nl
            jrv = jrv[inversevec0]
            jz = jz[inversevec0]

            lowfilt = ff2pv[jrv]
            invertiblevec = integersinsameorderbycolumn2(lowfilt, jcp)
            inversevec1 = Array{Int}(undef,nl)
            inversevec1[invertiblevec] = 1:nl
            jrv = jrv[inversevec1]
            jz = jz[inversevec1]
            
            translatorvecb = inversevec0[inversevec1]
            rt, ct, zt = transposeLighter(jrv, jcp, jz, nll)
            colsum = ct .- 1

            # for sloth (apologies) we'll leave some unsed stuff in row m+1
            pmhist = zeros(Int,m+1,m)
            fpi = zeros(Int,m+1,m)
            for p = 1:m
                for q = jcp[p]:jcp[p+1] - 1
                    pmhist[ff2pv[jrv[q]],p] += 1
                end
            end
            fpi[1,1:m] = jcp[1:m]
            for p = 1:m
                for q = 1:m
                    fpi[q+1,p] = fpi[q,p]+pmhist[q,p]
                end
            end

            #### reset ff2pv for next round
            ff2pv = fill(m+1, nl)

            oldclaw = Array{Int}(undef,m)
        end

        numpairs = 0
        
        for i = 1:m
            fill!(izfull, 0)
            lrange = cran(jcp, i)
            izfull[jrv[lrange]] = jz[lrange]

            for j = i+1:m
                dij = symmat[j,i]
                dij == 0 && continue
                if sd <= maxsd - 1
                    k = cran(jcp, j)
                    kk = jrv[k]
                    if1 = izfull[kk] .> 0
                    k = k[if1]
                    kk = kk[if1]
                    farfilt = jz[k]
                    claw = min.(izfull[kk], dij)
                    append!(r, k)
                    append!(z, min.(farfilt, claw))
                    if2 = (claw .>= farfilt) .& npsupp[k]
                    pflist[numpairs .+ (1:sum(if2))] .= (length(r)-length(k)+1:length(r))[if2]
                    npsupp[k[if2]] .= false
                    numpairs += sum(if2)
                    if sd == maxsd - 1
                        pmhist[i, j] += sum(if2)
                        ff2pv[k[if2]] .= i
                    end
                else
                    oldclaw[1:i-1] = minimum(symmat[1:i-1, [i, j]]; dims=2)
                    for l = 1:i-1
                        if fpi[l, j] < fpi[l+1, j]
                            ocl = oldclaw[l]
                            if ocl < dij
                                for k = fpi[l, j]:fpi[l+1, j]-1
                                    ## may have to reindex this
                                    kk = jrv[k]
                                    farfilt = jz[k]
                                    if zll[kk] <= ocl
                                        break
                                    elseif oldclaw[l] < min(farfilt, izfull[kk])
                                        claw = min(izfull[kk], dij)
                                        if claw >= farfilt
                                            if npsupp[k]
                                                push!(r, k)
                                                push!(z, farfilt)
                                                numpairs += 1
                                                pflist[numpairs] = length(r)
                                                npsupp[k] = false
                                                ff2pv[k] = i
                                            elseif oldclaw[ff2pv[k]] >= farfilt
                                                continue
                                            elseif saveface(ct, kk, colsum, farfilt, oldclaw, rt, zt)
                                                push!(r, k)
                                                push!(z, farfilt)
                                            end
                                        elseif (claw > 0) && saveface(ct, kk, colsum, claw, oldclaw, rt, zt)
                                            push!(r, k)
                                            push!(z, claw)
                                        end
                                    end
                                end
                            end
                        end
                    end
                    
                    k = fpi[i,j]:fpi[i+1,j]-1
                    kk = jrv[k]
                    farfilt = jz[k]
                    claw = min.(izfull[kk], dij)
                    if2 = (dij .>= farfilt) .& npsupp[k]
                    farfilt[.!if2] .= min.(claw[.!if2], farfilt[.!if2])
                    if3 = saveface.(Ref(ct), kk, Ref(colsum), farfilt, Ref(oldclaw), Ref(rt), Ref(zt))
                    pflist[numpairs .+ (1:sum(if2))] .= length(r) .+ cumsum(if2 .| if3)[if2];
                    numpairs += sum(if2)
                    append!(r, k[if2 .| if3])
                    append!(z, farfilt[if2 .| if3])
                    npsupp[k[if2]] .= false
                    ff2pv[k[if2]] .= i
                
                    k = fpi[i+1,j]:fpi[j,j]-1
                    kk = jrv[k]
                    if1 = izfull[kk] .> 0
                    k = k[if1]
                    kk = kk[if1]
                    farfilt = jz[k]
                    claw = min.(izfull[kk], dij)
                    if2 = (claw .>= farfilt) .& npsupp[k]
                    farfilt[.!if2] .= min.(claw[.!if2], farfilt[.!if2])
                    if3 = saveface.(Ref(ct), kk, Ref(colsum), farfilt, Ref(oldclaw), Ref(rt), Ref(zt))
                    pflist[numpairs .+ (1:sum(if2))] .= length(r) .+ cumsum(if2 .| if3)[if2];
                    numpairs += sum(if2)
                    append!(r, k[if2 .| if3])
                    append!(z, farfilt[if2 .| if3])
                    npsupp[k[if2]] .= false
                    ff2pv[k[if2]] .= i
                    
                    k = fpi[j,j]:fpi[j+1,j]-1
                    kk = jrv[k]
                    if1 = izfull[kk] .> 0
                    k = k[if1]
                    kk = kk[if1]
                    claw = min.(izfull[kk], dij)
                    if2 = saveface.(Ref(ct), kk, Ref(colsum), claw, Ref(oldclaw), Ref(rt), Ref(zt))
                    append!(r, k[if2])
                    append!(z, claw[if2])
                
                    k = fpi[j+1,j]:jcp[j+1]-1
                    kk = jrv[k]
                    if1 = izfull[kk] .> 0
                    k = k[if1]
                    kk = kk[if1]
                    farfilt = jz[k]
                    claw = min.(izfull[kk], dij)
                    if2 = (claw .>= farfilt) .& npsupp[k]
                    farfilt[.!if2] .= min.(claw[.!if2], farfilt[.!if2])
                    if3 = saveface.(Ref(ct), kk, Ref(colsum), farfilt, Ref(oldclaw), Ref(rt), Ref(zt))
                    pflist[numpairs .+ (1:sum(if2))] .= length(r) .+ cumsum(if2 .| if3)[if2];
                    numpairs += sum(if2)
                    append!(r, k[if2 .| if3])
                    append!(z, farfilt[if2 .| if3])
                    npsupp[k[if2]] .= false
                    ff2pv[k[if2]] .= i
                end
            end
            # update the column pattern and the total number of nonzeros
            # encountered per codim2 face
            c[i+1] = length(r) + 1
            if sd == maxsd
                colsum[jrv[cran(jcp,i)]] .+= 1
            end
        end

        deleteat!(pflist, numpairs+1:nl)
        if sd == maxsd r = translatorvecb[r] end
        firstv[sd] = c
        farfaces[sd] = r
        prepairs[sd] = pflist
        grain[sd] = z
        if isempty(farfaces[sd])
            for nextcard = sd+1:maxsd
                firstv[nextcard] = [1;1]
                farfaces[nextcard] = Array{Int}(undef,0)
                prepairs[nextcard] = Array{Int}(undef,0)
                grain[nextcard] = Array{Int}(undef,0)
            end
            break
        end
    end
    (farfaces=farfaces, firstv=firstv, grain=grain, prepairs=prepairs, symmat=symmat, nvl2ovl=vperm)
end


function saveface(ct::Vector{Int},kk::Int,colsum::Vector{Int},farfilt::Int,oldclaw::Vector{Int},rt::Vector{Int},zt::Vector{Int})::Bool
    for l = ct[kk]:colsum[kk]
        if (zt[l] >= farfilt) && (oldclaw[rt[l]] >= farfilt)
            return false
        end
    end
    true
end

function generate2faces(symmat::Matrix{Int})
    m = size(symmat,1)
    L = sum(tril(symmat, -1) .> 0)
    rowval = Array{Int}(undef,L)
    nzval = Array{Int}(undef,L)
    colptr = Array{Int}(undef,m+1)
    marker = 0
    colptr[1] = 1
    for i = 1:m
        colptr[i+1]=colptr[i]
        for j = (i+1):m
            if symmat[j,i] > 0
                colptr[i+1] += 1
                rowval[colptr[i+1]-1] = j
                nzval[colptr[i+1]-1] = symmat[j,i]
            end
        end
    end
    rowval,colptr,nzval
end

function generate3faces!(farfaces_cell, firstv_cell, grain_cell, prepairs_cell, m, symmat)
    grain = grain_cell[2]::Vector{Int}
    farfaces = farfaces_cell[2]::Vector{Int}
    firstv = firstv_cell[2]::Vector{Int}

    numverts = length(firstv)-1
    numedges = length(farfaces)
    numpairs = 0

    closefaces = Array{Int}(undef,numedges)
    for i = 1:m closefaces[cran(firstv,i)] .= i end
    iso = integersinsameorder(farfaces)
    closefaces_higsorted = Array{Int}(undef,numedges)
    grain_higsorted = Array{Int}(undef,numedges)
    closefaces_higsorted[iso] = closefaces
    grain_higsorted[iso] = grain

    firstv_hs = zeros(Int,m+1)
    for face in farfaces
        firstv_hs[face+1] += 1
    end
    firstv_hs[1] = 1
    for i = 2:m+1
        firstv_hs[i] = firstv_hs[i-1]+firstv_hs[i]
    end

    adist = Array{Int}(undef,m)
    idist = Array{Int}(undef,m)
    r = Array{Int}(undef,numedges)
    z = Array{Int}(undef,numedges)
    s = Array{Int}(undef,numedges)

    clawvec = Array{Int}(undef,m)
    ncheckedges = trues(numedges)

    for a = 1:m
        adist .= 0
        adist[crows(firstv,farfaces,a)] = crows(firstv,grain,a)
        for ip = cran(firstv,a)
            i = farfaces[ip]
            dai = grain[ip]
            idist .= 0
            idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
            idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
            for jp = cran(firstv,i)
                if ncheckedges[jp]
                    j = farfaces[jp]
                    dij = grain[jp]
                    # note this condition bakes in the req. that j be adjacent to a
                    if dij <= dai && dij <= adist[j]
                        numpairs += 1
                        ncheckedges[jp] = false
                        clawvec[1:i] .= 0
                        for lp = cran(firstv_hs,j)
                            l = closefaces_higsorted[lp]
                            if l >= i
                                break
                            elseif idist[l] != 0
                                clawvec[l] = min(idist[l],grain_higsorted[lp])
                            end
                        end
                        for kp = cran(firstv,j)
                            k = farfaces[kp]
                            djk = grain[kp]
                            dak = adist[k]
                            dik = idist[k]
                            # this bakes in req. that dik > 0
                            if dak < dij && dak<djk && dak < dik
                                dijk = min(dij,dik,djk)
                                keepface = true
                                for bp = cran(firstv_hs,k)
                                    b = closefaces_higsorted[bp]
                                    if b >= i
                                        break
                                    elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
                                        keepface = false
                                        break
                                    end
                                end
                                if keepface
                                    push!(r, kp)
                                    push!(z, dijk)
                                    push!(s, i)
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    holdi = 0
    for edge = findall(ncheckedges)
        i = closefaces[edge]
        j = farfaces[edge]
        dij = grain[edge]
        if i != holdi
            fill!(idist, 0)
            idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
            idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
            holdi = i
        end
        clawvec[1:i] .= 0
        for lp = cran(firstv_hs,j)
            l = closefaces_higsorted[lp]
            if l >= i
                break
            elseif idist[l] != 0
                clawvec[l] = min(idist[l],grain_higsorted[lp])
            end
        end
        #### a facsimile of above
        for kp = cran(firstv,j)
            k = farfaces[kp]
            dik = idist[k]
            if dik==0
                continue
            end
            djk = grain[kp]
            dijk = min(dij,dik,djk)
            keepface = true
            for bp = cran(firstv_hs,k)
                b = closefaces_higsorted[bp]
                if b >= i
                    break
                elseif min(clawvec[b],grain_higsorted[bp]) >= dijk
                    keepface = false
                    break
                end
            end
            if keepface
                push!(r, kp)
                push!(z, dijk)
                push!(s, i)
            end
        end
        ####
    end

    num3faces = length(r)
    holderlengths = length(r)
    deletionrange = (num3faces+1):holderlengths
    deleteat!(r,deletionrange)
    deleteat!(z,deletionrange)
    deleteat!(s,deletionrange)

    iso = integersinsameorder(s)
    r[iso] = r
    z[iso] = z
    fv3 = zeros(Int,numverts+1)
    fv3[1] = 1
    for face in s
        fv3[face+1]+=1
    end
    for i = 2:(numverts+1)
        fv3[i] = fv3[i-1]+fv3[i]
    end

    pairmarker = 0
    npes = trues(numedges)
    prepairs = Array{Int}(undef,numedges)
    for i = 1:num3faces
        edge = r[i]
        if npes[edge] && z[i] == grain[edge]
            npes[edge]=false
            pairmarker+=1
            prepairs[pairmarker]=i
        end
    end
    deleteat!(prepairs,(pairmarker+1):numedges)

    farfaces_cell[3]=r
    firstv_cell[3]=fv3
    grain_cell[3]=z
    prepairs_cell[3]=prepairs

    buffer1 = Array{Array{Int,1},1}(undef,1)
    buffer2 = Array{Array{Int,1},1}(undef,1)
    buffer1[1] = Array{Int}(undef,0)
    buffer2[1] = ones(Int,numverts+1)
    append!(farfaces_cell,buffer1)
    append!(grain_cell,buffer1)
    append!(prepairs_cell,buffer1)
    append!(firstv_cell,buffer1)

    return r,fv3,z,prepairs,numpairs
end

