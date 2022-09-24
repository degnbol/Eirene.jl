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

    fpi = Array{Int}(undef,0)
    ff2pv = Array{Int}(undef,0)
    pmhist = zeros(Int,m,m)

    for sd = 3:maxsd
        nl = length(farfaces[sd-1])
        nll = length(farfaces[sd-2])

        startlength = nl
        stepsize = min(10^7, ceil(Int, nl/4))

        npsupp = trues(nl)
        pflist = Array{Int}(undef,nl)
        jrv = farfaces[sd-1]
        jcp = firstv[sd-1]
        jz = grain[sd-1]
        zll = grain[sd-2]
        izfull = Array{Int}(undef,nll)
        r = Array{Int}(undef,startlength)
        z = Array{Int}(undef,startlength)
        c = Array{Int}(undef,m+1)
        c[1] = 1
        numpairs = [0]
        facecount = [0]
        if sd == maxsd-1
            ff2pv = Array{Int}(undef,nl)
            ff2pv .= m+1
        end
        if sd == maxsd
            #### sort j-matrix by grain
            alterweight = Array{Int}(undef,length(zll));
            maxweight = maximum(zll);
            for i = 1:length(alterweight)
                alterweight[i] = 1+maxweight-zll[i]
            end
            lowfilt = [alterweight[v] for v in jrv]
            invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
            inversevec0 = Array{Int}(undef,nl)
            inversevec0[invertiblevec]=1:nl
            jrv = [jrv[v] for v in inversevec0]
            jz = [jz[v] for v in inversevec0]

            lowfilt = [ff2pv[v] for v in jrv]
            invertiblevec = integersinsameorderbycolumn2(lowfilt,jcp)
            inversevec1 = Array{Int}(undef,nl)
            inversevec1[invertiblevec]=1:nl
            jrv = [jrv[v] for v in inversevec1]
            jz = [jz[v] for v in inversevec1]
            translatorvecb = [inversevec0[v] for v in inversevec1]
            inversevec0 = []; inversevec1 = []; lowfilt = []; invertiblevec = []
            #gc()
            rt, ct, zt = transposeLighter(jrv, jcp, jz, nll)
            colsum = ct .- 1

            # for sloth (apologies) we'll leave some unsed stuff in row m+1
            pmhist = zeros(Int,m+1,m)
            fpi = zeros(Int,m+1,m)
            processfpi!(pmhist,fpi,jcp,jrv,ff2pv,m)

            #### reset ff2pv for next round
            ff2pvold = copy(ff2pv)
            ff2pv = fill(m+1, nl)

            oldclaw = Array{Int}(undef,m)
        end

        for i = 1:m
            izfull .= 0
            lrange = cran(jcp,i)
            izfull[jrv[lrange]] = jz[lrange]

            for j = (i+1):m
                dij = symmat[j,i]
                dij == 0 && continue
                if sd < maxsd-1
                    process_sd_lt_maxsd!(i,j,dij,stepsize,
                                         facecount,numpairs,
                                         jrv,jcp,jz,
                                         r,z,pflist,
                                         izfull,ff2pv,pmhist,
                                         npsupp)
                elseif sd == maxsd-1
                    process_sd_onelt_maxsd_1!(i,j,dij,stepsize,
                                              facecount,numpairs,
                                              jrv,jcp,jz,
                                              r,z,pflist,
                                              izfull,ff2pv,pmhist,
                                              npsupp)
                else
                    for l = 1:(i-1)
                        oldclaw[l] = minimum(symmat[l,[i,j]])
                    end
                    process_maxsd_one2i!(i,j,dij,stepsize,
                                         facecount,numpairs,
                                         jrv,jcp,jz,
                                         r,z,pflist,
                                         oldclaw,zll,colsum,
                                         rt,ct,zt,
                                         izfull,ff2pv,
                                         pmhist,fpi,
                                         npsupp)

                    process_maxsd_i2i!(i,j,dij,stepsize,
                                       facecount,numpairs,
                                       jrv,jcp,jz,
                                       r,z,pflist,
                                       oldclaw,zll,colsum,
                                       rt,ct,zt,
                                       izfull,ff2pv,
                                       pmhist,fpi,
                                       npsupp)

                    process_maxsd_i2j!(i,j,dij,stepsize,
                                       facecount,numpairs,
                                       jrv,jcp,jz,
                                       r,z,pflist,
                                       oldclaw,zll,colsum,
                                       rt,ct,zt,
                                       izfull,ff2pv,
                                       pmhist,fpi,
                                       npsupp)

                    process_maxsd_j2j!(i,j,dij,stepsize,
                                       facecount,numpairs,
                                       jrv,jcp,jz,
                                       r,z,pflist,
                                       oldclaw,zll,colsum,
                                       rt,ct,zt,
                                       izfull,ff2pv,
                                       pmhist,fpi,
                                       npsupp)

                    process_maxsd_j2end!(i,j,dij,stepsize,
                                         facecount,numpairs,
                                         jrv,jcp,jz,
                                         r,z,pflist,
                                         oldclaw,zll,colsum,
                                         rt,ct,zt,
                                         izfull,ff2pv,
                                         pmhist,fpi,
                                         npsupp)
                end
            end
            # update the column pattern and the total number of nonzeros
            # encountered per codim2 face
            c[i+1] = facecount[1]+1
            if sd == maxsd
                colsum[jrv[cran(jcp,i)]] .+= 1
            end
        end
        delrange = c[end]:length(r)
        deleteat!(r,delrange)
        deleteat!(z,delrange)
        deleteat!(pflist,(numpairs[1]+1):nl)
        if sd == maxsd
            r = translatorvecb[r]
        end
        firstv[sd] = c
        farfaces[sd] = r
        prepairs[sd] = pflist
        grain[sd] = z
        if isempty(farfaces[sd])
            for nextcard = (sd+1):maxsd
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

function process_sd_lt_maxsd!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},pmhist::Array{Int,2},
        npsupp::BitArray{1})
    for k = cran(jcp,j)
        kk = jrv[k]
        farfilt = jz[k]
        if izfull[kk]>0
            claw = min(izfull[kk],dij)
            faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
            if claw >= farfilt && npsupp[k]
                pairupdate!(k,facecount,pflist,numpairs,npsupp,1)
            end
        end
    end
end

function process_sd_onelt_maxsd_1!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},pmhist::Array{Int,2},
        npsupp::BitArray{1})
    for k = cran(jcp,j)
        kk = jrv[k]
        farfilt = jz[k]
        if izfull[kk] > 0
            claw = min(izfull[kk],dij)
            faceupdate!(facecount,r,z,k,min(farfilt,claw),stepsize)
            if npsupp[k] && (claw >= farfilt)
                pairupdatedeluxe!(k,i,j,numpairs,facecount,pflist,ff2pv,npsupp,pmhist)
            end
        end
    end
end

function process_maxsd_one2i!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1})
    for l = 1:(i-1)
        if fpi[l,j]<fpi[l+1,j]
            ocl = oldclaw[l]
            if ocl < dij
                process_maxsd_one2i_subroutine!(i,j,dij,stepsize,
                                                facecount,numpairs,
                                                jrv,jcp,jz,
                                                r,z,pflist,
                                                oldclaw,zll,colsum,
                                               rt,ct,zt,
                                               izfull,ff2pv,
                                               pmhist,fpi,
                                               npsupp,
                                              l,ocl)
            end
        end
    end
end

function process_maxsd_one2i_subroutine!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1},
        l::Int,ocl::Int)
    for k = fpi[l,j]:(fpi[l+1,j]-1)
        ## may have to reindex this
        kk = jrv[k]
        farfilt = jz[k]
        if zll[kk] <= ocl
            break
        elseif oldclaw[l] < min(farfilt,izfull[kk])
            claw = min(izfull[kk],dij)
            if claw >= farfilt
                if npsupp[k]
                    faceupdate!(facecount,r,z,k,farfilt,stepsize)
                    pairupdate!(k,facecount,pflist,numpairs,npsupp,3)
                    ff2pv[k] = i
                elseif oldclaw[ff2pv[k]] >= farfilt
                    continue
                elseif saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
                    faceupdate!(facecount,r,z,k,farfilt,stepsize)
                end
            elseif (claw>0) && saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
                faceupdate!(facecount,r,z,k,claw,stepsize)
            end
        end
    end
end

function process_maxsd_i2i!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1})
    for k = fpi[i,j]:(fpi[i+1,j]-1)
        kk = jrv[k]
        farfilt = jz[k]
        if dij >= farfilt && npsupp[k]
            faceupdate!(facecount,r,z,k,farfilt,stepsize)
            pairupdate!(k,facecount,pflist,numpairs,npsupp,4)
            ff2pv[k] = i
        else
            farfilt = min(dij,farfilt)
            if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
                faceupdate!(facecount,r,z,k,farfilt,stepsize)
            end
        end
    end
end

function process_maxsd_i2j!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1})
    for k = fpi[i+1,j]:(fpi[j,j]-1)
        kk = jrv[k]
        if izfull[kk] > 0
            farfilt = jz[k]
            claw = min(izfull[kk],dij)
            if claw >= farfilt && npsupp[k]
                faceupdate!(facecount,r,z,k,farfilt,stepsize)
                pairupdate!(k,facecount,pflist,numpairs,npsupp,5)
                ff2pv[k] = i
            else
                farfilt = min(claw,farfilt)
                if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
                    faceupdate!(facecount,r,z,k,farfilt,stepsize)
                end
            end
        end
    end
end

function process_maxsd_j2j!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1})
    for k = fpi[j,j]:fpi[j+1,j]-1
        kk = jrv[k]
        if izfull[kk]>0
            claw = min(izfull[kk],dij)
            if saveface(ct,kk,colsum,claw,oldclaw,rt,zt)
                faceupdate!(facecount,r,z,k,claw,stepsize)
            end
        end
    end
end

function process_maxsd_j2end!(
        i::Int,j::Int,dij::Int,stepsize::Int,
        facecount::Array{Int,1},numpairs::Array{Int,1},
        jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
        r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
        oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
        rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
        izfull::Array{Int,1},ff2pv::Array{Int,1},
        pmhist::Array{Int,2},fpi::Array{Int,2},
        npsupp::BitArray{1})
    for k = fpi[j+1,j]:jcp[j+1]-1
        kk = jrv[k]
        if izfull[kk]>0
            farfilt = jz[k]
            claw = min(izfull[kk],dij)
            if claw >= farfilt && npsupp[k]
                faceupdate!(facecount,r,z,k,farfilt,stepsize)
                pairupdate!(k,facecount,pflist,numpairs,npsupp,6)
                ff2pv[k] = i
            else
                farfilt = min(claw,farfilt)
                if saveface(ct,kk,colsum,farfilt,oldclaw,rt,zt)
                    faceupdate!(facecount,r,z,k,farfilt,stepsize)
                end
            end
        end
    end
end

function pairupdate!(k::Int,facecount::Array{Int,1},pflist::Array{Int,1},numpairs::Array{Int,1},npsupp::BitArray{1},iterateNumber)
    numpairs[1] += 1
    pflist[numpairs[1]] = facecount[1]
    npsupp[k] = false
end

function pairupdatedeluxe!(k::Int,i::Int,j::Int,numpairs::Array{Int,1},facecount::Array{Int,1},pflist::Array{Int,1},ff2pv::Array{Int,1},npsupp::BitArray{1},pmhist::Array{Int,2})
    numpairs[1] += 1
    pmhist[i,j] += 1
    npsupp[k] = false
    pflist[numpairs[1]] = facecount[1]
    ff2pv[k] = i
end

function faceupdate!(facecount::Array{Int,1},r::Array{Int,1},z::Array{Int,1},k::Int,farfilt::Int,stepsize::Int)
    facecount[1] += 1
    if facecount[1] > length(r)
        append!(r,Array{Int}(undef,stepsize))
        append!(z,Array{Int}(undef,stepsize))
    end
    r[facecount] .= k
    z[facecount] .= farfilt
end

function faceupdatedeluxe!(facecount::Array{Int,1},r::Array{Int,1},z::Array{Int,1},k::Int,farfilt::Int,stepsize::Int,s::Array{Int,1},i::Int)
    facecount[1] += 1
    if facecount[1] > length(r)
        append!(r,Array{Int}(undef,stepsize))
        append!(z,Array{Int}(undef,stepsize))
        append!(s,Array{Int}(undef,stepsize))
    end
    r[facecount] .= k
    z[facecount] .= farfilt
    s[facecount] .= i
end

function saveface(ct::Array{Int,1},kk::Int,colsum::Array{Int,1},farfilt::Int,oldclaw::Array{Int,1},rt::Array{Int,1},zt::Array{Int,1})
    for l = ct[kk]:colsum[kk]
        if zt[l] >= farfilt && oldclaw[rt[l]] >= farfilt
            return false
        end
    end
    true
end

function processfpi!(pmhist::Array{Int,2},fpi::Array{Int,2},jcp::Array{Int,1},jrv::Array{Int,1},ff2pv::Array{Int,1},m::Integer)
    for p = 1:m
        for q = jcp[p]:jcp[p+1]-1
            pmhist[ff2pv[jrv[q]],p]+=1
        end
    end
    for p = 1:m
        fpi[1,p] = jcp[p]
        for q = 1:m
            fpi[q+1,p] = fpi[q,p]+pmhist[q,p]
        end
    end
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
    grain::Array{Int,1} = grain_cell[2]
    farfaces::Array{Int,1} = farfaces_cell[2]
    firstv::Array{Int,1} = firstv_cell[2]

    numverts = length(firstv)-1
    numedges = length(farfaces)
    stepsize = size(symmat,1)^2
    facecount = [0]
    numpairs = 0

    closefaces = Array{Int}(undef,numedges)
    for i = 1:m
        closefaces[cran(firstv,i)] .= i
    end
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
                                    faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
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
                faceupdatedeluxe!(facecount,r,z,kp,dijk,stepsize,s,i)
            end
        end
        ####
    end

    num3faces = facecount[1]
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

