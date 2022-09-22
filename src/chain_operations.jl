#!/usr/bin/env julia

"""
- the output array tid has size equal to the number of columns of the
input array, NOT the column-rank of the input array
-	the columns of M must be ordered by grain (ascending)
-	the first (rank of M) elements of tlab index the complete set
of nonzero columns in the reduced matrix
"""
function morselu!(
        Mrv::Array{Int,1},
        Mrowgrain::Array{Int,1},
        Mcp::Array{Int,1},
        Mcolgrain::Array{Int,1},
        lowlab::Array{Int,1},
        higlab::Array{Int,1},
        pplow::Array{Int,1},
        pphig::Array{Int,1},
        Mm::Int;
        storetransform::Bool=true,
        diagnostic::Bool=false)
    rowlab = higlab
    collab = lowlab
    
    Mm = [length(higlab)]
    Mn = [length(lowlab)]
    Mn0 = Mn[1]
    maxnz = Mcp[Mn[1]+1]
    
    maxnumpairs = min(Mm[1],Mn[1]); numjunpairs = [length(pplow)]; numsenpairs = [0]
    Sprows=Array{Int}(undef,maxnumpairs);Spcols=Array{Int}(undef,maxnumpairs);
    Jprows=Array{Int}(undef,maxnumpairs);Jpcols=Array{Int}(undef,maxnumpairs);
    Jprows[1:numjunpairs[1]]=pphig;Jpcols[1:numjunpairs[1]]=pplow
    comprows = convert(Array{Int,1},(numjunpairs[1]+1):Mm[1])
    compcols = convert(Array{Int,1},(numjunpairs[1]+1):Mn[1])
    
    Trv=Array{Int}(undef,0);Srv = Array{Int}(undef,0)
    Tcp=ones(Int,Mn[1]+1);Scp=ones(Int,Mn[1]+1)

    if diagnostic
        numsenpairsOLD = numsenpairs[1]
    end
    schurit4!(Mrv,Mcp,Mm,Mn,Mn0,
              rowlab,collab,
              Jprows,Jpcols,numjunpairs,
              Sprows,Spcols,numsenpairs,
              comprows,compcols,
              Trv,Tcp,Srv,Scp;
              updatetransform = storetransform)
    maxnz = max(maxnz,Mcp[Mn[1]+1])

    #gc()
    rowfilt = Array{Int}(undef,length(comprows)); colfilt = Array{Int}(undef,length(compcols))
    counter = 0
    while Mcp[Mn[1]+1]>1
        counter+=1
        for i = 1:Mm[1]
            rowfilt[i] = Mrowgrain[rowlab[i]]
        end
        for j = 1:Mn[1]
            colfilt[j] = Mcolgrain[collab[j]]
        end
        getPairsLightWrite2!(Mrv,Mcp,rowfilt,colfilt,Mm[1],Mn[1],Jprows,Jpcols,numjunpairs)
        comprows = intervalcomplementuniqueunsortedinput(Jprows[1:numjunpairs[1]],Mm[1])
        compcols = intervalcomplementuniqueunsortedinput(Jpcols[1:numjunpairs[1]],Mn[1])
        if diagnostic
            numsenpairsOLD = numsenpairs[1]
        end

        schurit4!(Mrv,Mcp,Mm,Mn,Mn0, rowlab,collab, Jprows,Jpcols,numjunpairs,
                  Sprows,Spcols,numsenpairs, comprows,compcols,
                  Trv,Tcp,Srv,Scp; updatetransform=storetransform)
        maxnz = max(maxnz,Mcp[Mn[1]+1])
    end
    lastSrowmarker = Scp[numsenpairs[1]+1]-1
    lastTrowmarker = Tcp[Mn[1]+1]-1
    deleteat!(Srv,(lastSrowmarker+1):length(Srv))
    deleteat!(Trv,(lastTrowmarker+1):length(Trv))
    deleteat!(Scp,(numsenpairs[1]+1):length(Scp))
    deleteat!(Tcp,(Mn[1]+2):length(Tcp))
    deleteat!(Sprows,(numsenpairs[1]+1):maxnumpairs)
    deleteat!(Spcols,(numsenpairs[1]+1):maxnumpairs)
    Tcp.+=lastSrowmarker
    append!(Scp,Tcp)
    append!(Srv,Trv[1:lastTrowmarker])
    tlab = Spcols[1:numsenpairs[1]]
    append!(tlab,collab[1:Mn[1]])
    return Srv,Scp,Sprows,Spcols,tlab,maxnz
end

function persistf2_core_vr(
        farfaces::Array{Array{Int,1},1},
        firstv::Array{Array{Int,1},1},
        prepairs::Array{Array{Int,1},1},
        grain::Array{Array{Int,1},1},
        maxsd::Integer)

    storetransform = true
    
    m = length(firstv[1])-1

    trv::Array{Array{Int,1},1} = Array{Array{Int,1},1}(undef,maxsd+1)
    tcp::Array{Array{Int,1},1} = Array{Array{Int,1},1}(undef,maxsd+1)
    phi::Array{Array{Int,1},1} = Array{Array{Int,1},1}(undef,maxsd+1)
    plo::Array{Array{Int,1},1} = Array{Array{Int,1},1}(undef,maxsd+1)
    tid::Array{Array{Int,1},1} = Array{Array{Int,1},1}(undef,maxsd+1)
    for i in [1, maxsd+1]
        tcp[i] = [1];
        trv[i] = Array{Int}(undef,0)
        tid[i] = Array{Int}(undef,0)
        phi[i] = Array{Int}(undef,0)
        plo[i] = Array{Int}(undef,0)
    end

    maxnzs = zeros(Int, maxsd+1)

    for sd = 2:maxsd
        if sd > length(farfaces)
            continue
        elseif sd > 2
            lowbasisnames = phi[sd-1]
        else
            lowbasisnames = Array{Int}(undef,0)
        end
        Mrv::Array{Int,1},
        Mcp::Array{Int,1},
        lowlab::Array{Int,1},
        higlab::Array{Int,1},
        Mm =
        filteredmatrixfromfarfaces(farfaces,firstv,prepairs,grain,sd,lowbasisnames)
        lowlabtemp = convert(Array{Int,1},1:length(lowlab))
        higlabtemp = convert(Array{Int,1},1:length(higlab))
        higfilttemp = grain[sd][higlab]
        lowfilttemp = grain[sd-1][lowlab]
        pplow = convert(Array,length(prepairs[sd]):-1:1)
        pphig = convert(Array,length(prepairs[sd]):-1:1)

        # NB: It is critical that the columns of the input array
        # be ordered according to filtration; in particular, the entries of
        # lowfiltemp should increase monotonically
        Srv,Scp,Sphigs,Splows,tlab,maxnz =
        morselu!(Mrv, higfilttemp, Mcp, lowfilttemp, lowlabtemp, higlabtemp, pplow, pphig, Mm,
                 storetransform = storetransform)
        trv[sd] = lowlab[Srv]
        tcp[sd] = Scp
        tid[sd] = lowlab[tlab]
        plo[sd] = lowlab[Splows]
        phi[sd] = higlab[Sphigs]
        maxnzs[sd]= maxnz
    end
    return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2!(D::Dict; maxsd=0, dictionaryoutput::Bool = true)
    farfaces = D["farfaces"]
    firstv = D["firstv"]
    prepairs = D["prepairs"]
    grain = D["grain"]
    if maxsd == 0 maxsd = length(farfaces) - 1 end

    trv,tcp,plo,phi,tid,maxnzs =
    persistf2_core_vr(farfaces, firstv, prepairs, grain, maxsd::Integer)
    if dictionaryoutput == true
        D["trv"] = trv
        D["tcp"] = tcp
        D["tid"] = tid
        D["plo"] = plo
        D["phi"] = phi
        D["maxnz"] = maxnzs
        return D
    else
        return trv, tcp, plo, phi, tid
    end
end

function boundarylimit_Lionly(trv,tcp,tid,sd)
    if isempty(tid[sd])
        Lirv = zeros(Int,0)
        Licp = ones(Int,1)
    else
        tsize			= length(tcp[sd])-1
        Lrv = copy(trv[sd])
        lowtranslator = zeros(Int,maximum(tid[sd]))
        lowtranslator[tid[sd]] = 1:length(tid[sd])
        Lrv .= lowtranslator[Lrv]
        Lirv,Licp   = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv,tcp[sd])
        Lirv,Licp   = transposeLighter(Lirv,Licp,tsize)
    end
    return Lirv,Licp
end

"""
- L, Li, and R are all indexed by tid, just like (trv,tcp)
- L (not Li) is the transpose of (trv,tcp)
- up to permutation and deletion of some zero rows, RLM is the fundamental
cycle matrix of the d-dimensional boundary operator with respect
to the doubly-minimal basis we have constructed, where M the submatrix of
the boundary with rows indexed by tid.
"""
function boundarylimit_core(brv,bcp,trv,tcp,tid,numl,nump,numnlpl)
    lowtranslator = zeros(Int,numl)
    lowtranslator[tid] = 1:numnlpl
    trv = lowtranslator[trv]
    Lirv,Licp	= morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(trv,tcp)
    Lirv,Licp   = transposeLighter(Lirv,Licp,numnlpl)
    Lrv,Lcp = transposeLighter(trv,tcp,numnlpl)
    brv,bcp = spmmF2silentLeft(Lrv,Lcp,brv,bcp,numnlpl)
    Rrv,Rcp = morseInverseF2orderedColsUnsortedRowsInSilentOut(brv,bcp)
    return Lrv,Lcp,Lirv,Licp,Rrv,Rcp
end

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])
    brv,bcp = maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd)
    numl = length(cp[sd-1])-1
    nump = length(phi[sd])
    numnlpl = numl-length(phi[sd-1])
    brv,bcp		= maxnonsingblock_cell(rv,cp,plo,phi,tid,sd)
    Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit(D::Dict,sd)
    # special case sd = 1 may possibly be degenerate
    trv = D["trv"];tcp=D["tcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
    if haskey(D,"farfaces")
        farfaces = D["farfaces"];firstv = D["firstv"]
        return boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd)
    else
        rv = D["rv"];cp = D["cp"]
        # Lrv,Lcp,Lirv,Licp,Rrv,Rcp = 
        # boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd)
        return boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd)
    end
end

function maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])

    lowtranslator = zeros(Int,numl)

    lowtranslator[tid[sd]] = 1:numnlpl
    brv = ff2aflight(farfaces,firstv,sd,phi[sd])
    brv = reshape(brv,length(brv))
    bcp = convert(Array{Int,1},1:sd:(nump*sd+1))
    supportedmatrix!(brv,bcp,tid[sd],1:nump,numl)
    brv .= lowtranslator[brv]
    append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. too long
    return brv,bcp
end

function maxnonsingblock_cell(rv,cp,plo,phi,tid,sd)
    numl = length(cp[sd-1]) - 1
    nump = length(phi[sd])
    numnlpl = numl - length(plo[sd-1])

    if isempty(phi[sd])
        brv = zeros(Int,0)
        bcp = ones(Int,numnlpl+1)
        return brv, bcp
    end

    lowtranslator = zeros(Int,numl)
    lowtranslator[tid[sd]] = 1:numnlpl
    dummy0 = zeros(Int,0)

    brv,dummy1,bcp,dummy2 = stackedsubmatrices(rv[sd],cp[sd],tid[sd],dummy0,phi[sd],numl)
    brv .= lowtranslator[brv]
    # note there's no +1 b/c bcp is already 1 elt. longer than the # of columns
    append!(bcp, fill(bcp[end], numnlpl-nump))
    return brv, bcp
end

function unpack!(D::Dict)
    maxsd = D["input"]["maxdim"]+2

    Lirv = Array{Array{Int,1},1}(undef,maxsd)
    Licp = Array{Array{Int,1},1}(undef,maxsd)
    Lrv  = Array{Array{Int,1},1}(undef,maxsd)
    Lcp  = Array{Array{Int,1},1}(undef,maxsd)
    Rrv  = Array{Array{Int,1},1}(undef,maxsd)
    Rcp  = Array{Array{Int,1},1}(undef,maxsd)

    N = maxsd-1
    if maxsd >= 2
        Lirv[maxsd],Licp[maxsd] = boundarylimit_Lionly(D["trv"],D["tcp"],D["tid"],maxsd)
    elseif maxsd == 1
        Lirv[maxsd] = Array{Int,1}(undef,0)
        Licp[maxsd] = ones(Int,1)
    end
    Lrv[maxsd] = Array{Int,1}(undef,0)
    Lcp[maxsd] = Array{Int,1}(undef,0)
    Rrv[maxsd] = Array{Int,1}(undef,0)
    Rcp[maxsd] = Array{Int,1}(undef,0)

    for i = 2:N
        Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = boundarylimit(D,i)
        if isempty(Lcp[i])
            println("ERROR MESSAGE IN unpack!: Lcp[i] = 0 and i = $(i)")
        end
    end

    Lirv[1]	= Array{Int}(undef,0)
    Lrv[1] = Array{Int}(undef,0)
    Rrv[1] = Array{Int}(undef,0)
    Licp[1] = ones(Int,1)
    Lcp[1] = ones(Int,1)
    Rcp[1] = ones(Int,1)

    D["Lrv"] = Lrv
    D["Lcp"] = Lcp
    D["Lirv"] = Lirv
    D["Licp"] = Licp
    D["Rrv"] = Rrv
    D["Rcp"] = Rcp

    # for each dimension one stores an array of arrays
    D["cyclerep"] = fill(Array{Array{Int,1},1}(undef,0),maxsd)

    for i = 2:maxsd
        dim = i-2
        m = nnzbars(D,dim=dim)
        cyclenames = barname2cyclename(D, 1:m, dim=dim)
        D["cyclerep"][i] = getcycle(D, cyclenames, dim=dim)
    end
    return
end

