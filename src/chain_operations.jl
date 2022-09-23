#!/usr/bin/env julia

function persistf2(
        farfaces::Array{Array{Int,1},1},
        firstv::Array{Array{Int,1},1},
        prepairs::Array{Array{Int,1},1},
        grain::Array{Array{Int,1},1})
    
    maxsd = length(farfaces) - 1
    m = length(firstv[1]) - 1

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
        lowbasisnames = sd > 2 ? phi[sd-1] : Array{Int}(undef,0)
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
            morselu!(Mrv, higfilttemp, Mcp, lowfilttemp, lowlabtemp, higlabtemp, pplow, pphig, Mm)
        trv[sd] = lowlab[Srv]
        tcp[sd] = Scp
        tid[sd] = lowlab[tlab]
        plo[sd] = lowlab[Splows]
        phi[sd] = higlab[Sphigs]
        maxnzs[sd] = maxnz
    end
    
    (trv=trv, tcp=tcp, plo=plo, phi=phi, tid=tid, maxnzs=maxnzs)
end

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
        Mm::Int)
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

    schurit4!(Mrv,Mcp,Mm,Mn,Mn0,
              rowlab,collab,
              Jprows,Jpcols,numjunpairs,
              Sprows,Spcols,numsenpairs,
              comprows,compcols,
              Trv,Tcp,Srv,Scp;
              updatetransform=true)
    maxnz = max(maxnz,Mcp[Mn[1]+1])

    #gc()
    rowfilt = Array{Int}(undef,length(comprows)); colfilt = Array{Int}(undef,length(compcols))
    counter = 0
    while Mcp[Mn[1]+1] > 1
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

        schurit4!(Mrv,Mcp,Mm,Mn,Mn0, rowlab,collab, Jprows,Jpcols,numjunpairs,
                  Sprows,Spcols,numsenpairs, comprows,compcols,
                  Trv,Tcp,Srv,Scp)
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
    Tcp .+= lastSrowmarker
    append!(Scp,Tcp)
    append!(Srv,Trv[1:lastTrowmarker])
    tlab = Spcols[1:numsenpairs[1]]
    append!(tlab,collab[1:Mn[1]])
    return Srv,Scp,Sprows,Spcols,tlab,maxnz
end

function getPairsLightWrite2!(
	rowval::Array{Int,1},
	colptr::Array{Int,1},
	rowfilt::Array{Int,1},
	colfilt::Array{Int,1},
	m::Integer,
	n::Integer,
	prows::Array{Int,1},
	pcols::Array{Int,1},
	numpairs::Array{Int,1})

	col2firstplace = zeros(Int,n)
	rowwisesum = zeros(Int,m)

	for j = 1:n
		firstplace = colptr[j]
		if firstplace < colptr[j+1]
			firstrow = rowval[firstplace]
			rowwisesum[firstrow]+=1
			if firstplace < colptr[j+1]-1
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newfilt = rowfilt[newrow]
					rowwisesum[newrow]+=1
					if newfilt > filt
						filt = newfilt
						firstplace = newplace
					end
				end
			end
			col2firstplace[j]=firstplace
		end
	end
	colwisesum = colsupportsum(colptr,n) # note allow extra on end
	colfiltptr = getcolptr2(colfilt,n)	# note allow extra on end
	colwisesumlinearized = integersinsameorderbycolumn2(colwisesum,colfiltptr)
	colnamesinorder = Array{Int}(undef,n)
	colnamesinorder[colwisesumlinearized]=1:n
	ncoveredsupp = trues(m)
	pairmarker = 0
	for jp = 1:n
		j = colnamesinorder[jp]
		if col2firstplace[j]>0
			firstplace = col2firstplace[j]
			firstrow = rowval[firstplace]
			if firstplace==colptr[j+1]-1
				if ncoveredsupp[firstrow]	#######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			else
				firstweight = rowwisesum[firstrow]
				filt = rowfilt[firstrow]
				for newplace = (firstplace+1):(colptr[j+1]-1)
					newrow = rowval[newplace]
					newweight = rowwisesum[newrow]
					if rowfilt[newrow]==filt && newweight <= firstweight && !ncoveredsupp[firstrow] && ncoveredsupp[newrow]
						firstweight = newweight
						firstrow = newrow
					end
				end
				if ncoveredsupp[firstrow] ######&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			end
			for ip = cran(colptr,j)
				ncoveredsupp[rowval[ip]]=false
			end
		end
	end

	numpairs[1] = pairmarker
	return
end

function colsupportsum(colptr,n::Int)
	x = Array{Int}(undef,n)
	@inbounds begin
	for i = 1:n
		x[i] = colptr[i+1]-colptr[i]
	end
	end
	return x
end

function getcolptr2(orderedpositiveintegerlist::Array{Int,1}, howfartolookbeforestopping::Int)
	#### please note: order can be ascending or descending
	v = orderedpositiveintegerlist
	if isempty(v)
		return []
	end
	colptr = Array{Int}(undef,length(v)+1)
	colptr[1] = 1
	transitioncounter = 1
	currentvalue = v[1]
	for i = 2:howfartolookbeforestopping
		if v[i]!=currentvalue
			transitioncounter+=1
			colptr[transitioncounter]=i
			currentvalue = v[i]
		end
	end
	colptr[transitioncounter+1]=howfartolookbeforestopping+1
	deleteat!(colptr,(transitioncounter+2):(length(v)+1))
	return colptr
end

function boundarylimit_Lionly(trv, tcp, tid, sd::Int)
    if isempty(tid[sd])
        Lirv = zeros(Int,0)
        Licp = ones(Int,1)
    else
        tsize = length(tcp[sd])-1
        Lrv = copy(trv[sd])
        lowtranslator = zeros(Int,maximum(tid[sd]))
        lowtranslator[tid[sd]] = 1:length(tid[sd])
        Lrv .= lowtranslator[Lrv]
        Lirv,Licp = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv,tcp[sd])
        Lirv,Licp = transposeLighter(Lirv,Licp,tsize)
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

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd::Int)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])
    brv,bcp = maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
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
    # note there's no +1 b/c bcp is already 1 elt. too long
    append!(bcp,fill(bcp[end],numnlpl-nump))
    return brv,bcp
end

function unpack!(grain, farfaces, firstv, trv, tcp, plo, phi, tid, maxsd::Int)
    Lirv = Array{Array{Int,1},1}(undef,maxsd)
    Licp = Array{Array{Int,1},1}(undef,maxsd)
    Lrv  = Array{Array{Int,1},1}(undef,maxsd)
    Lcp  = Array{Array{Int,1},1}(undef,maxsd)
    Rrv  = Array{Array{Int,1},1}(undef,maxsd)
    Rcp  = Array{Array{Int,1},1}(undef,maxsd)

    N = maxsd-1
    if maxsd >= 2
        Lirv[maxsd], Licp[maxsd] = boundarylimit_Lionly(trv, tcp, tid, maxsd)
    elseif maxsd == 1
        Lirv[maxsd] = Array{Int,1}(undef,0)
        Licp[maxsd] = ones(Int,1)
    end
    Lrv[maxsd] = Array{Int,1}(undef,0)
    Lcp[maxsd] = Array{Int,1}(undef,0)
    Rrv[maxsd] = Array{Int,1}(undef,0)
    Rcp[maxsd] = Array{Int,1}(undef,0)

    for i = 2:N
        Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,i)
        @assert !isempty(Lcp[i]) "unpack!: Lcp[i] = 0 and i = $(i)"
    end

    Lirv[1]	= Array{Int}(undef,0)
    Lrv[1] = Array{Int}(undef,0)
    Rrv[1] = Array{Int}(undef,0)
    Licp[1] = ones(Int,1)
    Lcp[1] = ones(Int,1)
    Rcp[1] = ones(Int,1)

    # for each dimension one stores an array of arrays
    cyclerep = fill(Array{Array{Int,1},1}(undef,0), maxsd)
    for sd = 2:maxsd
        dim = sd-2
        m = nnzbars(grain, plo, phi, tid; dim=dim)
        cyclenames = barname2cyclename(grain, plo, phi, tid, 1:m; dim=dim)
        cyclerep[sd] = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenames)
    end
    
    (cyclerep=cyclerep, Lrv=Lrv, Lcp=Lcp, Lirv=Lirv, Licp=Licp, Rrv=Rrv, Rcp=Rcp)
end

