#!/usr/bin/env julia

function persistf2(farfaces::Vector{Vector{Int}}, firstv::Vector{Vector{Int}},
        prepairs::Vector{Vector{Int}}, grain::Vector{Vector{Int}})
    
    maxsd = length(farfaces) - 1
    m = length(firstv[1]) - 1

    trv = Vector{Vector{Int}}(undef,maxsd+1)
    tcp = Vector{Vector{Int}}(undef,maxsd+1)
    phi = Vector{Vector{Int}}(undef,maxsd+1)
    plo = Vector{Vector{Int}}(undef,maxsd+1)
    tid = Vector{Vector{Int}}(undef,maxsd+1)
    for i in [1, maxsd+1]
        tcp[i] = Int[1]
        trv[i] = Int[]
        tid[i] = Int[]
        phi[i] = Int[]
        plo[i] = Int[]
    end

    maxnzs = zeros(Int, maxsd+1)

    for sd = 2:maxsd
        lowbasisnames = sd > 2 ? phi[sd-1] : Int[]
        Mrv, Mcp, lowlab, higlab, Mm =
            filteredmatrixfromfarfaces(farfaces,firstv,prepairs,grain,sd,lowbasisnames)
        higfilttemp = grain[sd][higlab]
        lowfilttemp = grain[sd-1][lowlab]
        # NB: It is critical that the columns of the input array
        # be ordered according to filtration; in particular, the entries of
        # lowfiltemp should increase monotonically
        Srv,Scp,Sphigs,Splows,tlab,maxnz =
        morselu!(Mrv, higfilttemp, Mcp, lowfilttemp, 
                 1:length(lowlab), 
                 1:length(higlab),
                 length(prepairs[sd]):-1:1, 
                 length(prepairs[sd]):-1:1)
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
function morselu!(Mrv::Vector{Int}, Mrowgrain::Vector{Int}, Mcp::Vector{Int},
        Mcolgrain::Vector{Int}, lowlab, higlab, pplow, pphig)
    rowlab = collect(higlab)
    collab = collect(lowlab)
    
    Mm = length(higlab)
    Mn0 = Mn = length(lowlab)
    maxnz = Mcp[Mn+1]
    
    maxnumpairs = min(Mm, Mn)
    numjunpairs = length(pplow)
    numsenpairs = 0
    Sprows = Vector{Int}(undef,maxnumpairs)
    Spcols = Vector{Int}(undef,maxnumpairs)
    Jprows = Vector{Int}(undef,maxnumpairs)
    Jpcols = Vector{Int}(undef,maxnumpairs)
    Jprows[1:numjunpairs] = pphig
    Jpcols[1:numjunpairs] = pplow
    comprows = numjunpairs+1 : Mm
    compcols = numjunpairs+1 : Mn
    
    Trv = Int[]
    Srv = Int[]
    Tcp = ones(Int, Mn+1)
    Scp = ones(Int, Mn+1)

    Mm, Mn, numsenpairs = schurit4!(Mrv,Mcp,Mm,Mn,Mn0, rowlab,collab, Jprows,Jpcols,numjunpairs,
              Sprows,Spcols,numsenpairs, comprows, compcols, Trv,Tcp,Srv,Scp)
    maxnz = max(maxnz,Mcp[Mn+1])

    rowfilt = Vector{Int}(undef, length(comprows))
    colfilt = Vector{Int}(undef, length(compcols))
    while Mcp[Mn+1] > 1
        rowfilt[1:Mm] = Mrowgrain[rowlab[1:Mm]]
        colfilt[1:Mn] = Mcolgrain[collab[1:Mn]]
        numjunpairs = getPairsLightWrite2!(Mrv,Mcp,rowfilt,colfilt,Mm,Mn,Jprows,Jpcols)
        comprows = intervalcomplementuniqueunsortedinput(Jprows[1:numjunpairs], Mm)
        compcols = intervalcomplementuniqueunsortedinput(Jpcols[1:numjunpairs], Mn)

        Mm, Mn, numsenpairs = schurit4!(Mrv, Mcp, Mm, Mn, Mn0, rowlab, collab, Jprows,Jpcols,numjunpairs,
                  Sprows,Spcols,numsenpairs, comprows, compcols, Trv,Tcp,Srv,Scp)
        maxnz = max(maxnz,Mcp[Mn+1])
    end
    lastSrowmarker = Scp[numsenpairs+1]-1
    lastTrowmarker = Tcp[Mn+1]-1
    deleteat!(Srv,(lastSrowmarker+1):length(Srv))
    deleteat!(Trv,(lastTrowmarker+1):length(Trv))
    deleteat!(Scp,(numsenpairs+1):length(Scp))
    deleteat!(Tcp,(Mn+2):length(Tcp))
    deleteat!(Sprows,(numsenpairs+1):maxnumpairs)
    deleteat!(Spcols,(numsenpairs+1):maxnumpairs)
    Tcp .+= lastSrowmarker
    append!(Scp, Tcp)
    append!(Srv, Trv[1:lastTrowmarker])
    tlab = Spcols[1:numsenpairs]
    append!(tlab, collab[1:Mn])
    Srv,Scp,Sprows,Spcols,tlab,maxnz
end

function getPairsLightWrite2!(rowval::Vector{Int}, colptr::Vector{Int},
        rowfilt::Vector{Int}, colfilt::Vector{Int}, m::Int, n::Int,
        prows::Vector{Int}, pcols::Vector{Int})

	col2firstplace = zeros(Int, n)
	rowwisesum = zeros(Int, m)

	for j = 1:n
		firstplace = colptr[j]
		if firstplace < colptr[j+1]
			firstrow = rowval[firstplace]
			rowwisesum[firstrow] += 1
			if firstplace < colptr[j+1] - 1
				filt = rowfilt[firstrow]
				for newplace = firstplace+1:colptr[j+1]-1
					newrow = rowval[newplace]
					newfilt = rowfilt[newrow]
					rowwisesum[newrow] += 1
					if newfilt > filt
						filt = newfilt
						firstplace = newplace
					end
				end
			end
			col2firstplace[j] = firstplace
		end
	end
	# note allow extra on end
	colwisesum = colptr[2:n+1] .- colptr[1:n]
	# note allow extra on end
	colfiltptr = getcolptr2(colfilt,n)
	colwisesumlinearized = integersinsameorderbycolumn2(colwisesum,colfiltptr)
	colnamesinorder = Array{Int}(undef,n)
	colnamesinorder[colwisesumlinearized] = 1:n
	ncoveredsupp = trues(m)
	pairmarker = 0
	for jp = 1:n
		j = colnamesinorder[jp]
		if col2firstplace[j] > 0
			firstplace = col2firstplace[j]
			firstrow = rowval[firstplace]
			if firstplace == colptr[j+1] - 1
				if ncoveredsupp[firstrow]	#&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
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
				if ncoveredsupp[firstrow] #&&  (pairmarker==0 || rowwisesum[firstrow]<10000)
					pairmarker+=1
					prows[pairmarker]=firstrow
					pcols[pairmarker]=j
				end
			end
            ncoveredsupp[rowval[cran(colptr,j)]] .= false
		end
	end

	pairmarker
end

"""
- v: ordered positive integer list. order can be ascending or descending.
- n: how far to look before stopping
"""
function getcolptr2(v::Vector{Int}, n::Int)
	isempty(v) && return []
	colptr = Array{Int}(undef, length(v)+1)
	colptr[1] = 1
	transitioncounter = 1
	currentvalue = v[1]
	for i = 2:n
		if v[i] != currentvalue
			transitioncounter += 1
			colptr[transitioncounter] = i
			currentvalue = v[i]
		end
	end
	colptr[transitioncounter+1] = n + 1
	deleteat!(colptr, transitioncounter+2:length(v)+1)
	colptr
end

function boundarylimit_Lionly(trv, tcp, tid, sd::Int)
    isempty(tid[sd]) && return zeros(Int,0), ones(Int,1)
    tsize = length(tcp[sd]) - 1
    Lrv = trv[sd]
    lowtranslator = zeros(Int, maximum(tid[sd]))
    lowtranslator[tid[sd]] = 1:length(tid[sd])
    Lrv .= lowtranslator[Lrv]
    Lirv, Licp = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Lrv, tcp[sd])
    Lirv, Licp = transposeLighter(Lirv, Licp, tsize)
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
    Lirv,Licp = morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(trv,tcp)
    Lirv,Licp = transposeLighter(Lirv,Licp,numnlpl)
    Lrv,Lcp = transposeLighter(trv,tcp,numnlpl)
    brv,bcp = spmmF2silentLeft(Lrv,Lcp,brv,bcp,numnlpl)
    Rrv,Rcp = morseInverseF2orderedColsUnsortedRowsInSilentOut(brv,bcp)
    return Lrv,Lcp,Lirv,Licp,Rrv,Rcp
end

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd::Int)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])
    brv, bcp = maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl - length(plo[sd-1])

    lowtranslator = zeros(Int,numl)

    lowtranslator[tid[sd]] = 1:numnlpl
    brv = ff2aflight(farfaces,firstv,sd,phi[sd])[:]
    bcp = collect(1:sd:(nump*sd+1))
    supportedmatrix!(brv, bcp, tid[sd], 1:nump, numl)
    # note there's no +1 b/c bcp is already 1 elt. too long
    return lowtranslator[brv], [bcp; fill(bcp[end], numnlpl-nump)]
end

function unpack!(grain, farfaces, firstv, trv, tcp, plo, phi, tid, maxsd::Int)
    Lirv = Vector{Vector{Int}}(undef,maxsd)
    Licp = Vector{Vector{Int}}(undef,maxsd)
    Lrv  = Vector{Vector{Int}}(undef,maxsd)
    Lcp  = Vector{Vector{Int}}(undef,maxsd)
    Rrv  = Vector{Vector{Int}}(undef,maxsd)
    Rcp  = Vector{Vector{Int}}(undef,maxsd)

    N = maxsd-1
    if maxsd >= 2
        Lirv[maxsd], Licp[maxsd] = boundarylimit_Lionly(trv, tcp, tid, maxsd)
    elseif maxsd == 1
        Lirv[maxsd] = Int[]
        Licp[maxsd] = ones(Int, 1)
    end
    Lrv[maxsd] = Int[]
    Lcp[maxsd] = Int[]
    Rrv[maxsd] = Int[]
    Rcp[maxsd] = Int[]

    for i = 2:N
        Lrv[i],Lcp[i],Lirv[i],Licp[i],Rrv[i],Rcp[i] = boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,i)
        @assert !isempty(Lcp[i]) "unpack!: Lcp[i] = 0 and i = $(i)"
    end

    Lirv[1]	= Int[]
    Lrv[1] = Int[]
    Rrv[1] = Int[]
    Licp[1] = ones(Int,1)
    Lcp[1] = ones(Int,1)
    Rcp[1] = ones(Int,1)

    # for each dimension one stores an array of arrays
    cyclerep = fill([Int[]], maxsd)
    for sd = 2:maxsd
        dim = sd-2
        m = nnzbars(grain, plo, phi, tid; dim=dim)
        cyclenames = barname2cyclename(grain, plo, phi, tid, 1:m; dim=dim)
        cyclerep[sd] = getcycle(farfaces,firstv,Lirv,Licp,Lrv,Lcp,Rrv,Rcp,plo,phi,tid,sd,cyclenames)
    end
    
    (cyclerep=cyclerep, Lrv=Lrv, Lcp=Lcp, Lirv=Lirv, Licp=Licp, Rrv=Rrv, Rcp=Rcp)
end

