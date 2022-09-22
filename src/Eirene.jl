#!/usr/bin/env julia
# WELCOME TO EIRENE!
#
# You should have received a copy of the GNU General Public License along with
# Eirene.  If not, please see <http://www.gnu.org/licenses/>.
#
# Eirene Library for Homological Algebra
# Copyright (C) 2016, 2017, 2018, 2019, 2020, 2021  Gregory Henselman
# www.gregoryhenselman.org
#
# Eirene is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Eirene is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Eirene.  If not, see <http://www.gnu.org/licenses/>.
#
# PLEASE HELP US DOCUMENT Eirene's recent work! Bibtex entries and
# contact information for teaching and outreach can be found at the
# Eirene homepage, http://gregoryhenselman.org/eirene.

__precompile__()

module Eirene

using Distances
using FileIO
using MultivariateStats
using SparseArrays
using LinearAlgebra
using Dates
using Statistics
using DelimitedFiles
using CSV
using Hungarian # Wasserstein distances

export eirene,
       eirenefilepath,
       barcode,
       classrep,
       wasserstein_distance # in wassterstein_distances.jl

include("wasserstein_distances.jl")
include("simplicial_constructions.jl")
include("schur_complements.jl")
include("inversion.jl")
include("multiplication.jl")
include("transposition.jl")
include("barcode_utils.jl")

##########################################################################################

####	CHAIN OPERATIONS

##########################################################################################

#=
notes on morselu!
- the output array tid has size equal to the number of columns of the
input array, NOT the column-rank of the input array
-	the columns of M must be ordered by grain (ascending)
-	the first (rank of M) elements of tlab index the complete set
of nonzero columns in the reduced matrix
=#
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
    schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
              rowlab,collab,
              Jprows,Jpcols,numjunpairs,
              Sprows,Spcols,numsenpairs,
              comprows,compcols,
              Trv,Tcp,Srv,Scp;
              updatetransform = storetransform,
              verbose = verbose
              )
    maxnz = max(maxnz,Mcp[Mn[1]+1])

    if verbose
        println("first shurr finished")
        if Mn[1]>0
            println([Mcp[Mn[1]] "nnz(M)" Mm[1] "Mm" Mn[1] "Mn" length(Mrv) "length(Mrv)"])
        else
            println("Mn = 0")
        end
    end
    #gc()
    rowfilt = Array{Int}(undef,length(comprows)); colfilt = Array{Int}(undef,length(compcols))
    counter = 0
    while Mcp[Mn[1]+1]>1
        if verbose
            println("starting block operation $(counter)")
        end
        counter+=1
        for i = 1:Mm[1]
            rowfilt[i] = Mrowgrain[rowlab[i]]
        end
        for j = 1:Mn[1]
            colfilt[j] = Mcolgrain[collab[j]]
        end
        getPairsLightWrite2!(Mrv,Mcp,rowfilt,colfilt,Mm[1],Mn[1],Jprows,Jpcols,numjunpairs,verbose=verbose)
        comprows = intervalcomplementuniqueunsortedinput(Jprows[1:numjunpairs[1]],Mm[1])
        compcols = intervalcomplementuniqueunsortedinput(Jpcols[1:numjunpairs[1]],Mn[1])
        if diagnostic
            numsenpairsOLD = numsenpairs[1]
        end

        schurit4!(		Mrv,Mcp,Mm,Mn,Mn0,
                  rowlab,collab,
                  Jprows,Jpcols,numjunpairs,
                  Sprows,Spcols,numsenpairs,
                  comprows,compcols,
                  Trv,Tcp,Srv,Scp;
                  updatetransform = storetransform,
                  verbose = verbose
                  )
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
        maxsd::Integer;
        record="cyclerep",
        verbose=false)

    farfaces::Array{Array{Int,1},1}
    firstv::Array{Array{Int,1},1}
    prepairs::Array{Array{Int,1},1}
    grain::Array{Array{Int,1},1}

    if record == "all" || record == "cyclerep"
        storetransform = true
    else
        storetransform = false
    end

    m = length(firstv[1])-1

    trv::Array{Array{Int,1},1}			=Array{Array{Int,1},1}(undef,maxsd+1);
    tcp::Array{Array{Int,1},1}			=Array{Array{Int,1},1}(undef,maxsd+1);
    phi::Array{Array{Int,1},1}			=Array{Array{Int,1},1}(undef,maxsd+1);
    plo::Array{Array{Int,1},1}			=Array{Array{Int,1},1}(undef,maxsd+1);
    tid::Array{Array{Int,1},1}			=Array{Array{Int,1},1}(undef,maxsd+1);
    for i in [1,maxsd+1]
        tcp[i]  =[1];
        trv[i]				=Array{Int}(undef,0)
        tid[i]				=Array{Int}(undef,0)
        phi[i]				=Array{Int}(undef,0)
        plo[i]				=Array{Int}(undef,0)
    end

    maxnzs =zeros(Int,maxsd+1);

    for sd = 2:maxsd
        if sd > length(farfaces)
            continue
        elseif sd>2
            lowbasisnames = phi[sd-1]
        else
            lowbasisnames = Array{Int}(undef,0)
        end
        Mrv::Array{Int,1},
        Mcp::Array{Int,1},
        lowlab::Array{Int,1},
        higlab::Array{Int,1},
        Mm =
        filteredmatrixfromfarfaces(farfaces,firstv,prepairs,grain,sd,lowbasisnames;verbose=verbose)
        lowlabtemp = convert(Array{Int,1},1:length(lowlab))
        higlabtemp = convert(Array{Int,1},1:length(higlab))
        higfilttemp = grain[sd][higlab]
        lowfilttemp = grain[sd-1][lowlab]
        if verbose
            println("Constructed Morse boundary operator, columns indexed by cells of dimension $(sd-1)")
        end
        pplow = convert(Array,length(prepairs[sd]):-1:1)
        pphig = convert(Array,length(prepairs[sd]):-1:1)

        # NB: It is critical that the columns of the input array
        # be ordered according to filtration; in particular, the entries of
        # lowfiltemp should increase monotonically
        Srv,Scp,Sphigs,Splows,tlab,maxnz =
        morselu!(
                 Mrv,
                 higfilttemp,
                 Mcp,
                 lowfilttemp,
                 lowlabtemp,
                 higlabtemp,
                 pplow,
                 pphig,
                 Mm,
                 storetransform = storetransform,
                 verbose = verbose)
        trv[sd] = lowlab[Srv]
        tcp[sd] = Scp
        tid[sd] = lowlab[tlab]
        plo[sd] = lowlab[Splows]
        phi[sd] = higlab[Sphigs]
        maxnzs[sd]= maxnz
    end
    return trv,tcp,plo,phi,tid,maxnzs
end

function persistf2!(D::Dict;
        maxsd=0,
        dictionaryoutput::Bool = true,
        verbose::Bool = false,
        record = "cyclerep")

    farfaces = D["farfaces"]
    firstv = D["firstv"]
    prepairs = D["prepairs"]
    grain = D["grain"]
    if maxsd == 0 maxsd = length(farfaces) - 1 end

    trv,tcp,plo,phi,tid,maxnzs =
    persistf2_core_vr(farfaces,firstv,prepairs,grain,maxsd::Integer;record=record,verbose = verbose)
    if dictionaryoutput == true
        D["trv"] = trv
        D["tcp"] = tcp
        D["tid"] = tid
        D["plo"] = plo
        D["phi"] = phi
        D["maxnz"] = maxnzs
        return D
    else
        return trv,tcp,plo,phi,tid
    end
end

function persistf2vr(s, maxsd;
        model = "vr",
        entryformat = "textfile",
        minrad			= -Inf,
        maxrad			= Inf,
        numrad			= Inf,
        nodrad  = [],
        fastop			= true,
        vscale			= "diagonal",
        pointlabels = [],
        verbose = false,
        record = "cyclerep")

    inputisfile = false
    filename = "user input julia array"
    if model == "pc" || model == "perseus_brips"
        pc = "genera"
    else
        pc = "n/a"
    end

    input = Dict(
                 "model"			=> model,
                 "genera"		=> copy(s),
                 "pc"			=> pc,
                 "source"		=> filename,
                 "maxdim" => maxsd-2,
                 "maxrad"		=> maxrad,
                 "minrad"		=> minrad,
                 "numrad"		=> numrad,
                 "nodrad"		=> nodrad,
                 "fastop"		=> fastop,
                 "record"	 => record,
                 )

    # Determine the number of points
    if model == "pc"
        numpoints = size(s,2)
    elseif model == "vr"
        numpoints = size(s,1)
        if !issymmetric(s)
            print("It appears the input matrix is not symmetric.  Only symmetric distance matrices are accepted when the <model> keyword argument has value \"clique\".")
            return
        end
    else
        print("Keyword argument <model> must be \"vr\", \"pc\", or \"cell\".  Please see documentaiton for further details.")
    end

    #### type the matrix
    s = convert(Array{Float64,2},s)
    if model == "pc"
        d = Distances.pairwise(Euclidean(), s, dims=2)
        if !isempty(nodrad)
            for i = 1:numpoints
                d[i,i] = nodrad[i]
            end
        end
    elseif model == "vr"
        d = convert(Array{Float64,2},s)
    end

    if fastop
        maxrad_alt = minimum(maximum(d,dims=1))
        maxrad_alt  = min(maxrad_alt,maxrad)
    else
        maxrad_alt = maxrad
    end

    d = minmaxceil(d, minrad = minrad, maxrad = maxrad_alt, numrad = numrad)

    # vfilt = diag(t)
    # recall that t will have the same order (NOT inverted) as d
    # <trueordercanonicalform> is a bit like <integersinsameorder>, just valid for floating point inputs, and with a bit more data in the output
    t,ocg2rad = trueordercanonicalform(d,factor=true)

    t = (1+maximum(t)).-t
    ocg2rad = reverse(ocg2rad,dims=1)

    if any(d.>maxrad_alt)
        t = t.-1
        deleteat!(ocg2rad,1)
    end

    vertices2keep = findall(diag(t).!=0)  # this step is necessary in order to cover the case where some vertices never enter the filtration
    t = t[vertices2keep,vertices2keep]

    #### Build the complex
    D = buildcomplex3(t, maxsd; verbose=verbose)
    D["ocg2rad"]=ocg2rad

    #### Compute persistence
    persistf2!(D; record=record)

    #### Store input data
    D["input"] = input
    # this covers the case where some vertices never enter the filtration
    D["nvl2ovl"] = vertices2keep[D["nvl2ovl"]]

    #### Store generators
    #gc()
    if record == "all" || record == "cyclerep"
        unpack!(D)
    end
    #gc()
    if record == "cyclerep"
        delete!(D,"trv")
        delete!(D,"tcp")
        delete!(D,"Lrv")
        delete!(D,"Lcp")
        delete!(D,"Rrv")
        delete!(D,"Rcp")
        delete!(D,"Lirv")
        delete!(D,"Licp")
        delete!(D,"prepairs")
    end
    
    return D
end

function boundarylimit_Lionly(trv,tcp,tid,sd;verbose=false)
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

function boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
    numl = length(farfaces[sd-1])
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])
    brv,bcp = maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
    numl = length(cp[sd-1])-1
    nump = length(phi[sd])
    numnlpl = numl-length(phi[sd-1])
    brv,bcp		= maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
    Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
    return boundarylimit_core(brv,bcp,trv[sd],tcp[sd],tid[sd],numl,nump,numnlpl)
end

function boundarylimit(D::Dict,sd)
    # special case sd = 1 may possibly be degenerate
    trv = D["trv"];tcp=D["tcp"];plo=D["plo"];phi=D["phi"];tid=D["tid"]
    if haskey(D,"farfaces")
        farfaces = D["farfaces"];firstv = D["firstv"]
        return boundarylimit_simplex(farfaces,firstv,trv,tcp,plo,phi,tid,sd;verbose=false)
    else
        rv = D["rv"];cp = D["cp"]
        # Lrv,Lcp,Lirv,Licp,Rrv,Rcp = boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
        return boundarylimit_cell(rv,cp,trv,tcp,plo,phi,tid,sd;verbose=false)
    end
end

function maxnonsingblock_simplex(farfaces,firstv,plo,phi,tid,sd;verbose=false)
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

function maxnonsingblock_cell(rv,cp,plo,phi,tid,sd;verbose=false)
    numl = length(cp[sd-1])-1
    nump = length(phi[sd])
    numnlpl = numl-length(plo[sd-1])

    if isempty(phi[sd])
        brv = zeros(Int,0)
        bcp =   ones(Int,numnlpl+1)
        return brv,bcp
    end

    lowtranslator = zeros(Int,numl)
    lowtranslator[tid[sd]] = 1:numnlpl
    dummy0 = zeros(Int,0)

    brv,dummy1,bcp,dummy2 = stackedsubmatrices(rv[sd],cp[sd],tid[sd],dummy0,phi[sd],numl)
    brv .= lowtranslator[brv]
    append!(bcp,fill(bcp[end],numnlpl-nump)) # <-- note there's no +1 b/c bcp is already 1 elt. longer than the # of columns
    return brv,bcp
end

function unpack!(D::Dict)
    maxsd = D["input"]["maxdim"]+2

    Lirv = Array{Array{Int,1},1}(undef,maxsd);  Licp = Array{Array{Int,1},1}(undef,maxsd)
    Lrv  = Array{Array{Int,1},1}(undef,maxsd);  Lcp  = Array{Array{Int,1},1}(undef,maxsd)
    Rrv  = Array{Array{Int,1},1}(undef,maxsd);  Rcp  = Array{Array{Int,1},1}(undef,maxsd)

    if D["input"]["record"] == "all"
        N = maxsd
    else
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
    end

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


##########################################################################################

####	COPY, SHIFT, INDEX, AND SLICE

##########################################################################################

cran(A::SparseMatrixCSC, j) = A.colptr[j]:(A.colptr[j+1]-1)
cran(colptr::Array, j::Int) = colptr[j]:(colptr[j+1]-1)
function cran(colptr::Array, J::Array{Int,1})
    m = nval(colptr,J)
    v = zeros(Int,m)
    c = 0
    for p=1:length(J)
        k = nval(colptr,J[p])
        v[c+1:c+k]=cran(colptr,J[p])
        c += k
    end
    return v
end
function cran(colptr::Array, J::UnitRange{Int})
    m = nval(colptr,J)
    v = zeros(Int,m)
    c = 0
    for p=1:length(J)
        k = nval(colptr,J[p])
        v[c+1:c+k]=cran(colptr,J[p])
        c += k
    end
    return v
end
cran(colptr::UnitRange, j) = colptr[j]

"""
For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column `j`.
"""
nval(colptr, j::Int) = colptr[j+1]-colptr[j]
"""
For a column sparse matrix with column pattern `colptr`, counts the number of values stored for column j, for each j in `J`.
"""
function nval(colptr, J)
    c = 0
    for p = 1:length(J)
        c += nval(colptr,J[p])
    end
    return c
end

crows(A::SparseMatrixCSC, j) = A.rowval[cran(A, j)]
crows(colptr::Array,rowval::Array, j) = rowval[cran(colptr, j)]

function findcol(cp, k)
    i = 1
    while cp[i]<=k
        i+=1
    end
    return i-1
end

function extend!(x::Array{Tv,1},n::Integer) where Tv
    if length(x)<n
        append!(x,Array{Tv}(undef,n-length(x)))
    end
end

function copycolumnsubmatrix(Arv::Array{Tv,1},Acp,columnindices) where Tv<:Integer
    allocationspace = 0
    for j in columnindices
        allocationspace+= Acp[j+1]-Acp[j]
    end
    Brv = Array{Tv}(undef,allocationspace)
    Bcp = Array{Int}(undef,length(columnindices)+1)
    Bcp[1]=1
    for jp = 1:length(columnindices)
        j = columnindices[jp]
        Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
        Brv[Bcp[jp]:(Bcp[jp+1]-1)]=Arv[Acp[j]:(Acp[j+1]-1)]
    end
    return Brv,Bcp
end



function extendcolumnlight!(rowval::Array{Ti,1},colptr::Array{Ti,1},v::Array{Ti},k::Ti,growthincrement::Ti) where Ti
    r = rowval
    c = colptr
    startpoint = copy(c[k+1])
    c[k+1]=c[k]+length(v)
    if length(r)<c[k+1]-1
        append!(r,Array{Int}(undef,max(growthincrement,length(v))))
    end
    r[startpoint:(c[k+1]-1)]=v
end

# colsinorder must be in sorted order
function supportedmatrix!(Mrowval::Array{Tv},Mcolptr::Array{Tv,1},rows1,colsinorder,Mm::Tv) where Tv<:Integer
    n = length(colsinorder)
    suppcol1 = falses(Mm)
    suppcol1[rows1].=true
    cpHolder = 1
    nz1 = 0
    for jp = 1:n
        for ip in cran(Mcolptr,colsinorder[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1+=1
                Mrowval[nz1]=i
            end
        end
        Mcolptr[jp]=cpHolder
        cpHolder = nz1+1
    end
    Mcolptr[n+1] = cpHolder
    deleteat!(Mcolptr,(n+2):length(Mcolptr))
    deleteat!(Mrowval,Mcolptr[end]:length(Mrowval))
end

function stackedsubmatrices(
        Mrowval::Array{Int,1},
        Mcolptr::Array{Int,1},
        rows1::Array{Int,1},
        rows2::Array{Int,1},
        cols::Array{Int,1},
        Mm::Int)

    n = length(cols)
    suppcol1 = falses(Mm)
    suppcol2 = falses(Mm)
    suppcol1[rows1].=true
    suppcol2[rows2].=true
    nz1 = 0; nz2 = 0
    for jp = 1:n
        for ip in cran(Mcolptr,cols[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1+=1
            elseif suppcol2[i]
                nz2+=1
            end
        end
    end
    rv1 = Array{Tv}(undef,nz1)
    rv2 = Array{Tv}(undef,nz2)
    cp1 = Array{Tv}(undef,n+1); cp1[1]=1
    cp2 = Array{Tv}(undef,n+1); cp2[1]=1
    nz1 = 0; nz2 = 0
    for jp = 1:n
        for ip in cran(Mcolptr,cols[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1+=1
                rv1[nz1]=i
            elseif suppcol2[i]
                nz2+=1
                rv2[nz2]=i
            end
        end
        cp1[jp+1] = nz1+1
        cp2[jp+1] = nz2+1
    end
    return rv1,rv2,cp1,cp2
end

##########################################################################################

####	MATRIX WEIGHTS AND FORMATTING

##########################################################################################

#=
The output of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil(N;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)
    S = copy(N) # NB: this is probably an important step; several pernicious problems in the development of rounding procedures turned out to be rooted in accidental rewriting of matrix entries
    return minmaxceil!(S;minrad=minrad,maxrad=maxrad,numrad=numrad)
end

#=
NB THIS FUNCTION MODIFIES ITS INPUT

The result of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
=#
function minmaxceil!(S;minrad=minimum(N),maxrad=maximum(N),numrad=Inf)

    if (minrad == Inf) || (maxrad == -Inf)
        return fill(Inf,size(S)...)
    end

    if minrad == "minedge"
        minrad = minimum(offdiagmin(S))
    end

    S[S.<minrad] .= minrad
    S[S.>maxrad] .= Inf

    fi = (LinearIndices(S))[findall(isfinite,S)] # stands for finite indices
    fv = S[fi]

    if isempty(fv)
        return S
    end

    if minrad == -Inf
        minrad =  minimum(fv)
    end
    if maxrad == Inf
        maxrad = maximum(fv)
    end
    if numrad == 1
        S[fi]	.= maxrad
        return S
    end
    if numrad == Inf
        return S
    end

    ran = range(minrad,maxrad,length=numrad)
    ran = Array{Float64}(ran)
    ran[end] = maxrad

    fvr = ceil2grid(fv,ran)

    S[fi]		= fvr

    return S
end

# ran should be an array in sorted order, with ran[end] >= maximum(A)
function ceil2grid(A,ran)
    if ran == "all"
        return A
    end
    B = copy(A)
    for j = 1:length(A)
        post = 1
        while ran[post] < A[j]
            post+=1
        end
        B[j] = ran[post]
    end
    return B
end

function checkceil2grid(numits)
    for p =	1:numits
        A = rand(70)
        ran = sort(rand(20))
        ran = ran*(1/maximum(ran)) # this guarantees that maximum(ran) > maximum(A)
        crct = crosscheckceil2grid(A,ran)
        if !crct
            println("error: please check checkceil2grid")
        end
    end
    return []
end

function crosscheckceil2grid(A, ran)
    B = ceil2grid(A, ran)
    for j = 1:length(A)
        q = findfirst(ran .> A[j])
        val = ran[q]
        if B[j] != val
            println("error: please check <checkceil2grid>")
            return false
        end
    end
    return true
end

function offdiagmin(S, i)
    if i == 1
        return(minimum(S[2:end,i]))
    elseif i == size(S,1)
        return minimum(S[1:end-1,i])
    else
        return min(minimum(S[1:i-1,i]),minimum(S[i+1:end,i]))
    end
end

function checkoffdiagmin(numits)
    for p = 1:numits
        S = rand(50,50)
        for q = 1:50
            checkval = minimum(deleteat!(S[:,q],q))
            if checkval != offdiagmin(S,q)
                return S
            end
        end
    end
    return zeros(Int,0)
end

function offdiagmean(S; defaultvalue=[])
    m,n = size(S)
    if m != n
        println("error in <offdiagmean>: input matrix should be square")
    end
    if isempty(defaultvalue)
        println("warning: no defaulvalue was set for <offdiagmin>; this parameter has been set to 0")
        defaultvalue = 0
    end
    if m == 1
        return defaultvalue
    end
    mu = zeros(m)
    for j = 1:m
        v = S[1:(j-1),j]
        u = S[(j+1):m,j]
        mu[j]  = mean(vcat(v[:],u[:]))
    end
    return mu
end

# NB: Assumes that the input array S has only finite entries.
# NB: The value for keyword <numrad> must be either a positive integer or Inf
function ordercanonicalform_3(
        S::Array{Tv};
        maxrad=Inf,
        minrad=-Inf,
        numrad=Inf,
        vscale="default",
        fastop::Bool=true) where Tv

    if size(S,1) == 0
        ocf = zeros(Int,0,0)
        ocg2rad		= zeros(Float64,0)
        return ocf,ocg2rad
    end

    # Format input
    symmat		= convert(Array{Float64},copy(S))
    m = size(symmat,1)

    if vscale == "default"
        for i = 1:m
            symmat[i,i] = minimum(symmat[:,i])
        end
    elseif typeof(vscale) <: Array
        vscale = convert(Array{Float64},copy(vscale))
        if length(vscale) != m
            print("Error: keyword <vscale> must take a value of <defualt>, <diagonal>, or <v>, where v is a vector of length equal to the number of vertices in the complex.")
            return
        end
        for i=1:m
            if offdiagmin(symmat,i) < vscale[i]
                print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
                return
            else
                symmat[i,i] = vscale[i]
            end
        end
    elseif vscale == "diagonal"
        vscale		=	Array{Float64,1}(m)
        for i=1:m
            vscale[i]	=	symmat[i,i]
        end
        # the following is in prnciple unnecessary, but it simplifies rounding
        for i=1:m
            if offdiagmin(symmat,i) < vscale[i]
                print("Error: the \"birth time\" assigned a vertex by keyword argument <vscale> may be no greater than that of any incident edge.  The value assigned to vertex $(i) by <vscale> is $(vscale[i]), and i is incident to an edge with birth time $(offdiagmin(symmat,i)).")
                return
            end
        end
    end

    # Deterime the public maxrad
    if maxrad == Inf
        publicmax = maximum(symmat)
    else
        publicmax = copy(maxrad)
    end

    # Deterime the public minrad
    publicmin	= minimum(symmat)
    publicmin	= max(publicmin,minrad)

    # It's important that this precede the other cases
    if publicmax < publicmin
        return zeros(Int,m,m),Array{Float64}(undef,0)
    end

    # This covers all remaining cases where numrad ==1
    if numrad == 1
        privatemax = minimum(maximum(symmat,1))
        ocf = zeros(Int,m,m)
        if fastop && publicmax >= privatemax
            index = findfirst(maximum(symmat,1),privatemax)
            ocf[index,:]=1
            ocf[:,index]=1
            for i = 1:m
                ocf[i,i]=1
            end
        else
            ocf[symmat.<=publicmax]=1
        end
        return ocf,[publicmax]
    end

    # If necessary, determine step size.  Recall we have already treated every case where numrad == 1.
    if numrad < Inf
        alpha = (publicmax-publicmin)/(numrad-1)
    elseif numrad == Inf
        alpha = Inf
    end

    # If one stops early, determine when to stop
    if fastop
        privatemax = minimum(maximum(symmat,1))
        privatemax = min(privatemax,publicmax)
        if numrad < Inf
            post = publicmin
            stepcounter = 1
            while post < privatemax
                stepcounter+=1
                if stepcounter == numrad
                    post = publicmax # must take this rather cumbersome step on account of numerical error.
                else
                    post+=alpha
                end
            end
            privatemax = post
        end
    else
        privatemax = publicmax
    end

    # Extract sortperm
    p = sortperm(vec(symmat),alg=MergeSort)

    # Compute the ocf
    val						= publicmin
    ocg2rad = Array{Float64}(undef,binom(m,2)+m) #the plus 1 covers values taken from the diagonal
    ocg2rad[1]				= val
    post = 1
    exceededmax = false
    ocf						= fill(-1,m,m) #will reverse the order on this after it's been filled
    stepcounter = 1
    for i = 1:length(p)
        if symmat[p[i]] <= val
            ocf[p[i]] = post
        else
            if numrad == Inf
                val = symmat[p[i]]
            else
                if symmat[p[i]] == Inf
                    val = Inf
                else
                    while symmat[p[i]] > val
                        stepcounter+=1
                        if stepcounter == numrad
                            val = publicmax # must take this rather cumbersome final step b/c of numerical error, since o/w it can and does happen that the (numrad)th grain value fails to equal publicmax
                        else
                            val+=alpha
                        end
                    end
                end
            end
            post+=1
            ocf[p[i]] = post
            ocg2rad[post]	= val
        end
        if val > privatemax
            ocf[p[i:end]] = post
            exceededmax		= true
            break
        end
    end
    if exceededmax
        cutoff = post
    else
        cutoff = post+1
    end
    deleteat!(ocg2rad,cutoff:length(ocg2rad))
    ocg2rad = reverse(ocg2rad,dims=1)
    ocf = cutoff - ocf
    return ocf,ocg2rad # additional outputs for diagnostic purposes -> #,privatemax,S,maxrad,publicmax,publicmin,maxrad
end

function ordercanonicalform_4(d;
        minrad=-Inf,
        maxrad=Inf,
        numrad=Inf, # note we have already performed the necessary rounding by this point
        fastop=true,
        vscale="diagonal")

    # round as necessary
    if minrad == "minedge" minrad = minimum(offdiagmin(d)) end
    d = minmaxceil(d,minrad=minrad,maxrad=maxrad,numrad=numrad)
    d[d.>maxrad] = maxrad+1; # we make this reassignment b/c ordercanonicalform_3 takes only arguments with finite entries

    (t,ocg2rad) = ordercanonicalform_3(d;
                                       minrad=minrad,
                                       maxrad=maxrad,
                                       numrad=Inf, # note we have already performed the necessary rounding by this point
                                       fastop=fastop,
                                       vscale=vscale)
    return t, ocg2rad
end

function ordercanonicalform(S::Array{Tv};
        minrad=-Inf,
        maxrad=Inf,
        numrad=Inf,
        fastop::Bool=true) where Tv

    symmat_float = convert(Array{Float64},copy(S))
    symmat = copy(S)
    m = size(symmat,1);
    convert(Tv,minrad)
    convert(Tv,maxrad)

    effectivemin = -maxrad
    effectivemax = -minrad
    symmat = -symmat

    for i = 1:m
        symmat[i,i] = Inf
    end
    if fastop
        maxmin = -Inf
        for j = 1:m
            holdmin = minimum(symmat[:,j])
            if holdmin > maxmin
                maxmin = holdmin
            end
        end
        if maxmin > effectivemin
            effectivemin = maxmin
        end
    end
    if numrad == 1
        for i = 1:m
            for j = (i+1):m
                sij = symmat[i,j]
                if sij>effectivemin
                    symmat[i,j] = 1
                    symmat[j,i] = 1
                else
                    symmat[i,j] = 0
                    symmat[j,i] = 0
                end
            end
        end
        for i = 1:m symmat[i,i] = 0 end
        ocg2rad = [1]
        return round.(Int,symmat),ocg2rad
    end
    numfilt = binom(m,2)
    for i = 1:m
        for j = (i+1):m
            sij = symmat[i,j]
            if sij >= effectivemax
                symmat[i,j] = Inf
                symmat[j,i] = Inf
            elseif sij < effectivemin
                # note that loose inequality
                # here could cause some bars
                # that disappear at the last
                # grain to appear to
                # live forever
                symmat[i,j] = -Inf
                symmat[j,i] = -Inf
            end
        end
    end
    for i = 1:m
        symmat[i,i] = -Inf
    end
    ocg2rad = zeros(Float64,binom(m,2))
    p = sortperm(symmat[:],alg=MergeSort)
    ordervalue = 0
    floatingvalue = symmat[p[1]]
    for i = 1:m^2
        ii = p[i]
        if symmat[ii] == floatingvalue
            symmat[ii] = ordervalue
        else
            ordervalue+=1
            floatingvalue = symmat[ii]
            ocg2rad[ordervalue]=symmat_float[ii]
            symmat[ii]=ordervalue
        end
    end
    deleteat!(ocg2rad,(ordervalue+1):binom(m,2))
    return round.(Int, symmat), ocg2rad
end

function trueordercanonicalform(M; rev=false, firstval=1, factor=false)
    m = length(M)
    perm = sortperm(M[:], rev=rev, alg=MergeSort)
    oca = zeros(Int, size(M)...)

    if m ==	0
        return zeros(Int,0), zeros(Int,0)
    end

    if factor
        numvals = 1
        post = 1
        for p = 1:m
            if M[perm[p]] != M[perm[post]]
                post = p
                numvals += 1
            end
        end
        oca2rad = Array{Float64}(undef,numvals)
        oca2rad[1] = M[perm[1]]
    end

    post = [1]
    k = [1]
    trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
    if factor
        return oca, oca2rad
    else
        return oca
    end
end

function trueordercanonicalform_subr!(M,perm,oca,oca2rad,post,k,m,factor)
    for p  = 1:m
        if M[perm[p]] != M[perm[post[1]]]
            post[1] = p
            k[1] =  k[1]+1
            if factor
                oca2rad[k[1]] = M[perm[post[1]]]
            end
        end
        oca[perm[p]] = k[1]
    end
end

function checktrueordercanonicalform(numits)
    for p = 1:numits
        m = rand(50:300,1)
        m = m[1]
        x = rand(m,m)
        if isodd(p)
            x = x+x';
        end
        ocf,val = trueordercanonicalform(x,factor=true)
        check1 = x == val[ocf]
        check2 = val == sort(unique(x))

        k = rand(1:100,1)
        k = k[1]
        ocau,valu	=   trueordercanonicalform(x,firstval=1,factor=true,rev=false)
        ocad,vald	=   trueordercanonicalform(x,firstval=1,factor=true,rev=true)
        check3 = ocau == (maximum(ocad)+1).-ocad
        check4 = x == valu[ocau]
        check5 = x == vald[ocad]
        check6 = length(valu) == maximum(ocau)
        check7 = length(vald) == maximum(ocad)

        if !all([check1, check2, check3, check4, check5, check6, check7])
            println("error: please check trueordercanonicalform")
            println([check1, check2, check3, check4, check5, check6, check7])
            return x
        end
    end
    return []
end

function checkoffdiagmean(numits)
    for p = 1:numits
        numpts = rand(50:100,1)
        numpts = numpts[1]
        S = rand(numpts,numpts)
        T = copy(S)
        mu = offdiagmean(T,defaultvalue = 0)
        u = zeros(numpts)
        for i = 1:numpts
            ran =   setdiff(1:numpts,i)
            u[i]= mean(S[ran,i])
        end
        if u[:] != mu[:]
            println("error: please check <offdiagmean>")
            return S,T,mu,u
        end
    end
    return []
end

function getstarweights(symmat)
    m = size(symmat,1)
    w = zeros(Int,m)
    getstartweights_subr2(symmat::Array{Int,2},w::Array{Int,1},m::Int)
    return w
end

function getstartweights_subr2(symmat::Array{Int,2},w::Array{Int,1},m::Int)
    s = copy(symmat)
    for i = 1:m
        s[i,i] = 0
    end
    l = Array{Int}(undef,m)
    lDown = Array{Int}(undef,m)
    val = Array{Array{Int,1}}(undef,m)
    supp = Array{Array{Int,1}}(undef,m)
    suppDown = Array{Array{Int,1}}(undef,m)
    for i = 1:m
        supp[i] = findall(!iszero,s[:,i])
        suppDown[i] = i .+findall(!iszero,s[(i+1):end,i])
        val[i] = s[supp[i],i]
        l[i] = length(supp[i])
        lDown[i] = length(suppDown[i])
    end

    for i = 1:m
        Si = supp[i]
        Vi = val[i]
        for jp = 1:l[i]
            j = Si[jp]
            dij = Vi[jp]
            Sj = suppDown[j]
            Vj = val[j]
            for kp = 1:lDown[j]
                k = Sj[kp]
                if k == i
                    continue
                end
                dkj = Vj[kp]
                dki = s[i,k]
                if dki >= dkj && dij >= dkj
                    w[i]+=1
                end
            end
        end
    end
    return w
end

##########################################################################################

####	COMBINATIONS, PERMUTATIONS, AND SET OPERATIONS

##########################################################################################


function integersinsameorder!(v::Array{Int,1},maxradue::Int)
    m = length(v)
    x = zeros(Int,maxradue)
    for i = 1:m
        x[v[i]]+=1
    end
    y = Array{Int}(undef,maxradue+1)
    y[1] = 1
    for i = 1:maxradue
        y[i+1]=y[i]+x[i]
    end
    for i = 1:length(v)
        u = v[i]
        v[i] = y[u]
        y[u]+=1
    end
    return v
end

function integersinsameorder(v::Array{Int,1})
    # Returns the permutation z on {1,...,length(v)} such z[i]<z[j] iff either
    # (a) v[i] < v[j], or
    # (b) v[i] = v[j] and i < j
    if isempty(v)
        z = Array{Int}(undef,0)
        return z
    else
        m = length(v)
        maxv = maximum(v)
        minv = minimum(v)
        minv = minv-1;
        x = zeros(Int,maxv-minv)
        z = Array{Int}(undef,length(v))
        for i = 1:m
            x[v[i]-minv]+=1
        end
        prevsum = 1
        for i = 1:length(x)
            sum = prevsum + x[i]
            x[i] = prevsum
            prevsum = sum
        end
        for i = 1:m
            u = v[i]
            z[i] = x[u-minv]
            x[u-minv]+=1
        end
        return z
    end
end

function integersinsameorder!(v::Array{Int,1})
    # Replaces v with the permutation z on {1,...,length(v)} such that z[i]<z[j] iff either
    # (a) v[i] < v[j], or
    # (b) v[i] = v[j] and i < j
    if isempty(v)
        z = Array{Int}(undef,0)
        return z
    else
        m = length(v)
        maxv = maximum(v)
        minv = minimum(v)
        minv = minv-1;
        x = zeros(Int,maxv-minv)
        z = Array{Int}(undef,length(v))
        for i = 1:m
            x[v[i]-minv]+=1
        end
        prevsum = 1
        for i = 1:length(x)
            sum = prevsum + x[i]
            x[i] = prevsum
            prevsum = sum
        end
        for i = 1:m
            u = v[i]
            v[i] = x[u-minv]
            x[u-minv]+=1
        end
    end
end


#=
- In beta; should be compared with integersinsameorderbycolumn3.  See
/Users/greghenselman/Google Drive/GregDirectory/julia_gd/Julia/workshop/workshop_Oct2017.jl
-   Functionally equivalent to integersinsameorderbycolumn; returns a
permutation z on {1,...,length(v)} so that for all j
- cran(colptr,j) maps to cran(colptr,j), and
- crows(colptr,v[z],j) is an array in sorted order
=#
function integersinsameorderbycolumn2(v::Array{Int,1},colptr)
    numcols = length(colptr)-1
    m = length(v)
    v = v.-(minimum(v)-1)
    x = zeros(Int,maximum(v))
    z = Array{Int}(undef,length(v))
    for j = 1:numcols
        if colptr[j] == colptr[j+1]
            continue
        end
        for i = colptr[j]:(colptr[j+1]-1)
            x[v[i]]+=1
        end
        maxv = v[colptr[j]];   minv = maxv
        for i = (colptr[j]+1):(colptr[j+1]-1)
            if v[i] > maxv
                maxv = v[i]
            elseif v[i] < minv
                minv = v[i]
            end
        end
        prevsum = colptr[j]
        for i = minv:maxv
            sum = prevsum + x[i]
            x[i] = prevsum
            prevsum = sum
        end
        for i = colptr[j]:(colptr[j+1]-1)
            u = v[i]
            z[i] = x[u]
            x[u]+=1
        end
        for i = minv:maxv
            x[i] = 0
        end
    end
    return z
end

##########################################################################################

####	USER-FRIENDLY UTILITIES

##########################################################################################

function boundarymatrix(C; dim=1, rows="a", cols="a")
    crr = complexrank(C,dim=dim-1)
    crc = complexrank(C,dim=dim)
    if rows == "a"
        rows = 1:crr#complexrank(C,dim=dim-1)
    end
    if cols == "a"
        cols = 1:crc#complexrank(C,dim=dim)
    end
    if empteval(maximum,cols,0) > crc
        println("error: keyword argument <cols> contains an integer greater than the rank of the complex in dimension <dim>")
        return
    elseif empteval(maximum,rows,0) > crc
        print("error: keyword argument <rows> contains an integer greater than the rank of the complex in dimension <dim-1>")
        return
    end
    if isempty(rows) || isempty(cols)
        rv = zeros(Int,0)
        cp = ones(Int,length(cols)+1)
        return rv,cp
    end
    ncols = length(cols)
    nrows = length(rows)
    sd = dim+1;
    if haskey(C,"farfaces")
        rv = ff2aflight(C,dim+1,cols)
        rv  = reshape(rv,length(rv))
        cp  = convert(Array{Int,1},1:sd:(ncols*sd+1))
        cols = 1:ncols
    else
        rv = C["rv"][sd]
        cp = C["cp"][sd]
    end
    rv,dummrv,cp,dummycp = stackedsubmatrices(
                                              rv,
                                              cp,
                                              rows,
                                              Array{Int}(undef,0),
                                              cols,
                                              max(empteval(maximum,rows,0),empteval(maximum,rv,0))
                                              )
    return rv,cp
end

function boundarymatrices(C::Dict)
    haskey(C,"farfaces") ? ff2complex(C["farfaces"], C["firstv"]) : C["rv"], C["cp"]
end

function classrep(D::Dict; dim=1, class=1, format="vertex x simplex")
    if any(class .> nnzbars(D,dim=dim))
        print("error: the value for keyword argument <class> has an integer greater than the number of nonzero bars in the specified dimension")
        return
    elseif !(0 <= dim <= D["input"]["maxdim"])
        print("error: barcodes were not computed in the specified dimension")
        return
    end
    if !haskey(D, "farfaces")
        format = "index"
    end

    if format == "vertex x simplex"
        return classrep_faces(D, dim=dim, class=class)
    elseif format == "vertex"
        return classrep_vertices(D, dim=dim, class=class)
    elseif format == "index"
        return classrep_cells(D, dim=dim, class=class)
    end
end

function classrep_cells(D::Dict; dim=1, class=1)
    sd = dim+2
    if haskey(D, "cyclerep")
        return D["cyclerep"][sd][class]
    else
        cyclename = barname2cyclename(D,class;dim=dim)
        return getcycle(D,sd,class)
    end
end

function classrep_faces(D::Dict; dim=1, class=1)
    sd = dim+2
    rep = classrep_cells(D, dim=dim, class=class)
    vrealization = vertexrealization(D, sd-1, rep)
    D["nvl2ovl"][vrealization]
end

function classrep_vertices(D::Dict; dim=1, class=1)
    sd = dim+2
    rep = classrep_cells(D, dim=dim, class=class)
    vertices = incidentverts(D, sd-1, rep)
    D["nvl2ovl"][vertices]
end

function cyclevertices(D::Dict; dim=1, cycle=1)
    sd = dim+2
    rep = getcycle(D, sd, cycle)
    vertices = incidentverts(D, sd-1, rep)
    D["nvl2ovl"][vertices]
end

function barcode(D::Dict; dim=1, ocf=false)
    if haskey(D,:perseusjlversion)
        return barcode_perseus(D,dim=dim)
    elseif haskey(D,:barcodes)
        return D[:barcodes][dim+1]
    elseif !haskey(D,"cp") & !haskey(D,"farfaces")
        print("unrecognized object:")
        display(D)
        return
    elseif dim > D["input"]["maxdim"]
        maxdim = D["input"]["maxdim"]
        println("error: barcodes were not computed in dimensions greater than $(maxdim).")
        return
    end
    sd = dim+2
    plo = D["plo"][sd]
    phi = D["phi"][sd]
    tid = D["tid"][sd]
    lg = D["grain"][sd-1]
    hg  = D["grain"][sd]

    mortalprimagrain = lg[plo]
    mortalultragrain = hg[phi]

    finind = findall(mortalprimagrain .!= mortalultragrain)
    numfin = length(finind)
    numinf = length(tid)-length(plo)

    mortalprimagrain = mortalprimagrain[finind]
    mortalultragrain = mortalultragrain[finind]

    mortalran = 1:numfin
    evergrran = numfin+1:numfin+numinf
    finran = 1:(2*numfin+numinf)
    # stands for finite range; these are the linear
    # indices of the barcode array that take finite
    # values
    evrgrbran = length(tid)-numinf+1:length(tid)
    # stands for evergreen birth range; this satisfies
    # tid[ebergrbran] = [array of evergreen classes]

    bc = zeros(Int, numfin+numinf, 2)
    bc[mortalran,1]    .= mortalprimagrain
    bc[mortalran,2]    .= mortalultragrain
    bc[evergrran,1] = lg[tid[evrgrbran]]

    if !ocf
        bcc = copy(bc)
        bc = Array{Float64}(bc)
        bc[finran] = D["ocg2rad"][bcc[finran]]
        bc[evergrran,2]    .= Inf
    else
        bc = length(D["ocg2rad"]).-bc
    end

    return bc
end


##########################################################################################

####	MISC

##########################################################################################


function binom(x, y)
    k = 1;
    for i = x:-1:(x-y+1)
        k = i*k
    end
    for i = 2:y
        k = k/i
    end
    k = convert(Int, k)
    return k
end


function addinteger!(v::Array{Tv,1},k::Int) where Tv
    for i = 1:length(v)
        v[i]+=k
    end
end

function sparseadjacencymatrix(A;inputis = "adjacencymatrix")
    if inputis == "adjacencymatrix"
        if size(A,1) != size(A,2)
            print("Error: unless the <inputis> keywork argument has value 'edges', the input array must be square.")
            return
        end
        m = size(A,1)
        rv = Array{Int}(undef,0)
        cp = zeros(Int,m+1)
        cp[1] = 1
        for i = 1:m
            adjverts = findall(!iszero,A[:,i])
            append!(rv,adjverts)
            cp[i+1] = cp[i]+length(adjverts)
        end
        return rv, cp
    elseif inputis == "edges"
        if size(A,1) != 2
            print("Error: when the <inputis> keywork argument has value 'edges', the input array must have exactly two rows.")
        end
        m = maximum(A)
        adjmat = falses(m,m)
        for i = 1:size(A,2)
            adjmat[A[1,i],A[2,i]]=true
        end
        adjmat[findall(transpose(adjmat))] .= true
        return sparseadjacencymatrix(adjmat)
    end
end

function hopdistance_sparse(rv,cp)
    m = length(cp)-1
    H = zeros(Int,m,m)
    for i = 1:m
        c = 0
        metnodes = falses(m)
        metnodes[i] = true
        fringenodes = falses(m)
        fringenodes[i] = true
        fringelist = [i]

        while !isempty(fringelist)
            c+=1
            for j in fringelist
                for k in crows(cp,rv,j)
                    if !metnodes[k]
                        metnodes[k] = true
                        fringenodes[k] = true
                        H[k,i] = c
                    end
                end
            end
            fringelist = findall(fringenodes)
            fringenodes[:].= false
        end
        H[.!metnodes,i].=m+1
    end
    return H
end

function hopdistance(rv,cp)
    return hopdistance_sparse(rv,cp)
end

function hopdistance(A;inputis = "fulladj")
    if inputis == "fulladj"
        rv,cp = sparseadjacencymatrix(A)
    elseif inputis == "edges"
        rv,cp = sparseadjacencymatrix(A,inputis="edges")
    else
        rv,cp = A
    end
    return hopdistance_sparse(rv,cp)
end


##########################################################################################

####	VIETORIS-RIPS CONSTRUCTION

##########################################################################################

function buildcomplex3(symmat::Array{Tv},maxsd; dictionaryoutput = true, verbose = false) where Tv

    grain = Array{Array{Int,1}}(undef,maxsd+1)
    farfaces = Array{Array{Int,1}}(undef,maxsd+1)
    prepairs = Array{Array{Int,1}}(undef,maxsd+1)
    firstv = Array{Array{Int,1}}(undef,maxsd+1)

    farfaces[maxsd+1] = Array{Int}(undef,0)
    firstv[maxsd+1] = ones(Int,1)
    grain[maxsd+1] = Array{Int}(undef,0)
    prepairs[maxsd+1] = Array{Int}(undef,0)

    m = size(symmat,1)
    w = vec(offdiagmean(symmat,defaultvalue=0)) # modified 02/12/2018

    vperm = sortperm(-w,alg=MergeSort)
    symmat = symmat[vperm,vperm]

    farfaces[1] = convert(Array,1:m)
    firstv[1] = convert(Array,1:(m+1))
    grain[1] = diag(symmat)
    prepairs[1] = Array{Int}(undef,0)

    r,c,z = generate2faces(symmat)
    farfaces[2] = r
    firstv[2] = c
    grain[2] = z
    prepairs[2] = Array{Int}(undef,0)

    if maxsd == 3
        generate3faces!(farfaces,firstv,grain,prepairs,m,symmat)
        if dictionaryoutput == true
            return Dict{String,Any}(
                                 "farfaces" => farfaces,
                                 "firstv" => firstv,
                                 "grain" => grain,
                                 "prepairs" => prepairs,
                                 "symmat" => symmat,
                                 "nvl2ovl"=>vperm)
        else
            return farfaces, firstv, grain, prepairs, symmat, vperm
        end
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
        zll= grain[sd-2]
        izfull = Array{Int}(undef,nll)
        r = Array{Int}(undef,startlength)
        z = Array{Int}(undef,startlength)
        c = Array{Int}(undef,m+1)
        c[1]=1
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
                    process_sd_lt_maxsd!(
                                         i::Int,j::Int,dij::Int,stepsize::Int,
                                         facecount::Array{Int,1},numpairs::Array{Int,1},
                                         jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                         r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                         izfull::Array{Int,1},ff2pv::Array{Int,1},pmhist::Array{Int,2},
                                        npsupp::BitArray{1})
                elseif sd == maxsd-1
                    process_sd_onelt_maxsd_1!(
                                              i::Int,j::Int,dij::Int,stepsize::Int,
                                              facecount::Array{Int,1},numpairs::Array{Int,1},
                                              jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                              r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                              izfull::Array{Int,1},ff2pv::Array{Int,1},pmhist::Array{Int,2},
                                             npsupp::BitArray{1})
                else
                    for l = 1:(i-1)
                        oldclaw[l] = minimum(symmat[l,[i,j]])
                    end
                    process_maxsd_one2i!(
                                         i::Int,j::Int,dij::Int,stepsize::Int,
                                         facecount::Array{Int,1},numpairs::Array{Int,1},
                                         jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                         r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                         oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
                                         rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
                                         izfull::Array{Int,1},ff2pv::Array{Int,1},
                                         pmhist::Array{Int,2},fpi::Array{Int,2},
                                        npsupp::BitArray{1})

                    process_maxsd_i2i!(
                                       i::Int,j::Int,dij::Int,stepsize::Int,
                                       facecount::Array{Int,1},numpairs::Array{Int,1},
                                       jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                       r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                       oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
                                       rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
                                       izfull::Array{Int,1},ff2pv::Array{Int,1},
                                       pmhist::Array{Int,2},fpi::Array{Int,2},
                                      npsupp::BitArray{1})

                    process_maxsd_i2j!(
                                       i::Int,j::Int,dij::Int,stepsize::Int,
                                       facecount::Array{Int,1},numpairs::Array{Int,1},
                                       jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                       r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                       oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
                                       rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
                                       izfull::Array{Int,1},ff2pv::Array{Int,1},
                                       pmhist::Array{Int,2},fpi::Array{Int,2},
                                      npsupp::BitArray{1})

                    process_maxsd_j2j!(
                                       i::Int,j::Int,dij::Int,stepsize::Int,
                                       facecount::Array{Int,1},numpairs::Array{Int,1},
                                       jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                       r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                       oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
                                       rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
                                       izfull::Array{Int,1},ff2pv::Array{Int,1},
                                       pmhist::Array{Int,2},fpi::Array{Int,2},
                                      npsupp::BitArray{1})

                    process_maxsd_j2end!(
                                         i::Int,j::Int,dij::Int,stepsize::Int,
                                         facecount::Array{Int,1},numpairs::Array{Int,1},
                                         jrv::Array{Int,1},jcp::Array{Int,1},jz::Array{Int,1},
                                         r::Array{Int,1},z::Array{Int,1},pflist::Array{Int,1},
                                         oldclaw::Array{Int,1},zll::Array{Int,1},colsum::Array{Int,1},
                                         rt::Array{Int,1},ct::Array{Int,1},zt::Array{Int,1},
                                         izfull::Array{Int,1},ff2pv::Array{Int,1},
                                         pmhist::Array{Int,2},fpi::Array{Int,2},
                                        npsupp::BitArray{1})
                end
            end
            # update the column pattern and the total number of nonzeros
            # encountered per codim2 face
            c[i+1] = facecount[1]+1
            if sd == maxsd
                colsum[jrv[cran(jcp,i)]].+=1
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
    #gc()
    if dictionaryoutput == true
        D = Dict{String,Any}(
                             "farfaces" => farfaces,
                             "firstv" => firstv,
                             "grain" => grain,
                             "prepairs" => prepairs,
                             "symmat" => symmat,
                             "nvl2ovl"=> vperm)
        return D
    else
        return farfaces,firstv,grain,prepairs,symmat,vperm
    end
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
                process_maxsd_one2i_subroutine!(
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
        kk = jrv[k]	## may have to reindex this
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
                elseif oldclaw[ff2pv[k]]>=farfilt
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
        if izfull[kk]>0
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
    for k = fpi[j,j]:(fpi[j+1,j]-1)
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
    for k = fpi[j+1,j]:(jcp[j+1]-1)
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
    numpairs[1]+=1
    pflist[numpairs[1]] = facecount[1]
    npsupp[k]=false
end

function pairupdatedeluxe!(k::Int,i::Int,j::Int,numpairs::Array{Int,1},facecount::Array{Int,1},pflist::Array{Int,1},ff2pv::Array{Int,1},npsupp::BitArray{1},pmhist::Array{Int,2})
    numpairs[1]+=1
    pmhist[i,j]+=1
    npsupp[k]=false
    pflist[numpairs[1]]=facecount[1]
    ff2pv[k] = i
end

function faceupdate!(facecount::Array{Int,1},r::Array{Int,1},z::Array{Int,1},k::Int,farfilt::Int,stepsize::Int)
    facecount[1]+=1
    if facecount[1]>length(r)
        append!(r,Array{Int}(undef,stepsize))
        append!(z,Array{Int}(undef,stepsize))
    end
    r[facecount].= k
    z[facecount].= farfilt
end

function faceupdatedeluxe!(facecount::Array{Int,1},r::Array{Int,1},z::Array{Int,1},k::Int,farfilt::Int,stepsize::Int,s::Array{Int,1},i::Int)
    facecount[1]+=1
    if facecount[1]>length(r)
        append!(r,Array{Int}(undef,stepsize))
        append!(z,Array{Int}(undef,stepsize))
        append!(s,Array{Int}(undef,stepsize))
    end
    r[facecount].= k
    z[facecount].= farfilt
    s[facecount].= i
end

function saveface(ct::Array{Int,1},kk::Int,colsum::Array{Int,1},farfilt::Int,oldclaw::Array{Int,1},rt::Array{Int,1},zt::Array{Int,1})
    keep = true
    for l = ct[kk]:colsum[kk]
        if  zt[l]>= farfilt && oldclaw[rt[l]]>=farfilt
            keep = false
            break
        end
    end
    return keep
end

function processfpi!(pmhist::Array{Int,2},fpi::Array{Int,2},jcp::Array{Int,1},jrv::Array{Int,1},ff2pv::Array{Int,1},m::Integer)
    for p = 1:m
        for q = jcp[p]:(jcp[p+1]-1)
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

function generate2faces(symmat)
    m = size(symmat,1)
    if issparse(symmat)
        return symmat
    else
        L = 0
        for i = 1:m
            for j = (i+1):m
                if symmat[j,i]>0
                    L+=1
                end
            end
        end
        rowval = Array{Int}(undef,L)
        nzval = Array{Int}(undef,L)
        colptr = Array{Int}(undef,m+1)
        marker = 0
        colptr[1] = 1
        for i = 1:m
            colptr[i+1]=colptr[i]
            for j = (i+1):m
                if symmat[j,i]>0
                    colptr[i+1]+=1
                    rowval[colptr[i+1]-1] = j
                    nzval[colptr[i+1]-1] = symmat[j,i]
                end
            end
        end
    end
    return rowval,colptr,nzval
end

function generate3faces!(farfaces_cell, firstv_cell, grain_cell, prepairs_cell, m, symmat)
    grain::Array{Int,1} = grain_cell[2]
    farfaces::Array{Int,1} = farfaces_cell[2]
    firstv::Array{Int,1} = firstv_cell[2]

    numverts = length(firstv)-1
    numedges = length(farfaces)
    stepsize = size(symmat,1)^2
    facecount= [0]
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
    for i = 2:(m+1)
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
        adist[:].=0
        adist[crows(firstv,farfaces,a)] = crows(firstv,grain,a)
        for ip = cran(firstv,a)
            i = farfaces[ip]
            dai = grain[ip]
            idist[:].=0
            idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
            idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
            for jp = cran(firstv,i)
                if ncheckedges[jp]
                    j = farfaces[jp]
                    dij = grain[jp]
                    if dij <= dai && dij <= adist[j] # note this condition bakes in the req. that j be adjacent to a
                        numpairs+=1
                        ncheckedges[jp] = false
                        clawvec[1:i] .= 0
                        for lp = cran(firstv_hs,j)
                            l = closefaces_higsorted[lp]
                            if l >= i
                                break
                            elseif idist[l]!=0
                                clawvec[l] = min(idist[l],grain_higsorted[lp])
                            end
                        end
                        for kp = cran(firstv,j)
                            k = farfaces[kp]
                            djk = grain[kp]
                            dak = adist[k]
                            dik = idist[k]
                            if dak < dij && dak<djk && dak<dik	# this bakes in req. that dik>0
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
            idist[:].=0
            idist[crows(firstv,farfaces,i)]	= crows(firstv,grain,i)
            idist[crows(firstv_hs,closefaces_higsorted,i)] = crows(firstv_hs,grain_higsorted,i)
            holdi = i
        end
        clawvec[1:i] .= 0
        for lp = cran(firstv_hs,j)
            l = closefaces_higsorted[lp]
            if l >= i
                break
            elseif idist[l]!=0
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
    r[iso]=r
    z[iso]=z
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

##########################################################################################

####	MAIN

##########################################################################################

function inputVmodel2defaultgeneraformat(s, model)
    if typeof(s) == String
        if in(model,["vr","pc"])
            entryformat = "textfile"
        elseif in(model,["complex"])
            entryformat = "sp"
        end
    else
        entryformat = "n/a"
    end
    return entryformat
end

"""
Computes the persistent homology of a filtered complex.
"""
function eirene(s;
        model		= "vr",
        maxdim = 1,
        minrad		= -Inf,
        maxrad		= Inf,
        numrad		= Inf,
        nodrad = [],
        fastop		= true,
        vscale		= "default",
        record		= "cyclerep",
        entryformat = inputVmodel2defaultgeneraformat(s, model),
        pointlabels	= [],
        verbose		= false)

    if model in ["vr","pc"]
        maxsd = maxdim+2
        return persistf2vr(s, maxsd;
                        model = model,
                        minrad = minrad,
                        maxrad = maxrad,
                        numrad = numrad,
                        nodrad = nodrad,
                        fastop = fastop,
                        record = record,
                        entryformat = entryformat,
                        pointlabels = pointlabels,
                        verbose = verbose)
    elseif model == "complex"
        return persistf2complex(s; maxdim=maxdim, entryformat=entryformat, record = record)
    else
        println("Error: the only valid values for keyword <model> are \"vr\", \"pc\", and \"complex\".")
        println("user input:")
        println(model)
    end
end

function genera2autoformat(rv, dp, dv, ev)
    if typeof(rv) == Array{Array{Int}}
        return "segmented complex"
    elseif !isempty(dp)
        return "dp"
    elseif !isempty(dv)
        return "dv"
    elseif !isempty(ev)
        return "ev"
    end
end

function diagonalentries(x)
    if size(x,1) != size(x,2)
        println("error: d should be square")
        return
    end
    m = size(x,2)
    v = zeros(m)
    for p = 1:m
        v[p]	= x[p,p]
    end
    return v
end

function offdiagmin(d::Array{Tv}) where Tv
    if size(d,1) != size(d,2)
        println("error: d should be square")
        return
    end
    v = zeros(Tv,size(d,2))
    for p = 1:size(d,2)
        val1 = empteval(minimum,d[1:p-1,p],Inf)
        val2 = empteval(minimum,d[p+1:end,p],Inf)
        v[p] = min(val1,val2)
    end
    return v
end

function iudsymmat(m)
    x	= rand(m,m)
    for p = 1:m
        for q = 1:p-1
            x[q,p] = x[p,q]
        end
    end
    return x
end

function vrmat(C::Dict)
    if C["input"]["model"] != "vr"
        println("error: <vrmat> only applies to vietoris-rips complexes")
    end
    nvl2ovl = C["nvl2ovl"]
    numpts = length(nvl2ovl)
    ovl2nvl = Array{Int}(undef,numpts)
    ovl2nvl[nvl2ovl] =   1:numpts
    symmat = copy(C["symmat"])
    symmat = symmat[ovl2nvl,ovl2nvl]
    s = Array{Float64}(undef,symmat)
    for p = 1:length(s)
        if symmat[p]	==   0
            s[p] = Inf
        else
            s[p] = C["ocg2rad"][symmat[p]]
        end
    end
    return s
end

function vertexlifemat(d; model="rand", scale=1/2)
    if model == "pc"
        s = Distances.pairwise(Euclidean(), d, dims=2)
    elseif model == "vr"
        s = copy(d)
    elseif model == "rand"
        s = iudsymmat(d)
    end
    v = offdiagmin(s)
    if typeof(scale) <: Number
        for	p = 1:size(s,2)
            s[p,p] = v[p]*scale
        end
    elseif scale == "rand"
        for p = 1:size(s,2)
            r = rand(1)
            s[p,p] = r[1]*v[p]
        end
    else
        println("error: scale must be either a scalar or the string \"rand\"")
        return
    end
    return s
end

function ceil2grid(M; origin=0, stepsize=1, numsteps=Inf)
    if stepsize <  0
        println("error in function <roundentries>: stepsize must be positive")
        return
    end
    if numsteps <  1
        println("error in function <roundentries>: numsteps must be positive")
        return
    end

    N = copy(M)
    N = Array{Float64}(N) # conversion
    N = (N .- origin)./stepsize
    N = ceil.(N)
    N = N.*stepsize.+origin

    N[N.<origin]	.= -Inf
    if numsteps < Inf
        maxval = origin+numsteps*stepsize
        N[N.>maxval].= Inf
    end
    N
end


end;
