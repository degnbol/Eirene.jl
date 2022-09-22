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

export eirene, barcode, classrep

include("wasserstein_distances.jl")
include("simplicial_constructions.jl")
include("schur_complements.jl")
include("inversion.jl")
include("multiplication.jl")
include("transposition.jl")
include("barcode_utils.jl")
include("vietorisrips.jl")
include("chain_operations.jl")
include("matrix_utils.jl")
include("perm_utils.jl")
include("index_utils.jl")

"""
Computes the persistent homology of a filtered complex.
"""
function eirene(d::Matrix{Float64}, maxdim; minrad=-Inf, maxrad=Inf, numrad=Inf, nodrad=[])
    numpoints = size(d, 1)
    @assert issymmetric(d)
    
    maxrad = min(maxrad, minimum(maximum(d,dims=1)))
    
    d = minmaxceil(d, minrad=minrad, maxrad=maxrad, numrad=numrad)

    # recall that t will have the same order (NOT inverted) as d
    # <trueordercanonicalform> is a bit like <integersinsameorder>, 
    # just valid for floating point inputs, and with a bit more data in the 
    # output
    t, ocg2rad = trueordercanonicalform(d, factor=true)

    t = (1 + maximum(t)) .- t
    ocg2rad = reverse(ocg2rad, dims=1)

    if any(d .> maxrad)
        t .-= 1
        deleteat!(ocg2rad, 1)
    end

    # this step is necessary in order to cover the case where some vertices 
    # never enter the filtration
    vertices2keep = findall(diag(t) .!= 0)
    t = t[vertices2keep, vertices2keep]

    #### Build the complex
    D = buildcomplex3(t, maxdim+2)
    D["ocg2rad"] = ocg2rad

    #### Compute persistence
    persistf2!(D)

    # this covers the case where some vertices never enter the filtration
    D["nvl2ovl"] = vertices2keep[D["nvl2ovl"]]

    #### Store generators
    #gc()
    unpack!(D)
    #gc()
    return D
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

    N[N .< origin] .= -Inf
    if numsteps < Inf
        maxval = origin + numsteps * stepsize
        N[N .> maxval] .= Inf
    end
    N
end


function classrep(D::Dict; dim=1, class=1, format="vertex x simplex")
    if any(class .> nnzbars(D, dim=dim))
        error("The value for keyword argument <class> has an integer greater than the number of nonzero bars in the specified dimension")
    elseif !(0 <= dim <= D["input"]["maxdim"])
        error("Barcodes were not computed in the specified dimension")
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
    bc[mortalran,1] .= mortalprimagrain
    bc[mortalran,2] .= mortalultragrain
    bc[evergrran,1] = lg[tid[evrgrbran]]
    
    if !ocf
        bcc = copy(bc)
        bc = Array{Float64}(bc)
        bc[finran] = D["ocg2rad"][bcc[finran]]
        bc[evergrran,2]    .= Inf
    else
        bc = length(D["ocg2rad"]).-bc
    end
    
    bc
end


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


end;
