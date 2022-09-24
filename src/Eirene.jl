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
using LinearAlgebra
using Statistics: mean

export eirene, barcode, classrep

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
function eirene(pointcloud::Matrix{Float64}, maxdim::Int; minrad=-Inf, maxrad=Inf, numrad=Inf, nodrad=[])
    @assert size(pointcloud,2) == 3 "points in R3 expected, with points along the first dimension."
    numpoints = size(pointcloud,1)
    d = pairwise(Euclidean(), pointcloud, dims=1)
    
    maxrad = min(maxrad, minimum(maximum(d,dims=1)))
    
    d = minmaxceil(d, minrad=minrad, maxrad=maxrad, numrad=numrad)

    # recall that t will have the same order (NOT inverted) as d
    # <trueordercanonicalform> is a bit like <integersinsameorder>, 
    # just valid for floating point inputs, and with a bit more data in the 
    # output
    t, ocg2rad = trueordercanonicalform(d)
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
    C = buildcomplex3(t, maxdim+2)
    # this covers the case where some vertices never enter the filtration
    nvl2ovl = vertices2keep[C.nvl2ovl]
    
    #### Compute persistence
    P = persistf2(C.farfaces, C.firstv, C.prepairs, C.grain)
    
    #### Store generators
    R = unpack!(C.grain, C.farfaces, C.firstv, P.trv, P.tcp, P.plo, P.phi, P.tid, maxdim+2)
    
    barcodes = [barcode(maxdim, C.grain, P.plo, P.phi, P.tid, ocg2rad; dim=dim) for dim in 1:maxdim]
    representatives = [[classrep(maxdim, C.farfaces, C.firstv, R.cyclerep, nvl2ovl, C.grain, P.plo, P.phi, P.tid; class=i, dim=dim) for i in 1:size(barcodes[dim],1)] for dim in 1:maxdim]
    return barcodes, representatives
end

function barcode(maxdim::Int, grain, plo, phi, tid, ocg2rad; dim::Int=1)
    @assert dim <= maxdim
    sd = dim+2
    plo = plo[sd]
    phi = phi[sd]
    tid = tid[sd]
    lg = grain[sd-1]
    hg = grain[sd]
    
    mortalprimagrain = lg[plo]
    mortalultragrain = hg[phi]
    
    finind = findall(mortalprimagrain .!= mortalultragrain)
    numfin = length(finind)
    numinf = length(tid)-length(plo)
    
    mortalprimagrain = mortalprimagrain[finind]
    mortalultragrain = mortalultragrain[finind]
    
    mortalran = 1:numfin
    evergrran = numfin+1:numfin+numinf
    finran = 1:2*numfin+numinf
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
    
    bcc = copy(bc)
    bc = Array{Float64}(bc)
    bc[finran] = ocg2rad[bcc[finran]]
    bc[evergrran,2] .= Inf
    
    bc
end

function classrep(maxdim::Int, farfaces, firstv, cyclerep, nvl2ovl, grain, plo, phi, tid; dim::Int=1, class::Int=1)
    @assert !any(class .> nnzbars(grain, plo, phi, tid, dim=dim)) "The value for keyword argument <class> has an integer greater than the number of nonzero bars in the specified dimension"
    @assert 0 <= dim <= maxdim "Barcodes were not computed in the specified dimension"
    classrep_faces(farfaces, firstv, cyclerep, nvl2ovl; dim=dim, class=class)
end

function classrep_faces(farfaces, firstv, cyclerep, nvl2ovl; dim::Int=1, class::Int=1)
    sd = dim+2
    vrealization = vertexrealization(farfaces, firstv, sd-1, cyclerep[sd][class])
    nvl2ovl[vrealization]
end


# run an example to precompile typed versions of each function called
eirene(rand(100, 3), 2; minrad=0)

end;


using CSV, DataFrames

fnames = readdir(expanduser("~/protTDA/data/GASS/xyzChain"); join=true)
dfs = CSV.read.(fnames, DataFrame)

xyzs = [Matrix(df[!, [:x, :y, :z]]) for df in dfs]

@time Eirene.eirene.(xyzs[1:8], 2; minrad=0.);

@time eirene.(pairwise.(Ref(Euclidean()), xyzs[1:8], dims=1); maxdim=2, minrad=0.);



@time for (fname, xyz) in zip(basename.(fnames), xyzs[1:4])
    bs, rs = Eirene.eirene(xyz, 2; minrad=0.)
    
    open("lars/$fname", "w") do io
        obj = Dict("H$i" => (barcode=b, representatives=r) for (i,(b,r)) ∈ enumerate(zip(bs, rs)))
        JSON.print(io, obj)
    end
end

