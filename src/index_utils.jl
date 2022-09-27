#!/usr/bin/env julia
using SparseArrays

cran(A::SparseMatrixCSC, j::Int) = A.colptr[j]:A.colptr[j+1]-1
cran(colptr::Vector{Int}, j::Int) = colptr[j]:colptr[j+1]-1
function cran(colptr::Vector{Int}, J::Union{UnitRange{Int},Vector{Int}})
    vcat(range.(colptr[J], colptr[J.+1].-1)...)
end
cran(colptr::UnitRange, j) = colptr[j]

crows(A::SparseMatrixCSC, j::Int) = A.rowval[A.colptr[j]:A.colptr[j+1]-1]
crows(colptr::Vector{Int}, rowval::Vector{Int}, j) = rowval[cran(colptr, j)]

function extend!(x::Vector{Int},n::Int)
    if length(x) < n
        append!(x,Vector{Int}(undef,n-length(x)))
    end
end

function copycolumnsubmatrix(Arv::Vector{Int},Acp::Vector{Int},columnindices::Vector{Int})
    fro = Acp[columnindices]
    to = Acp[columnindices .+ 1]
    Brv = Arv[vcat(range.(fro, to .- 1)...)]
    Bcp = Int[1; to .- fro] |> cumsum
    return Brv, Bcp
end

"""
colsinorder must be in sorted order.
"""
function supportedmatrix!(Mrowval::Vector{Int},Mcolptr::Vector{Int},rows1,colsinorder,Mm::Int)
    n = length(colsinorder)
    suppcol1 = falses(Mm)
    suppcol1[rows1] .= true
    cpHolder = 1
    nz1 = 0
    for jp = 1:n
        for ip in cran(Mcolptr, colsinorder[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1 += 1
                Mrowval[nz1] = i
            end
        end
        
        Mcolptr[jp] = cpHolder
        cpHolder = nz1 + 1
    end
    Mcolptr[n+1] = cpHolder
    deleteat!(Mcolptr, n+2:length(Mcolptr))
    deleteat!(Mrowval, Mcolptr[end]:length(Mrowval))
end

function stackedsubmatrices(
        Mrowval::Vector{Int},
        Mcolptr::Vector{Int},
        rows1::Vector{Int},
        rows2::Vector{Int},
        cols::Vector{Int},
        Mm::Int)::Tuple{Vector{Int},Vector{Int},Vector{Int},Vector{Int}}
    
    suppcol1 = falses(Mm)
    suppcol2 = falses(Mm)
    suppcol1[rows1] .= true
    suppcol2[rows2] .= true
    
    Mrowvals = [Mrowval[cran(Mcolptr, cols[jp])] for jp = 1:length(cols)]
    rv1 = [_m[suppcol1[_m]] for _m in Mrowvals]
    rv2 = [_m[suppcol2[_m]] for _m in Mrowvals]
    cp1 = cumsum(Int[1; length.(rv1)])
    cp2 = cumsum(Int[1; length.(rv2)])
    rv1 = vcat(rv1...)
    rv2 = vcat(rv2...)
    return rv1,rv2,cp1,cp2
end

