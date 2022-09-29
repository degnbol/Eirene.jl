#!/usr/bin/env julia
using SparseArrays

cran(A::SparseMatrixCSC, j::Int) = A.colptr[j]:A.colptr[j+1]-1
cran(colptr::Vector{Int}, j::Int) = colptr[j]:colptr[j+1]-1
function cran(colptr::Vector{Int}, J::Union{UnitRange{Int},Vector{Int}})
    vcat(range.(view(colptr, J), view(colptr, J .+ 1) .- 1)...)
end
cran(colptr::UnitRange, j) = colptr[j]

crows(A::SparseMatrixCSC, j::Int) = A.rowval[A.colptr[j]:A.colptr[j+1]-1]
crows(colptr::Vector{Int}, rowval::Vector{Int}, j) = view(rowval, cran(colptr, j))

function extend!(x::Vector{Int},n::Int)
    if length(x) < n
        append!(x, Vector{Int}(undef, n-length(x)))
    end
end

function copycolumnsubmatrix(Arv::Vector{Int},Acp::Vector{Int},columnindices)
    fro = view(Acp, columnindices)
    to = view(Acp, columnindices .+ 1)
    Brv = Arv[vcat(range.(fro, to .- 1)...)]
    Bcp = Int[1; to .- fro] |> cumsum
    return Brv, Bcp
end

"""
colsinorder must be in sorted order.
"""
function supportedmatrix!(Mrowval::Vector{Int}, Mcolptr::Vector{Int}, rows1, colsinorder,Mm::Int)
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

function stackedsubmatrices(Mrowval::Vector{Int}, Mcolptr::Vector{Int}, rows1, rows2, cols, Mm::Int)
    n = length(cols)
    suppcol1 = falses(Mm)
    suppcol2 = falses(Mm)
    suppcol1[rows1] .= true
    suppcol2[rows2] .= true
    nz1 = 0
    nz2 = 0
    for jp = 1:n
        for ip in cran(Mcolptr,cols[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1 += 1
            elseif suppcol2[i]
                nz2 += 1
            end
        end
    end
    rv1 = Array{Int}(undef,nz1)
    rv2 = Array{Int}(undef,nz2)
    cp1 = Array{Int}(undef,n+1)
    cp2 = Array{Int}(undef,n+1)
    cp1[1] = 1
    cp2[1] = 1
    nz1 = 0
    nz2 = 0
    for jp = 1:n
        for ip in cran(Mcolptr, cols[jp])
            i = Mrowval[ip]
            if suppcol1[i]
                nz1 += 1
                rv1[nz1] = i
            elseif suppcol2[i]
                nz2 += 1
                rv2[nz2] = i
            end
        end
        cp1[jp+1] = nz1 + 1
        cp2[jp+1] = nz2 + 1
    end
    rv1, rv2, cp1, cp2
end

