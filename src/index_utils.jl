#!/usr/bin/env julia
using SparseArrays

cran(A::SparseMatrixCSC, j) = A.colptr[j]:A.colptr[j+1]-1
cran(colptr::Array, j::Int) = colptr[j]:colptr[j+1]-1
function cran(colptr::Array, J::Union{UnitRange{Int},Vector{Int}})
    vcat(range.(colptr[J], colptr[J.+1].-1)...)
end
cran(colptr::UnitRange, j) = colptr[j]

crows(A::SparseMatrixCSC, j) = A.rowval[A.colptr[j]:A.colptr[j+1]-1]
crows(colptr::Array, rowval::Array, j) = rowval[cran(colptr, j)]

function extend!(x::Array{Int,1},n::Int)
    if length(x) < n
        append!(x,Array{Int}(undef,n-length(x)))
    end
end

function copycolumnsubmatrix(Arv::Array{Int,1},Acp,columnindices)
    allocationspace = 0
    for j in columnindices
        allocationspace+= Acp[j+1]-Acp[j]
    end
    Brv = Array{Int}(undef,allocationspace)
    Bcp = Array{Int}(undef,length(columnindices)+1)
    Bcp[1] = 1
    for jp = 1:length(columnindices)
        j = columnindices[jp]
        Bcp[jp+1] = Bcp[jp]+Acp[j+1]-Acp[j]
        Brv[Bcp[jp]:(Bcp[jp+1]-1)] = Arv[Acp[j]:(Acp[j+1]-1)]
    end
    return Brv, Bcp
end


function extendcolumnlight!(rowval::Array{Ti,1},colptr::Array{Ti,1},v::Array{Ti},k::Ti,growthincrement::Ti) where Ti
    r = rowval
    c = colptr
    startpoint = copy(c[k+1])
    c[k+1]=c[k]+length(v)
    if length(r)<c[k+1]-1
        append!(r,Array{Int}(undef,max(growthincrement,length(v))))
    end
    r[startpoint:(c[k+1]-1)] = v
end

# colsinorder must be in sorted order
function supportedmatrix!(Mrowval::Array{Int},Mcolptr::Array{Int,1},rows1,colsinorder,Mm::Int)
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
    rv1 = Array{Int}(undef,nz1)
    rv2 = Array{Int}(undef,nz2)
    cp1 = Array{Int}(undef,n+1); cp1[1]=1
    cp2 = Array{Int}(undef,n+1); cp2[1]=1
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

