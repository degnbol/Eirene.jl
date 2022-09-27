#!/usr/bin/env julia

"""
Returns the permutation z on {1,...,length(v)} such z[i]<z[j] iff either
(a) v[i] < v[j], or
(b) v[i] = v[j] and i < j
"""
function integersinsameorder(v::Array{Int,1})
    isempty(v) && return Array{Int}(undef,0)
    m = length(v)
    maxv = maximum(v)
    minv = minimum(v)
    minv = minv-1;
    x = zeros(Int, maxv - minv)
    z = Array{Int}(undef, length(v))
    for i = 1:m
        x[v[i]-minv] += 1
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
    z
end

"""
- In beta; should be compared with integersinsameorderbycolumn3.  See
/Users/greghenselman/Google Drive/GregDirectory/julia_gd/Julia/workshop/workshop_Oct2017.jl
-   Functionally equivalent to integersinsameorderbycolumn; returns a
permutation z on {1,...,length(v)} so that for all j
- cran(colptr,j) maps to cran(colptr,j), and
- crows(colptr,v[z],j) is an array in sorted order
"""
function integersinsameorderbycolumn2(v::Array{Int,1}, colptr)
    numcols = length(colptr)-1
    m = length(v)
    v .-= minimum(v) - 1
    x = zeros(Int, maximum(v))
    z = Array{Int}(undef,length(v))
    for j = 1:numcols
        colptr[j] != colptr[j+1] || continue
        vs = v[colptr[j]:colptr[j+1]-1]
        for _v in vs x[_v] += 1 end
        maxv = maximum(vs)
        minv = minimum(vs)
        prevsum = colptr[j]
        for i = minv:maxv
            sum = prevsum + x[i]
            x[i] = prevsum
            prevsum = sum
        end
        for i = colptr[j]:(colptr[j+1]-1)
            u = v[i]
            z[i] = x[u]
            x[u] += 1
        end
        for i = minv:maxv
            x[i] = 0
        end
    end
    z
end

