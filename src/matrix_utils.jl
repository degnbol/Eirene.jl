#!/usr/bin/env julia

"""
The output of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
"""
function minmaxceil(N::Matrix{Float64}; minrad=minimum(N),maxrad=maximum(N), numrad=Inf)
    # NB: this is probably an important step; several pernicious problems in 
    # the development of rounding procedures turned out to be rooted in 
    # accidental rewriting of matrix entries
    S = copy(N)
    return minmaxceil!(S; minrad=minrad, maxrad=maxrad, numrad=numrad)
end

"""
NB THIS FUNCTION MODIFIES ITS INPUT

The result of this function is obtatined by *first* rounding all values below
minrad up to minrad, *then* rounding all values above maxrad up to Inf, *then*
rounding the entries valued in [minrad,maxrad] to the nearest greater element of
[minrad:numrad:maxrad], where by definition [minrad:Inf:maxrad] is the closed
inteveral containing all reals between minrad and maxrad.
If minrad == -Inf, then we set it to minimum(N) before performing any operations.
If maxrad == Inf, then we set it to maximum(N) before performing any operations.
"""
function minmaxceil!(S::Matrix{Float64}; minrad=minimum(N), maxrad=maximum(N), numrad=Inf)
    if (minrad == Inf) || (maxrad == -Inf)
        return fill(Inf, size(S)...)
    end

    S[S .< minrad] .= minrad
    S[S .> maxrad] .= Inf

    # stands for finite indices
    fi = (LinearIndices(S))[findall(isfinite, S)]
    fv = S[fi]

    isempty(fv) && return S

    if minrad == -Inf
        minrad = minimum(fv)
    end
    if maxrad == Inf
        maxrad = maximum(fv)
    end
    if numrad == 1
        S[fi] .= maxrad
        return S
    end
    numrad == Inf && return S

    ran = range(minrad, maxrad, length=numrad) .|> Float64
    ran[end] = maxrad

    S[fi] = ceil2grid(fv, ran)
    S
end

# ran should be an array in sorted order, with ran[end] >= maximum(A)
function ceil2grid(A, ran)
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


function offdiagmean(S::Matrix{Int}; default=0)
    m,n = size(S)
    @assert m == n "offdiagmean: matrix arg should be square"
    m != 1 || return default
    mu = zeros(m)
    for j = 1:m
        v = S[1:(j-1),j]
        u = S[(j+1):m,j]
        mu[j] = mean(vcat(v[:],u[:]))
    end
    return mu
end

function trueordercanonicalform(M::Matrix{Float64})
    m = length(M)
    m > 0 || return zeros(Int,0), zeros(Int,0)
    perm = sortperm(M[:], alg=MergeSort)
    oca = zeros(Int, size(M)...)

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

    trueordercanonicalform_subr!(M,perm,oca,oca2rad,[1],[1],m)
    return oca, oca2rad
end

function trueordercanonicalform_subr!(M, perm, oca, oca2rad, post, k, m::Int)
    for p  = 1:m
        if M[perm[p]] != M[perm[post[1]]]
            post[1] = p
            k[1] = k[1]+1
            oca2rad[k[1]] = M[perm[post[1]]]
        end
        oca[perm[p]] = k[1]
    end
end


