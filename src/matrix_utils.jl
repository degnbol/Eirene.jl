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


function offdiagmin(S, i)
    if i == 1
        return(minimum(S[2:end,i]))
    elseif i == size(S,1)
        return minimum(S[1:end-1,i])
    else
        return min(minimum(S[1:i-1,i]),minimum(S[i+1:end,i]))
    end
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


