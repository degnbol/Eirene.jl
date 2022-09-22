#!/usr/bin/env julia

function transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Am) where Ti
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    # Attach destination matrix
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(undef,Am+1)
    Crowval = Array{Ti}(undef,Annz)
    # Compute the column counts of C and store them shifted forward by one in
    # Ccolptr
    Ccolptr[1:end] .= 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and
    # Cnzval, tracking write positions in Ccolptr
    for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr
    # shift, but the first colptr remains incorrect
    Ccolptr[1] = 1

    return Crowval, Ccolptr
end

function transposeLighter(Arowval::Array{Ti},Acolptr::Array{Ti},Anzval::Array{Tv},Am::Integer) where {Tv, Ti}
    Annz = Acolptr[end]-1
    An = length(Acolptr)-1
    Cm = An
    Cn = Am
    Ccolptr = Array{Ti}(undef,Am+1)
    Crowval = Array{Ti}(undef,Annz)
    Cnzval = Array{Tv}(undef,Annz)
    # Compute the column counts of C and store them shifted forward by one in Ccolptr
    Ccolptr[1:end] .= 0
    @inbounds for k in 1:Annz
        Ccolptr[Arowval[k]+1] += 1
    end
    # From these column counts, compute C's column pointers
    # and store them shifted forward by one in Ccolptr
    countsum = 1
    @inbounds for k in 2:(Cn+1)
        overwritten = Ccolptr[k]
        Ccolptr[k] = countsum
        countsum += overwritten
    end
    # Distribution-sort the row indices and nonzero values into Crowval and Cnzval,
    # tracking write positions in Ccolptr
    @inbounds for Aj in 1:An
        for Ak in Acolptr[Aj]:(Acolptr[Aj+1]-1)
            Ai = Arowval[Ak]
            Ck = Ccolptr[Ai+1]
            Crowval[Ck] = Aj
            Cnzval[Ck] = Anzval[Ak]
            Ccolptr[Ai+1] += 1
        end
    end
    # Tracking write positions in Ccolptr as in the last block fixes the colptr shift,
    # but the first colptr remains incorrect
    Ccolptr[1] = 1
    return Crowval, Ccolptr, Cnzval
end

#=
Accepts an mxn integer array, M.  To M we implicitly associate an array N,
as follows. If M has no nonzero entries, then N is the zeron-one array such that
supp(N[:,i]) = M[:,i].  Here we regard M[:,i] as a set.  I *believe* one ignores
duplicate entries, but have not checked. A similar interpretation holds when M
has zero entries - one simply discards the zeros.
Let S denote the support of N, r denote row02row1translator, and c denote
col02col1translator.  This function returns the data specifying a sparse
zero-one matrix whose support is
{ [c[j],r[i]] : [i,j] \in N and neither c[j] nor r[i] are zero }.
=#
function presparsefull2unsortedsparsetranspose(
        M::Array{Tv,2},
        row02row1translator,
        col02col1translator;
        verbose::Bool=false) where Tv<:Integer

    Mm,Mn = size(M)

    if Mn == 0
        rowval1 = Array{Int}(undef,0)
        if isempty(row02row1translator)
            colptr1 = ones(Int,1)
        else
            colptr1 = ones(Int,1+maximum(row02row1translator))
        end
        return rowval1,colptr1,Mn
    end

    for i = 1:(Mm*Mn)
        if M[i]>0
            M[i] = row02row1translator[M[i]]
        end
    end
    m0 = maximum(row02row1translator)
    rowcounter = zeros(Int,m0)
    for k in M  #allow for zero values
        if k>0
            rowcounter[k]+=1
        end
    end
    colptr1 = Array{Int}(undef,m0+1)
    colptr1[1]=1
    for i = 1:m0
        colptr1[i+1]=colptr1[i]+rowcounter[i]
    end
    rowval1 = Array{Int}(undef,colptr1[end]-1)
    placer = copy(colptr1)
    if verbose
        coverageTestVec = trues(colptr1[end]-1)
    end
    for j = 1:Mn
        jval = col02col1translator[j]
        for i = 1:Mm
            row = M[i,j]
            if row > 0
                if verbose
                    coverageTestVec[placer[row]]=false
                end
                rowval1[placer[row]]=jval
                placer[row]+=1
            end
        end
    end
    if verbose
        if any(coverageTestVec)
            print("please refer to the coverageTestVec")
            sleep(4)
        end
        println([length(rowval1) "length(output of presparsefull2unsortedsparsetranspose)"])
    end
    #gc()
    return rowval1,colptr1,Mn
end

