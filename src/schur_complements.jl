#!/usr/bin/env julia

function copycolind2colind!(rowvalA::Vector{Int}, colptrA::Vector{Int}, columnindices, rowvalB::Vector{Int},colptrB::Vector{Int}, startingDestination::Int, growthIncrement::Int)
    numnewrows = sum(colptrA[columnindices .+ 1] .- colptrA[columnindices])

    if length(rowvalB) < colptrB[startingDestination] + numnewrows - 1
        append!(rowvalB, Vector{Int}(undef,numnewrows+growthIncrement))
    end
    if length(colptrB) < startingDestination + length(columnindices)
        append!(colptrB, Vector{Int}(undef,startingDestination + length(columnindices)))
    end
    colptrB[1] = 1
    for i = 1:length(columnindices)
        # the index of colptr pointing to end of this new column
        k = startingDestination + i
        col = columnindices[i]
        colptrB[k] = colptrB[k-1] + colptrA[col+1] - colptrA[col]
        rowvalB[colptrB[k-1]:colptrB[k]-1] = view(rowvalA, colptrA[col]:colptrA[col+1]-1)
    end
end

function finddownstreamelements_embeddedupperunitriangularmatrix(
        Mrv, Mcp, Mm::Int, initialelements, prows::Vector{Int}, pcols::Vector{Int})
    @assert length(prows) == length(pcols) "length of p doesn't match length of q"
    length(prows) == 0 && return Int[]
    n = length(prows)
    rowtranslator = Vector{Int}(undef, Mm)
    rowtranslator[prows] .= 1:n
    prowsupp = falses(Mm)
    prowsupp[prows] .= true
    downstreamsupport = falses(n)
    downstreamsupport[view(rowtranslator, initialelements .& prowsupp)] .= true
    for jp = n:-1:1
        ran = cran(Mcp, pcols[jp])
        for ip in ran
            rawrow = Mrv[ip]
            if prowsupp[rawrow] && downstreamsupport[rowtranslator[rawrow]]
                rawrows = view(Mrv, ran)
                downstreamsupport[view(rowtranslator, view(rawrows, view(prowsupp, rawrows)))] .= true
                break
            end
        end
    end
    findall(downstreamsupport)
end

function schurit4!(Mrv::Vector{Int},Mcp::Vector{Int},Mm::Int,Mn::Int,Mn0::Int,
        rowlab,collab,Jprows,Jpcols,numjunpairs::Int,Sprows,Spcols,numsenpairs::Int,
        comprows,compcols,Trv,Tcp,Srv,Scp)
    println("schurit4!")
    @time begin
    Mm0 = Mm
    Mm = length(comprows)
    Mn = length(compcols)

    copycolind2colind!(Trv,Tcp,Jpcols[numjunpairs:-1:1],Srv,Scp,numsenpairs+1,0)

    Sprows[numsenpairs .+ (1:numjunpairs)] .= rowlab[view(Jprows, numjunpairs:-1:1)]
    Spcols[numsenpairs .+ (1:numjunpairs)] .= collab[view(Jpcols, numjunpairs:-1:1)]

    numsenpairs += numjunpairs

    rowsum = falses(Mm0)
    for j = compcols[1:Mn]
        rowsum[view(Mrv, cran(Mcp, j))] .= true
    end
    
    keptlist = finddownstreamelements_embeddedupperunitriangularmatrix(
        Mrv, Mcp, Mm0, rowsum, Jprows[1:numjunpairs], Jpcols[1:numjunpairs])
    keptmarker = length(keptlist)
    
    prows = view(Jprows, keptlist)
    pcols = view(Jpcols, keptlist)
    
    Arv,Crv,Acp,Ccp = stackedsubmatrices(Mrv,Mcp,prows,comprows,pcols,Mm0)
    Brv,Drv,Bcp,Dcp = stackedsubmatrices(Mrv,Mcp,prows,comprows,compcols,Mm0)
    Lrv,Lcp = copycolumnsubmatrix(Trv,Tcp,pcols)
    Rrv,Rcp = copycolumnsubmatrix(Trv,Tcp,compcols)
    translator = Vector{Int}(undef,Mm0)
    translator[prows] = 1:keptmarker
    Arv .= translator[Arv]
    Brv .= translator[Brv] 
    translator[comprows] = 1:length(comprows)
    Crv .= translator[Crv]
    Drv .= translator[Drv]

    Airv, Aicp = morseInverseF2orderedColsUnsortedRowsInSilentOut(Arv,Acp)
    Brv, Bcp = spmmF2silentLeft(Airv,Aicp,Brv,Bcp,keptmarker)

    # repurposing rowsum as col-to-row translator
    translator[1:keptmarker] .= view(collab, pcols)
    rowlab[1:Mm] = rowlab[view(comprows, 1:Mm)]
    collab[1:Mn] = collab[view(compcols, 1:Mn)]
    
    extend!(Mcp, Mn+1)
    deleteat!(Mcp, Mn+2:length(Mcp))
    Mcp[1] = 1

    rowSupp = zeros(Int, Mm)
    rowList = Array{Int}(undef, Mm)
    rowvalCj = BitArray(undef, Mm)
    for i in 1:Mn
        if length(Mrv) < Mcp[i] + Mm
            append!(Mrv, Vector{Int}(undef, Mm))
        end
        eyerange = cran(Dcp, i)
        newrowscounter = 0
        for kp = eyerange
            row = Drv[kp]
            rowSupp[row] = i
            newrowscounter += 1
            rowList[newrowscounter] = row
            rowvalCj[row] = true
        end
        if Bcp[i] < Bcp[Mn+1]
            for jp in Bcp[i]:(Bcp[i+1] - 1)
                j = Brv[jp]
                for kp in Ccp[j]:(Ccp[j+1] - 1)
                    row = Crv[kp]
                    if rowSupp[row] != i
                        rowSupp[row] = i
                        newrowscounter += 1
                        rowList[newrowscounter] = row
                        rowvalCj[row] = true
                    else
                        rowvalCj[row] = !rowvalCj[row]
                    end
                end
            end
        end
        ii = i+1
        Mcp[ii] = Mcp[i]
        for k = 1:newrowscounter
            row = rowList[k]
            if rowvalCj[row]
                Mrv[Mcp[ii]] = row
                Mcp[ii] += 1
            end
        end
    end
    blockprodsumsilenticolsleftWrite2!(Lrv,Lcp,Brv,Bcp,Rrv,Rcp,Mn0,Mn,Trv,Tcp,0,translator)
    end
    Mm, Mn, numsenpairs
end

