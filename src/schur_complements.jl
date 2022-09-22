#!/usr/bin/env julia

function copycolind2colind!(rowvalA::Array{Tv,1},colptrA::Array{Tv,1},columnindices,rowvalB::Array{Tv,1},colptrB::Array{Tv,1},startingDestination::Integer,growthIncrement::Integer)  where Tv<:Integer
    numnewrows = sum(colptrA[columnindices .+ 1] .- colptrA[columnindices])

    if length(rowvalB)<colptrB[startingDestination]+numnewrows-1
        append!(rowvalB,Array{Tv}(undef,numnewrows+growthIncrement))
    end
    if length(colptrB)<startingDestination+length(columnindices)
        append!(colptrB,Array{Tv}(undef,startingDestination+length(columnindices)))
    end
    colptrB[1] = 1
    for i = 1:length(columnindices)
        k = startingDestination+i #the index of colptr pointing to end of this new column
        col = columnindices[i]
        colptrB[k]=colptrB[k-1]+colptrA[col+1]-colptrA[col]
        rowvalB[colptrB[k-1]:(colptrB[k]-1)]=rowvalA[colptrA[col]:(colptrA[col+1]-1)]
    end
end

function finddownstreamelements_embeddedupperunitriangularmatrix(
        Mrv, Mcp, Mm, initialelements::Array{Int,1}, prows::Array{Int,1}, pcols::Array{Int,1})
    if length(prows) != length(pcols)
        print("length of p doesn't match length of q")
        return
    elseif length(prows) == 0
        return Array{Int}(undef,0)
    end
    n = length(prows)
    rowtranslator = Array{Int}(undef,Mm)
    rowtranslator[prows] .= 1:n
    prowsupp = falses(Mm)
    prowsupp[prows] .= true
    downstreamsupport = falses(n)
    for i = 1:length(initialelements)
        row = initialelements[i]
        if prowsupp[row]
            downstreamsupport[rowtranslator[row]] = true
        end
    end
    for jp = n:-1:1
        j = pcols[jp]
        ran = cran(Mcp,j)
        for ip in ran
            rawrow = Mrv[ip]
            if prowsupp[rawrow] && downstreamsupport[rowtranslator[rawrow]]
                for kp in ran
                    rawrow = Mrv[kp]
                    if prowsupp[rawrow]
                        downstreamsupport[rowtranslator[rawrow]] = true
                    end
                end
                break
            end
        end
    end
    counter = sum(downstreamsupport)
    downstreamelements = Array{Int}(undef,counter)
    counter = 0
    for i = 1:n
        if downstreamsupport[i]
            counter+=1
            downstreamelements[counter] = i
        end
    end
    downstreamelements
end

function schurit4!(	Mrv,Mcp,Mm,Mn,Mn0,
        rowlab,collab,
        Jprows,Jpcols,numjunpairs,
        Sprows,Spcols,numsenpairs,
        comprows,compcols,
        Trv,Tcp,Srv,Scp;
        updatetransform = true,
        diagonstic = false)

    Mm0 = copy(Mm[1])
    Mm[1] = length(comprows)
    Mn[1] = length(compcols)

    copycolind2colind!(Trv,Tcp,Jpcols[numjunpairs[1]:-1:1],Srv,Scp,numsenpairs[1]+1,0)

    topspot = numsenpairs[1]+numjunpairs[1]
    for i = 1:numjunpairs[1]
        placementlocation = numsenpairs[1]+i
        extractionlocation = numjunpairs[1]-i+1
        Sprows[placementlocation] = rowlab[Jprows[extractionlocation]]
        Spcols[placementlocation] = collab[Jpcols[extractionlocation]]
    end

    numsenpairs[1] += numjunpairs[1]

    rowsum = falses(Mm0)
    for jp = 1:Mn[1]
        j = compcols[jp]
        for ip = cran(Mcp, j)
            rowsum[Mrv[ip]] = true
        end
    end

    keptlist = finddownstreamelements_embeddedupperunitriangularmatrix(
        Mrv,Mcp,Mm0,findall(rowsum),Jprows[1:numjunpairs[1]],Jpcols[1:numjunpairs[1]])

    keptmarker = length(keptlist)
    prows = Array{Int}(undef,keptmarker)
    pcols = Array{Int}(undef,keptmarker)
    for i = 1:keptmarker
        keptindex = keptlist[i]
        prows[i] = Jprows[keptindex]
        pcols[i] = Jpcols[keptindex]
    end

    Arv,Crv,Acp,Ccp = stackedsubmatrices(Mrv,Mcp,prows,comprows,pcols,Mm0)
    Brv,Drv,Bcp,Dcp = stackedsubmatrices(Mrv,Mcp,prows,comprows,compcols,Mm0)
    Lrv,Lcp = copycolumnsubmatrix(Trv,Tcp,pcols)
    Rrv,Rcp = copycolumnsubmatrix(Trv,Tcp,compcols)

    translator = Array{Int}(undef,Mm0)
    translator[prows] = 1:keptmarker
    Arv .= translator[Arv]
    Brv .= translator[Brv] 
    translator[comprows] = 1:length(comprows)
    Crv .= translator[Crv]
    Drv .= translator[Drv]

    Airv,Aicp = morseInverseF2orderedColsUnsortedRowsInSilentOut(Arv,Acp)
    Brv,Bcp = spmmF2silentLeft(Airv,Aicp,Brv,Bcp,keptmarker)

    for j = 1:keptmarker
        # repurposing rowsum as col-to-row translator
        translator[j]=collab[pcols[j]]
    end
    collabcopy = copy(collab)
    for i = 1:Mm[1]
        rowlab[i] = rowlab[comprows[i]]
    end
    for j = 1:Mn[1]
        collab[j] = collab[compcols[j]]
    end

    blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Mm[1],Mn[1],Mrv,Mcp,0)

    if updatetransform
        blockprodsumsilenticolsleftWrite2!(Lrv,Lcp,Brv,Bcp,Rrv,Rcp,Mn0,Mn[1],Trv,Tcp,0,translator)
    end
end

