#!/usr/bin/env julia

function spmmF2silentLeft(Arowval::Vector{Int},Acolptr::Vector{Int},Browval::Vector{Int},Bcolptr::Vector{Int},Am)
    mA = Am
    nB = length(Bcolptr)-1
    rowvalA = Arowval; colptrA = Acolptr
    rowvalB = Browval; colptrB = Bcolptr
    preallocationIncrement = colptrA[end]+colptrB[end]

    colptrC = Vector{Int}(undef,nB+1)
    colptrC[1] = 1
    rowSupp = zeros(Int,mA)
    rowList = Vector{Int}(undef,mA)
    rowvalCj = Vector{Bool}(undef,mA)
    rowvalC = Vector{Int}(undef,preallocationIncrement)
    for i in 1:nB
        newrowscounter = 0
        eyerange = cran(Bcolptr,i)
        newrowscounter = length(eyerange)
        for l=1:newrowscounter
            ll = rowvalB[eyerange[l]]
            rowList[l] = ll
            rowSupp[ll] = i
            rowvalCj[ll] = true
        end
        for jp in colptrB[i]:(colptrB[i+1] - 1)
            j = rowvalB[jp]
            for kp in colptrA[j]:(colptrA[j+1] - 1)
                k = rowvalA[kp]
                if rowSupp[k] != i
                    rowSupp[k] = i
                    newrowscounter +=1
                    rowList[newrowscounter] = k
                    rowvalCj[k] = true
                else
                    rowvalCj[k] = !rowvalCj[k]
                end
            end
        end
        nzRows = findall(rowvalCj[rowList[1:newrowscounter]])
        colptrC[i+1] = colptrC[i]+length(nzRows)

        if colptrC[i+1] > length(rowvalC) + 1
            append!(rowvalC, Vector{Int}(undef, preallocationIncrement))
        end
        rowvalC[colptrC[i]:colptrC[i+1]-1] = sort(rowList[nzRows])
    end
    deleteat!(rowvalC,colptrC[end]:length(rowvalC))
    rowvalC, colptrC
end


function blockprodsumsilenticolsleftWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement,col2silenti)
    extend!(Mcp, Dn + 1)
    deleteat!(Mcp,Dn+2:length(Mcp))
    Mcp[1] = 1

    rowSupp = zeros(Int, Dm)
    rowList = Vector{Int}(undef, Dm)
    rowvalCj = BitArray(undef, Dm)
    for i in 1:Dn
        if length(Mrv) < Mcp[i] + Dm
            extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
        end
        eyerange = cran(Dcp,i)
        newrowscounter = 0
        for kp = eyerange
            row = Drv[kp]
            rowSupp[row] = i
            newrowscounter += 1
            rowList[newrowscounter] = row
            rowvalCj[row] = true
        end
        for jp in Bcp[i]:(Bcp[i+1] - 1)
            j = Brv[jp]
            row = col2silenti[j]
            if rowSupp[row] != i
                rowSupp[row] = i
                newrowscounter +=1
                rowList[newrowscounter] = row
                rowvalCj[row] = true
            else
                rowvalCj[row] = !rowvalCj[row]
            end
        end
        for jp in Bcp[i]:(Bcp[i+1] - 1)
            j = Brv[jp]
            for kp in Ccp[j]:(Ccp[j+1] - 1)
                row = Crv[kp]
                if rowSupp[row] != i
                    rowSupp[row] = i
                    newrowscounter +=1
                    rowList[newrowscounter] = row
                    rowvalCj[row] = true
                else
                    rowvalCj[row] = !rowvalCj[row]
                end
            end
        end
        ii = i+1
        Mcp[ii] = Mcp[i]
        for k = 1:newrowscounter
            row = rowList[k]
            if rowvalCj[row]
                Mrv[Mcp[ii]] = row
                Mcp[ii]+=1
            end
        end
    end
end

