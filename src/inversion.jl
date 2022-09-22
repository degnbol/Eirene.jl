#!/usr/bin/env julia

function morseInverseF2orderedColsUnsortedRowsInSilentOut(Arowval::Array{Tv,1},Acolptr::Array{Tv,1}) where Tv<:Integer
    mA = length(Acolptr)-1
    colptrA = Acolptr
    rowvalA = Arowval
    preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(undef,mA+1); colptrC[1]=1
    rowSupp = zeros(Tv, mA)
    rowList = Array{Tv}(undef,mA)
    rowvalCj = Array{Bool}(undef,mA)
    rowvalC = Array{Tv}(undef,mA)
    totalrowscounter = 0
    onepast = 0
    for i in 1:mA
        if colptrC[i]+mA > length(rowvalC)+1
            append!(rowvalC, Array{Int}(undef,preallocationIncrement))
        end
        if colptrA[i]+1 == colptrA[i+1]
            colptrC[i+1] = colptrC[i]
        elseif colptrA[i]+2==colptrA[i+1]
            if Arowval[colptrA[i]]<i
                k = Arowval[colptrA[i]]
            else
                k = Arowval[colptrA[i]+1]
            end
            ccap = colptrC[i]+colptrC[k+1]-colptrC[k]-1
            rowvalC[colptrC[i]:ccap]= rowvalC[colptrC[k]:(colptrC[k+1]-1)]
            rowvalC[ccap+1]=k
            colptrC[i+1]=ccap+2
        else
            eyerange = cran(Acolptr,i)
            newRowsCounter = length(eyerange)
            for l=1:newRowsCounter
                ll = Arowval[eyerange[l]]
                rowList[l] = ll
                rowSupp[ll] = i
                rowvalCj[ll]=true
            end
            rowvalCj[i] = false ## note we have to make this correction
            for jp in eyerange
                j = rowvalA[jp]
                if j < i
                    for kp in colptrC[j]:(colptrC[j+1] - 1)
                        k = rowvalC[kp]
                        if rowSupp[k] != i
                            rowSupp[k] = i
                            newRowsCounter +=1
                            rowList[newRowsCounter] = k
                            rowvalCj[k] = true
                        else
                            rowvalCj[k] = !rowvalCj[k]
                        end
                    end
                end
            end
            marker = colptrC[i]
            for l = 1:newRowsCounter
                if rowvalCj[rowList[l]]
                    rowvalC[marker]=rowList[l]
                    marker+=1
                end
            end
            colptrC[i+1]=marker
        end

    end
    deleteat!(rowvalC,colptrC[end]:length(rowvalC))
    return rowvalC, colptrC
end

function morseInverseF2orderedColsUnsortedRowsSilentInSilentOut(Arowval::Array{Tv,1},Acolptr::Array{Tv,1}) where Tv<:Integer
    mA = length(Acolptr)-1
    colptrA = Acolptr
    rowvalA = Arowval
    preallocationIncrement = colptrA[end]

    colptrC = Array{Tv}(undef,mA+1); colptrC[1]=1
    rowSupp = zeros(Tv, mA)
    rowList = Array{Tv}(undef,mA)
    rowvalCj = Array{Bool}(undef,mA)
    rowvalC = Array{Tv}(undef,mA)
    totalrowscounter = 0
    onepast = 0
    for i in 1:mA
        if colptrC[i]+mA > length(rowvalC)+1
            append!(rowvalC, Array{Int}(undef,preallocationIncrement))
        end
        if colptrA[i] == colptrA[i+1]
            colptrC[i+1] = colptrC[i]
        elseif colptrA[i]+1==colptrA[i+1]
            k = Arowval[colptrA[i]]
            ccap = colptrC[i]+colptrC[k+1]-colptrC[k]
            rowvalC[colptrC[i]:(ccap-1)]= crows(colptrC,rowvalC,k)
            rowvalC[ccap]=k
            colptrC[i+1]=ccap+1
        else
            eyerange = cran(Acolptr,i)
            newRowsCounter = length(eyerange)
            for l=1:newRowsCounter
                ll = Arowval[eyerange[l]]
                rowList[l] = ll
                rowSupp[ll] = i
                rowvalCj[ll]=true
            end
            for jp in eyerange
                for kp in cran(colptrC,rowvalA[jp])
                    k = rowvalC[kp]
                    if rowSupp[k] != i
                        rowSupp[k] = i
                        newRowsCounter +=1
                        rowList[newRowsCounter] = k
                        rowvalCj[k] = true
                    else
                        rowvalCj[k] = !rowvalCj[k]
                    end
                end
            end
            marker = colptrC[i]
            for l = 1:newRowsCounter
                if rowvalCj[rowList[l]]
                    rowvalC[marker]=rowList[l]
                    marker+=1
                end
            end
            colptrC[i+1]=marker
        end
    end
    deleteat!(rowvalC,colptrC[end]:length(rowvalC))
    return rowvalC, colptrC
end

