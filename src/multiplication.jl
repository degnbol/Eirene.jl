#!/usr/bin/env julia

function addcol!(
        oddfloods::BitArray{1},shoreline::Array{Int,1},watermark::Int,
        peakcounter::Array{Int,1},Acp::Array{Int,1},Arv::Array{Int,1},
        flippedlist::Array{Int,1},j::Int)
    for i = cran(Acp,j)
        ii = Arv[i]
        if shoreline[ii] != watermark
            peakcounter[1]+=1
            flippedlist[peakcounter]=ii
            shoreline[ii] = watermark
            oddfloods[ii] = true
        else
            oddfloods[ii] = !oddfloods[ii]
        end
    end
end

function spmmF2(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am) where Tv<:Integer
    mA = Am
    nB = length(Bcolptr)-1
    rowvalA = Arowval; colptrA = Acolptr
    rowvalB = Browval; colptrB = Bcolptr
    preallocationIncrement = colptrA[end]+colptrB[end]

    colptrC = Array{Tv}(undef,nB+1)
    colptrC[1] = 1
    rowSupp = zeros(Tv, mA)
    rowList = Array{Tv}(undef,mA)
    rowvalCj = Array{Bool}(undef,mA)
    rowvalC = Array{Tv}(undef,preallocationIncrement)
    for i in 1:nB
        newrowscounter = 0
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

        if colptrC[i+1] > length(rowvalC)+1
            append!(rowvalC, Array{Int}(undef,preallocationIncrement))
        end
        rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
    end
    deleteat!(rowvalC,colptrC[end]:length(rowvalC))
    return rowvalC, colptrC
end

function spmmF2_testfun()
    N =	1000; # matrix dimension
    n = 100;  # number of samples
    for p = 1:n
        m1 =	sprand(N,N,0.01); rv1 =	m1.rowval;  cp1 =	m1.colptr;
        m2 =	sprand(N,N,0.01);	rv2 =	m2.rowval;		cp2 =	m2.colptr;
        if !isempty(rv1)
            m1m2 = ceil.(Int,m1)*ceil.(Int,m2)
            m1m2 = mod.(m1m2,2)
            qrv = m1m2.rowval;
            qcp = m1m2.colptr;
            prv,pcp = spmmF2(rv1,cp1,rv2,cp2,maximum(rv1))
            for q = 1:N
                if crows(qcp,qrv,q) != sort(crows(pcp,prv,q))
                    print("error 1: please see spmmF2_testfun")
                    return
                end
            end
        end
        prv,pcp = spmmF2(zeros(Int,0),ones(Int,N+1),rv2,cp2,1)
        if (prv !=	zeros(Int,0)) || (pcp != ones(Int,N+1))
            print("error 2: please see spmmF2_testfun")
            return
        end
    end
    print("test complete - no errors detected")
end

function spmmF2silentLeft(Arowval::Array{Tv,1},Acolptr::Array{Tv,1},Browval::Array{Tv,1},Bcolptr::Array{Tv,1},Am) where Tv<:Integer
    mA = Am
    nB = length(Bcolptr)-1
    rowvalA = Arowval; colptrA = Acolptr
    rowvalB = Browval; colptrB = Bcolptr
    preallocationIncrement = colptrA[end]+colptrB[end]

    colptrC = Array{Tv}(undef,nB+1)
    colptrC[1] = 1
    rowSupp = zeros(Tv,mA)
    rowList = Array{Tv}(undef,mA)
    rowvalCj = Array{Bool}(undef,mA)
    rowvalC = Array{Tv}(undef,preallocationIncrement)
    for i in 1:nB
        newrowscounter = 0
        eyerange = cran(Bcolptr,i)
        newrowscounter = length(eyerange)
        for l=1:newrowscounter
            ll = rowvalB[eyerange[l]]
            rowList[l] = ll
            rowSupp[ll] = i
            rowvalCj[ll]=true
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

        if colptrC[i+1] > length(rowvalC)+1
            append!(rowvalC, Array{Int}(undef,preallocationIncrement))
        end
        rowvalC[colptrC[i]:(colptrC[i+1]-1)] = sort(rowList[nzRows])
    end
    deleteat!(rowvalC,colptrC[end]:length(rowvalC))
    return rowvalC, colptrC
end

function blockprodsumWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement)
    extend!(Mcp,Dn+1)
    deleteat!(Mcp,(Dn+2):length(Mcp))
    Mcp[1]=1

    rowSupp = zeros(Int, Dm)
    rowList = Array{Int}(undef, Dm)
    rowvalCj = BitArray(undef,Dm)
    for i in 1:Dn
        if length(Mrv)<Mcp[i]+Dm
            extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
        end
        eyerange = cran(Dcp,i)
        newrowscounter = 0
        for kp = eyerange
            row = Drv[kp]
            rowSupp[row]=i
            newrowscounter+=1
            rowList[newrowscounter]=row
            rowvalCj[row] = true
        end
        if Bcp[i]<Bcp[Dn+1]
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
        end
        ii = i+1
        Mcp[ii]=Mcp[i]
        for k = 1:newrowscounter
            row = rowList[k]
            if rowvalCj[row]
                Mrv[Mcp[ii]]=row
                Mcp[ii]+=1
            end
        end
    end
end

function blockprodsumsilenticolsleftWrite2!(Crv,Ccp,Brv,Bcp,Drv,Dcp,Dm,Dn,Mrv,Mcp,preallocationIncrement,col2silenti)
    extend!(Mcp,Dn+1)
    deleteat!(Mcp,(Dn+2):length(Mcp))
    Mcp[1]=1

    rowSupp = zeros(Int, Dm)
    rowList = Array{Int}(undef,Dm)
    rowvalCj = BitArray(undef,Dm)
    for i in 1:Dn
        if length(Mrv)<Mcp[i]+Dm
            extend!(Mrv,length(Mrv)+Dm+preallocationIncrement)
        end
        eyerange = cran(Dcp,i)
        newrowscounter = 0
        for kp = eyerange
            row = Drv[kp]
            rowSupp[row]=i
            newrowscounter+=1
            rowList[newrowscounter]=row
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
        Mcp[ii]=Mcp[i]
        for k = 1:newrowscounter
            row = rowList[k]
            if rowvalCj[row]
                Mrv[Mcp[ii]]=row
                Mcp[ii]+=1
            end
        end
    end
end

