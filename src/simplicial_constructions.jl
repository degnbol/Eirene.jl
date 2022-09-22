#!/usr/bin/env julia

ocff2of(grain::Array{Int}, ocg2rad::Array{Int}) = ocg2rad[grain]
ocff2of(grain::Array{Int}, ocg2rad::Array{Float64}) = ocg2rad[grain]
function ocff2of(grain::Array{Array{Int,1},1}, ocg2rad::Array{Float64})
    [ocg2rad[g] for g in grain]
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(farfaces, firstv, facecardinality::Int, facenames)
    numfaces = length(facenames)
    m = length(firstv[2]) - 1
    preallocationspace = 0
    loci::Array{Int,1} = copy(facenames)
    vrealization = Array{Int}(undef,facecardinality,numfaces)
    post0::Int = 1
    post1::Int = 1

    for sd = facecardinality:-1:1
        cp::Array{Int,1} = firstv[sd]
        for i = 1:numfaces
            locus = loci[i]
            if cp[post0] > locus
                while cp[post0] > locus
                    post0-=1
                end
                post1 = post0+1
            elseif cp[post1] <= locus
                while cp[post1] <= locus
                    post1+=1
                end
                post0 = post1-1
            end
            loci[i] = farfaces[sd][locus]
            vrealization[sd,i] = post0
        end
    end
    return vrealization
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict, facecardinality::Int, facenames)
    return vertexrealization(D["farfaces"], D["firstv"], facecardinality, facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# all versions of <vertexrealization> must be interpreted with
# respect to the PERMUTED labeling scheme
function vertexrealization(D::Dict; dim::Int=1, class::Int=1)
    sd = dim+2
    facecard = dim+1

    if haskey(D,"cyclerep")
        rep = D["cyclerep"][sd][class]
    else
        cyclename = barname2cyclename(D,class; dim=dim)
        rep = getcycle(D,sd,cyclename)
    end

    return vertexrealization(D::Dict, facecard, rep)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(farfaces, firstv, facecardinality, facenames)
    numfaces::Int = length(facenames)
    m::Int = length(firstv[2]) - 1
    preallocationspace = 0
    vsupp = falses(m)
    loci::Array{Int,1} = copy(facenames)
    post0::Int = 1
    post1::Int = 1
    for sd = facecardinality:-1:1
        cp::Array{Int,1} = firstv[sd]
        for i = 1:numfaces
            locus = loci[i]
            if cp[post0] > locus
                while cp[post0] > locus
                    post0-=1
                end
                post1 = post0+1
            elseif cp[post1] <= locus
                while cp[post1] <= locus
                    post1+=1
                end
                post0 = post1-1
            end
            loci[i] = farfaces[sd][locus]
            vsupp[post0] = true
        end
    end
    return findall(vsupp)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict, facecardinality, facenames)
    return incidentverts(D["farfaces"],D["firstv"],facecardinality,facenames)
end

# NB: eirene permutes vertex labels prior to calculation;
# the output of <incidentverts> must be interpreted with respect
# to the PERMUTED labeling scheme
function incidentverts(D::Dict; dim::Int=1, class::Int=1)
    facecardinality = dim+1

    if haskey(D,"cyclerep")
        rep = D["cyclerep"][dim+2][class]
    else
        cyclename = barname2cyclename(D,class;dim=dim)
        rep = getcycle(D,facecardinality,class)
    end
    return incidentverts(D,facecardinality,rep)
end

function buildclosefromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr;facecard = size(lclosefaces,1)+1)
    m = length(hcolptr)-1
    n = length(hrowval)
    hclosefaces = Array{Int}(undef,facecard,n)
    if n == 0
        return hclosefaces
    else
        rowdepth = facecard-1
        rosettacol = Array{Int}(undef,maximum(lrowval))
        for i = 1:m
            rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
            for j = cran(hcolptr,i)
                farface = hrowval[j]
                for k = 1:rowdepth
                    hclosefaces[k,j]=rosettacol[lclosefaces[k,farface]]
                end
                hclosefaces[facecard,j] = rosettacol[lrowval[farface]]
            end
        end
        return hclosefaces
    end
end

"""
PLEASE NOTE: COLUMNS MUST BE IN SORTED ORDER FOR THIS TO WORK PROPERLY.
"""
function buildallfromclose(lrowval,lcolptr,lclosefaces,hrowval,hcolptr,selectedcolumnindices)
    m = length(hcolptr)-1
    numhigs = length(hrowval)
    numselected = length(selectedcolumnindices)
    rowdepth = size(lclosefaces,1)
    sd = rowdepth+1
    hclosefaces = Array{Int}(undef,sd+1,numselected)
    if numselected == 0
        return hclosefaces
    end
    rosettacol = Array{Int}(undef,maximum(lrowval))
    columnsupp = falses(numhigs)
    columnsupp[selectedcolumnindices].=true
    columnmarker = 0
    for i = 1:m
        rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
        for j = cran(hcolptr,i)
            if columnsupp[j]
                columnmarker+=1
                farface = hrowval[j]
                for k = 1:rowdepth
                    hclosefaces[k,columnmarker]=rosettacol[lclosefaces[k,farface]]
                end
                hclosefaces[sd,columnmarker] = rosettacol[lrowval[farface]]
                hclosefaces[sd+1,columnmarker] = farface
            end
        end
    end
    return hclosefaces
end

function buildclosefaces!(lrowval,lcolptr,lclosefaces,lfarfaces,hrowval,hcolptr,destinationmatrix)
    m = length(hcolptr)-1
    n = length(hrowval)
    rowdepth = size(lclosefaces,1)
    sd = rowdepth+1
    rosettacol = Array{Int}(undef,maximum(lrowval))
    for i = 1:m
        rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
        for j = cran(hcolptr,i)
            farface = hrowval[j]
            for k = 1:rowdepth
                destinationmatrix[k,j]=rosettacol[lclosefaces[farface]]
            end
            destinationmatrix[sd,j] = rosettacol[lrowval[farface]]
        end
    end
    for j = 1:n
        for i = 1:sd
            lclosefaces[i,j]=destinationmatrix[i,j]
        end
    end
end

function buildclosefromfar(farfaces,firstv,sd)
    m = length(firstv[1])-1
    n = length(farfaces[sd])
    # destinationmatrix = Array{Int}(undef,sd,n)
    if sd == 1
        return Array{Int}(undef,0,m)
    end
    lclosefaces = Array{Int}(undef,1,firstv[2][end]-1)
    for i = 1:m
        lclosefaces[cran(firstv[2],i)]=i
    end
    if sd == 2
        return lclosefaces'
    end
    for i = 3:sd
        lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
    end
    return lclosefaces
end

function buildclosefromfar(farfaces,firstv,sd,columnsinorder)
    m = length(firstv[1])-1
    n = length(farfaces[sd])
    if sd == 1
        return Array{Int}(undef,0,m)
    end
    lclosefaces = Array{Int}(undef,1,firstv[2][end]-1)
    for i = 1:m
        lclosefaces[cran(firstv[2],i)].=i
    end
    if sd == 2
        return lclosefaces[columnsinorder]'
    end
    for i = 3:(sd-1)
        lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
    end
    lclosefaces = buildclosefromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;facecard = sd-1)
    return lclosefaces
end

function buildallfromfar(farfaces,firstv,sd,columnsinorder)
    m = length(firstv[1])-1
    n = length(farfaces[sd])
    # destinationmatrix = Array{Int}(undef,sd,n)
    if sd == 1
        return Array{Int}(undef,0,m)
    end
    lclosefaces = Array{Int}(undef,1,firstv[2][end]-1)
    for i = 1:m
        lclosefaces[cran(firstv[2],i)].=i
    end
    if sd == 2
        return vcat(lclosefaces[columnsinorder]',farfaces[sd][columnsinorder]')
    end
    for i = 3:(sd-1)
        lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
        #gc()
    end
    lclosefaces = buildallfromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder)
    #gc()
    return lclosefaces
end

function ff2boundary(farfaces,firstv;sd=1)
    rv = Array{Int}(undef,0)
    cp = [1]
    if sd == 1
        rv = Array{Int}(undef,0)
        cp = ones(Int,length(farfaces[1])+1)
    else
        n = length(farfaces[sd])
        rv = ff2aflight(farfaces,firstv,sd,1:n)
        rv = vec(rv)
        cp = convert(Array{Int,1},sd+1:sd:(1+n*sd))
        prepend!(cp,[1])
    end
    return rv,cp
end

function ff2complex(farfaces,firstv;maxsd = length(farfaces))
    Nrv = fill(Array{Int}(undef,0),maxsd)
    Ncp = fill(Array{Int}(undef,0),maxsd)
    Nrv		= convert(Array{Array{Int,1}},Nrv)
    Ncp		= convert(Array{Array{Int,1}},Ncp)
    Nrv[1] = Array{Int}(undef,0)
    Ncp[1]	= fill(1,length(farfaces[1])+1)
    for sd = 2:maxsd
        Nrv[sd],Ncp[sd] = ff2boundary(farfaces,firstv,sd=sd)
    end
    return Nrv,Ncp
end

function eirened2complex(C)
    if in(C["input"]["model"],["pc","vr"])
        rv,cp = boundarymatrices(C)
        fv = ocff2of(C["grain"],C["ocg2rad"])
    elseif in(C["input"]["model"],["complex"])
        rv = C["rv"]
        cp = C["cp"]
    else
        println("Error: the value of C[\"input\"][\"model\"] must be \"pc\", \"vr\", or \"complex\".")
    end
    fv = ocff2of(C["grain"],C["ocg2rad"])
    return rv,cp,fv
end


function ff2aflight_sc2(farfaces,firstv,columns)
    sd = 2
    if isempty(farfaces[sd])
        return Array{Int}(undef,2,0)
    end
    f0faces::Array{Int,1} = farfaces[sd]
    colptr::Array{Int,1} = firstv[2]
    columnpost::Int   = 1
    columnpostp1::Int = 2
    faces::Array{Int,2} = Array{Int}(undef,2,length(columns))

    for fp = 1:length(columns)
        f0 = columns[fp]
        if f0 >= colptr[columnpostp1]
            while f0 >= colptr[columnpostp1]
                columnpostp1+=1
            end
            columnpost = columnpostp1-1
        elseif f0 < colptr[columnpost]
            while f0 < colptr[columnpost]
                columnpost-=1
            end
            columnpostp1 = columnpost+1
        end
        faces[1,fp] = columnpost
        faces[2,fp] = f0faces[f0]
    end
    return faces
end

function ff2aflight_sc3(farfaces,firstv,columns)
    sd = 3

    if isempty(farfaces[sd])
        return Array{Int}(undef,3,0)
    end

    fcfaces::Array{Int,2} = buildclosefromfar(farfaces,firstv,sd-1,1:length(farfaces[2]))

    f0faces::Array{Int,1} = farfaces[sd]
    f1faces::Array{Int,1} = farfaces[sd-1]

    fvscm0::Array{Int,1}  = firstv[sd]
    fvscm1::Array{Int,1}  = firstv[sd-1]
    fvscm2::Array{Int,1}  = firstv[sd-2]

    holdi=[1];holdip1=[2]
    t1::Array{Int,1} = Array{Int}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)

    faces::Array{Int,2} = Array{Int}(undef,3,length(columns))
    for fp = 1:length(columns)
        f0 = columns[fp]
        f1 = f0faces[f0]
        f2 = f1faces[f1]
        f3 = fcfaces[f1]
        updatetranslator!(f0::Int,fvscm0::Array{Int,1} ,holdi::Array{Int,1},holdip1::Array{Int,1},t1::Array{Int,1},fvscm1::Array{Int,1},f1faces::Array{Int,1})
        faces[1,fp] = t1[f3]
        faces[2,fp] = t1[f2]
        faces[3,fp] = f1
    end
    return faces
end

function ff2aflight_scgt3(farfaces,firstv,sd,columns)
    isempty(farfaces[sd]) && return Array{Int}(undef,sd,0)

    f0faces::Array{Int,1} = farfaces[sd]
    f1faces::Array{Int,1} = farfaces[sd-1]
    f2faces::Array{Int,1} = farfaces[sd-2]
    fcfaces::Array{Int,2} = buildallfromfar(farfaces,firstv,sd-2,1:(firstv[sd-2][end]-1))

    fvscm0::Array{Int,1}  = firstv[sd]
    fvscm1::Array{Int,1}  = firstv[sd-1]
    fvscm2::Array{Int,1}  = firstv[sd-2]
    fvscm3::Array{Int,1}  = firstv[sd-3]

    holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
    t1::Array{Int,1} = Array{Int}(undef,fvscm2[end]-1);t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
    t2::Array{Int,1} = Array{Int}(undef,fvscm3[end]-1);t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

    scm0::Int = sd; scm1::Int = sd-1; scm2::Int = sd-2
    faces::Array{Int,2} = Array{Int}(undef,sd,length(columns))

    for fp = 1:length(columns)
        f0 = columns[fp]
        f1 = f0faces[f0]
        f2 = f1faces[f1]
        updatetranslator!(f0::Int,fvscm0::Array{Int,1},holdi::Array{Int,1},holdip1::Array{Int,1},t1::Array{Int,1},fvscm1::Array{Int,1},f1faces::Array{Int,1})
        updatetranslator!(f1::Int,fvscm1::Array{Int,1},holdj::Array{Int,1},holdjp1::Array{Int,1},t2::Array{Int,1},fvscm2::Array{Int,1},f2faces::Array{Int,1})
        for i = 1:scm2
            faces[i,fp] = t1[t2[fcfaces[i,f2]]]
        end
        faces[scm1,fp] = t1[f2]
        faces[scm0,fp] = f1
    end
    faces
end

function updatetranslator!(f0,firstv0,holdi,holdip1,t,firstv1,farfaces1)
    if firstv0[holdip1[1]] <= f0
        while firstv0[holdip1[1]] <= f0
            holdip1[1]+=1
        end
        holdi[1] = holdip1[1]-1
        t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
    elseif firstv0[holdi[1]] > f0
        while firstv0[holdi[1]] > f0
            holdi[1]-=1
        end
        holdip1[1] = holdi[1]+1
        t[crows(firstv1,farfaces1,holdi[1])]=cran(firstv1,holdi[1])
    end
end

function ff2aflight_subr!(
        columns::UnitRange{Int},f0faces::Array{Int,1},f1faces::Array{Int,1},f2faces::Array{Int,1},fcfaces::Array{Int,2},fvscm0::Array{Int,1},
        fvscm1::Array{Int,1},fvscm2::Array{Int,1},holdi::Array{Int,1},holdip1::Array{Int,1},holdj::Array{Int,1},
        holdjp1::Array{Int,1},t1::Array{Int,1},t2::Array{Int,1},faces::Array{Int,2},
        scm0::Int,scm1::Int,scm2::Int)
    for fp = 1:length(columns)
        f0 = columns[fp]
        f1 = f0faces[f0]
        f2 = f1faces[f1]
        updatetranslator!(f0::Int,fvscm0::Array{Int,1} ,holdi::Array{Int,1},holdip1::Array{Int,1},t1::Array{Int,1},fvscm1::Array{Int,1},f1faces::Array{Int,1})
        updatetranslator!(f1::Int,fvscm1::Array{Int,1},holdj::Array{Int,1},holdjp1::Array{Int,1},t2::Array{Int,1},fvscm2::Array{Int,1},f2faces::Array{Int,1})
        for i = 1:scm2
            faces[i,fp] = t1[t2[fcfaces[i,f2]]]
        end
        faces[scm1,fp] = t1[f2]
        faces[scm0,fp] = f1
    end
end

function ff2aflight(farfaces, firstv, sd, columns)
    if sd == 1
        return Array{Int}(undef,0,length(columns))
    elseif sd == 2
        return ff2aflight_sc2(farfaces,firstv,columns)
    elseif sd == 3
        return ff2aflight_sc3(farfaces,firstv,columns)
    else
        return ff2aflight_scgt3(farfaces,firstv,sd,columns)
    end
end

function ff2aflight(D::Dict, sd, columns)
    farfaces = D["farfaces"]; firstv = D["firstv"]
    faces = ff2aflight(farfaces,firstv,sd,columns)
    return faces
end

#=

NB
- Input argument <grain> must be arranged least to greatest

OUTPUTS
- higlab
The concatenated vector [pphigs,nphigs)]
- lowlab
The concatenated vector [pplows,nplows[perm]], where perm is a permutation such
that the entries of lowgrain[nplows[perm]] appear in ascending order, numer-
ically.
- Mrv, Mcp, Mm
Sparse matrix representation of transpose(D[lowlab,higlab]), where D is
submatrix of the total boundary operator indexed by cells of dimension sd-1
(along the columns) and sd-2 (along the rows).

=#
function filteredmatrixfromfarfaces(
        farfaces,
        firstv,
        prepairs,
        grain,
        sd::Integer,
        lowbasisnames::Array{Int,1})

    numhigs = length(farfaces[sd])
    numlows = length(farfaces[sd-1])
    numppair= length(prepairs[sd])

    pphigs = prepairs[sd]
    pplows = farfaces[sd][pphigs]
    lpls = lowbasisnames
    hphs = farfaces[sd+1][prepairs[sd+1]]
    nplows = intervalcomplementuniqueunsortedinput(vcat(lpls,pplows),numlows)
    nphigs = intervalcomplementuniqueunsortedinput(vcat(hphs,pphigs),numhigs)

    numnhph = numhigs-length(hphs)
    Ml = numlows - length(lpls)
    Mh = numhigs - length(hphs)

    higtranslator = zeros(Int,numnhph)
    lowtranslator = zeros(Int,numlows)
    lowtranslator[pplows] = 1:numppair

    if sd > 1
        npfilt = grain[sd-1][nplows]
        nporder = integersinsameorder(npfilt)
        addinteger!(nporder,numppair)
    else
        npfilt = zeros(Int,0)
        nporder =	zeros(Int,0)
    end

    lowtranslator[nplows] = nporder
    higsinpointorder = intervalcomplementuniqueunsortedinput(hphs,numhigs)
    lowlab = Array{Int}(undef,Ml)
    lowlab[1:numppair]=pplows
    lowlab[nporder]=nplows
    higlab = vcat(pphigs,nphigs)

    if verbose
        comparisonsuppvec = trues(numhigs)
        comparisonsuppvec[hphs]=false
        comparisonvec=findall(comparisonsuppvec)
        differencecounter = 0
        for i = 1:length(higsinpointorder)
            if higsinpointorder[i]!=comparisonvec[i]
                differencecounter+=1
            end
        end
        if differencecounter>0
            print(["hi ho comparison vec" differencecounter])
            print(length(higsinpointorder))
            print(length(comparisonvec))
            print(comparisonvec[1:20])
            print(higsinpointorder[1:20])
            sleep(5)
        end
    end
    ppsupp = falses(numhigs)
    ppsupp[pphigs].=true
    ppmarker = 0
    nppmarker = numppair
    for i = 1:numnhph
        hig = higsinpointorder[i]
        if ppsupp[hig]
            ppmarker+=1
            higtranslator[i]=ppmarker
        else
            nppmarker+=1
            higtranslator[i]=nppmarker
        end
    end
    allfaces = buildallfromfar(farfaces,firstv,sd,higsinpointorder;verbose = verbose)
    if verbose
        print("done building allfromfar")
    end
    Mrv,Mcp,Mm = presparsefull2unsortedsparsetranspose(allfaces,lowtranslator,higtranslator;verbose=verbose)
    higtranslator = [];npfilt = [];ppsupp = [];allfaces = []
    #gc()
    if verbose && length(Mrv)>(Mcp[end]-1)
        print("There was the thought that Mrv should have no extra elements")
        sleep(3)
    end
    return Mrv,Mcp,lowlab,higlab,Mm
end

function getmaxdim(farfaces)
    l = length(farfaces)
    maxdim = l
    for i = 1:l
        if length(farfaces[i]) == 0
            maxdim = i
            break
        end
    end
    return maxdim
end

function grain2maxsd(grain)
    c = 0
    for i = 1:length(grain)
        if !isempty(grain[i])
            c = i
        end
    end
    return c
end

function skelcount(numvertices,maxsdinality)
    c = 0
    for i = 1:maxsdinality
        c += binom(numvertices,i)
    end
    return c
end

