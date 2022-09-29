#!/usr/bin/env julia

ocff2of(grain::Vector{Int}, ocg2rad::Vector{Int}) = ocg2rad[grain]
ocff2of(grain::Vector{Int}, ocg2rad::Vector{Float64}) = ocg2rad[grain]
function ocff2of(grain::Vector{Vector{Int}}, ocg2rad::Vector{Float64})
    [ocg2rad[g] for g in grain]
end

"""
NB: eirene permutes vertex labels prior to calculation;
all versions of <vertexrealization> must be interpreted with
respect to the PERMUTED labeling scheme
"""
function vertexrealization(farfaces::Vector{Vector{Int}}, firstv::Vector{Vector{Int}}, facecardinality::Int, facenames::Vector{Int})::Matrix{Int}
    numfaces = length(facenames)
    m = length(firstv[2]) - 1
    loci = copy(facenames)
    vrealization = Matrix{Int}(undef, facecardinality, numfaces)
    post0 = 1
    post1 = 1

    for sd = facecardinality:-1:1
        cp = firstv[sd]
        for i = 1:numfaces
            locus = loci[i]
            if cp[post0] > locus
                while cp[post0] > locus
                    post0 -= 1
                end
                post1 = post0 + 1
            elseif cp[post1] <= locus
                while cp[post1] <= locus
                    post1 += 1
                end
                post0 = post1-1
            end
            loci[i] = farfaces[sd][locus]
            vrealization[sd,i] = post0
        end
    end
    vrealization
end

function buildclosefromclose(lrowval::Vector{Int}, lcolptr, lclosefaces, hrowval, hcolptr; facecard::Int=size(lclosefaces,1)+1)
    m = length(hcolptr) - 1
    n = length(hrowval)
    hclosefaces = Matrix{Int}(undef, facecard, n)
    n == 0 && return hclosefaces
    rowdepth = facecard - 1
    rosettacol = Vector{Int}(undef, maximum(lrowval))
    for i = 1:m
        rosettacol[lrowval[cran(lcolptr,i)]]=cran(lcolptr,i)
        for j = cran(hcolptr,i)
            farface = hrowval[j]
            for k = 1:rowdepth
                hclosefaces[k,j] = rosettacol[lclosefaces[k,farface]]
            end
            hclosefaces[facecard,j] = rosettacol[lrowval[farface]]
        end
    end
    hclosefaces
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
    hclosefaces = Matrix{Int}(undef, sd+1, numselected)
    numselected == 0 && return hclosefaces
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
    hclosefaces
end

function buildclosefromfar(farfaces, firstv, sd, columnsinorder)::Matrix{Int}
    m = length(firstv[1]) - 1
    n = length(farfaces[sd])
    sd == 1 && return Matrix{Int}(undef,0,m)
    lclosefaces = Matrix{Int}(undef, 1, firstv[2][end]-1)
    for i = 1:m
        lclosefaces[cran(firstv[2],i)] .= i
    end
    sd == 2 && return lclosefaces[columnsinorder]'
    for i = 3:sd-1
        lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
    end
    lclosefaces = buildclosefromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder;facecard = sd-1)
end

function buildallfromfar(farfaces,firstv,sd,columnsinorder)::Matrix{Int}
    m = length(firstv[1])-1
    n = length(farfaces[sd])
    if sd == 1
        return Matrix{Int}(undef,0,m)
    end
    lclosefaces = Matrix{Int}(undef, 1, firstv[2][end]-1)
    for i = 1:m
        lclosefaces[cran(firstv[2],i)] .= i
    end
    if sd == 2
        return vcat(lclosefaces[columnsinorder]', farfaces[sd][columnsinorder]')
    end
    for i = 3:(sd-1)
        lclosefaces = buildclosefromclose(farfaces[i-1],firstv[i-1],lclosefaces,farfaces[i],firstv[i];facecard=i-1)
    end
    lclosefaces = buildallfromclose(farfaces[sd-1],firstv[sd-1],lclosefaces,farfaces[sd],firstv[sd],columnsinorder)
end

function ff2aflight_sc2(farfaces::Vector{Vector{Int}}, firstv::Vector{Vector{Int}}, columns)
    sd = 2
    isempty(farfaces[sd]) && return Matrix{Int}(undef,2,0)
    
    f0faces = farfaces[sd]
    colptr = firstv[2]
    columnpost = 1
    columnpostp1 = 2
    faces = Matrix{Int}(undef, 2, length(columns))

    for fp = 1:length(columns)
        f0 = columns[fp]
        if f0 >= colptr[columnpostp1]
            while f0 >= colptr[columnpostp1]
                columnpostp1 += 1
            end
            columnpost = columnpostp1 - 1
        elseif f0 < colptr[columnpost]
            while f0 < colptr[columnpost]
                columnpost -= 1
            end
            columnpostp1 = columnpost + 1
        end
        faces[1,fp] = columnpost
        faces[2,fp] = f0faces[f0]
    end
    faces
end

function ff2aflight_sc3(farfaces::Vector{Vector{Int}}, firstv::Vector{Vector{Int}}, columns)
    sd = 3
    isempty(farfaces[sd]) && return Array{Int}(undef,3,0)

    fcfaces = buildclosefromfar(farfaces, firstv, sd-1, 1:length(farfaces[2]))

    f0faces = farfaces[sd]
    f1faces = farfaces[sd-1]

    fvscm0  = firstv[sd]
    fvscm1  = firstv[sd-1]
    fvscm2  = firstv[sd-2]

    holdi = [1]
    holdip1 = [2]
    t1 = Vector{Int}(undef, fvscm2[end]-1)
    t1[crows(fvscm1,f1faces,1)] = cran(fvscm1,1)

    faces = Matrix{Int}(undef, 3, length(columns))
    for fp = 1:length(columns)
        f0 = columns[fp]
        f1 = f0faces[f0]
        f2 = f1faces[f1]
        f3 = fcfaces[f1]
        updatetranslator!(f0,fvscm0, holdi, holdip1, t1, fvscm1, f1faces)
        faces[1,fp] = t1[f3]
        faces[2,fp] = t1[f2]
        faces[3,fp] = f1
    end
    faces
end

function ff2aflight_scgt3(farfaces::Vector{Vector{Int}}, firstv::Vector{Vector{Int}}, sd::Int, columns)
    isempty(farfaces[sd]) && return Matrix{Int}(undef, sd, 0)

    f0faces = farfaces[sd]
    f1faces = farfaces[sd-1]
    f2faces = farfaces[sd-2]
    fcfaces = buildallfromfar(farfaces,firstv,sd-2,1:(firstv[sd-2][end]-1))::Matrix{Int}

    fvscm0 = firstv[sd]
    fvscm1 = firstv[sd-1]
    fvscm2 = firstv[sd-2]
    fvscm3 = firstv[sd-3]

    holdi=[1];holdip1=[2];holdj=[1];holdjp1=[2]
    t1 = Vector{Int}(undef,fvscm2[end]-1)
    t1[crows(fvscm1,f1faces,1)]=cran(fvscm1,1)
    t2 = Vector{Int}(undef,fvscm3[end]-1)
    t2[crows(fvscm2,f2faces,1)]=cran(fvscm2,1)

    faces = Matrix{Int}(undef, sd, length(columns))
    for fp = 1:length(columns)
        f0 = columns[fp]
        f1 = f0faces[f0]
        f2 = f1faces[f1]
        updatetranslator!(f0,fvscm0,holdi,holdip1,t1,fvscm1,f1faces)
        updatetranslator!(f1,fvscm1,holdj,holdjp1,t2,fvscm2,f2faces)
        for i = 1:sd-2
            faces[i,fp] = t1[t2[fcfaces[i,f2]]]
        end
        faces[sd-1, fp] = t1[f2]
        faces[sd, fp] = f1
    end
    faces
end

function updatetranslator!(f0::Int, firstv0::Vector{Int}, holdi::Vector{Int}, holdip1::Vector{Int}, t::Vector{Int}, firstv1::Vector{Int}, farfaces1::Vector{Int})
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


function ff2aflight(farfaces, firstv, sd, columns)
    if sd == 1
        return Matrix{Int}(undef, 0, length(columns))
    elseif sd == 2
        return ff2aflight_sc2(farfaces,firstv,columns)
    elseif sd == 3
        return ff2aflight_sc3(farfaces,firstv,columns)
    else
        return ff2aflight_scgt3(farfaces,firstv,sd,columns)
    end
end

"""
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
"""
function filteredmatrixfromfarfaces(farfaces, firstv, prepairs, grain, sd::Int, lowbasisnames::Vector{Int})
    numhigs = length(farfaces[sd])
    numlows = length(farfaces[sd-1])
    numppair = length(prepairs[sd])

    pphigs = prepairs[sd]
    pplows = farfaces[sd][pphigs]
    lpls = lowbasisnames
    hphs = farfaces[sd+1][prepairs[sd+1]]
    nplows = intervalcomplementuniqueunsortedinput(vcat(lpls,pplows), numlows)
    nphigs = intervalcomplementuniqueunsortedinput(vcat(hphs,pphigs), numhigs)

    numnhph = numhigs - length(hphs)

    higtranslator = zeros(Int, numnhph)
    lowtranslator = zeros(Int, numlows)
    lowtranslator[pplows] = 1:numppair

    if sd > 1
        npfilt = grain[sd-1][nplows]
        nporder = integersinsameorder(npfilt)
        nporder .+= numppair
    else
        npfilt = zeros(Int, 0)
        nporder = zeros(Int, 0)
    end

    lowtranslator[nplows] = nporder
    higsinpointorder = intervalcomplementuniqueunsortedinput(hphs,numhigs)
    lowlab = Vector{Int}(undef, numlows - length(lpls))
    lowlab[1:numppair] = pplows
    lowlab[nporder] = nplows
    higlab = vcat(pphigs,nphigs)

    ppsupp = falses(numhigs)
    ppsupp[pphigs] .= true
    ppmarker = 0
    nppmarker = numppair
    for i = 1:numnhph
        hig = higsinpointorder[i]
        if ppsupp[hig]
            ppmarker += 1
            higtranslator[i] = ppmarker
        else
            nppmarker += 1
            higtranslator[i] = nppmarker
        end
    end
    allfaces = buildallfromfar(farfaces,firstv,sd,higsinpointorder)
    Mrv, Mcp, Mm = presparsefull2unsortedsparsetranspose(allfaces, lowtranslator, higtranslator)
    return Mrv, Mcp, lowlab, higlab, Mm
end

"""
Get the values in range 1:n not found in v.
- v: unique positive integers
- n: interval endpoint
"""
function intervalcomplementuniqueunsortedinput(v::Vector{Int}, n::Int)
    complementsupport = trues(n)
    complementsupport[v] .= false
    findall(complementsupport)
end



