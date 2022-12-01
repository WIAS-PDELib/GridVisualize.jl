"""
$(SIGNATURES)

Create customized distinguishable colormap for interior regions.
For this we use a kind of pastel colors.
"""
region_cmap(n)=distinguishable_colors(max(5,n),
                                      [RGB(0.85,0.6,0.6), RGB(0.6,0.85,0.6),RGB(0.6,0.6,0.85)],
                                      lchoices = range(70, stop=80, length=5),
                                      cchoices = range(25, stop=65, length=15),
                                      hchoices = range(20, stop=360, length=15)
                                      )

"""
$(SIGNATURES)

Create customized distinguishable colormap for boundary regions.

These use fully saturated colors.
"""
bregion_cmap(n)=distinguishable_colors(max(5,n),
                                      [RGB(1.0,0.0,0.0), RGB(0.0,1.0,0.0), RGB(0.0,0.0,1.0)],
                                      lchoices = range(50, stop=75, length=10),
                                      cchoices = range(75, stop=100, length=10),
                                      hchoices = range(20, stop=360, length=30)
                                      )



"""
$(SIGNATURES)

Create RGB color from color name string.
"""
function Colors.RGB(c::String)
    c64=Colors.color_names[c]
    RGB(c64[1]/255,c64[2]/255, c64[3]/255)
end


"""
$(SIGNATURES)

Create RGB color from color name symbol.
"""
Colors.RGB(c::Symbol)=Colors.RGB(String(c))

"""
$(SIGNATURES)

Create RGB color from tuple
"""
Colors.RGB(c::Tuple)=Colors.RGB(c...)


"""
$(SIGNATURES)

Create color tuple from  color description (e.g. string)
"""
rgbtuple(c)=rgbtuple(Colors.RGB(c))


"""
$(SIGNATURES)

Create color tuple from  RGB color.
"""
rgbtuple(c::RGB)=(red(c),green(c),blue(c))




"""
$(SIGNATURES)

Extract visible tetrahedra - those intersecting with the planes
`x=xyzcut[1]` or `y=xyzcut[2]`  or `z=xyzcut[3]`. 

Return corresponding points and facets for each region for drawing as mesh (Makie,MeshCat)
or trisurf (pyplot)
"""
function extract_visible_cells3D(coord,cellnodes,cellregions,nregions,xyzcut;
                                 primepoints=zeros(0,0),Tp=SVector{3,Float32},Tf=SVector{3,Int32})
    
    function take(coord,simplex,xyzcut)
        all_lt=@MVector ones(Bool,3)
        all_gt=@MVector ones(Bool,3)
        for idim=1:3
            for inode=1:4
                c=coord[idim,simplex[inode]]-xyzcut[idim]
                all_lt[idim]=all_lt[idim] && (c<0.0)
                all_gt[idim]=all_gt[idim] && (c>0.0)
            end
        end
        tke=false
        tke=tke  ||   (!all_lt[1])  &&  (!all_gt[1]) && (!all_gt[2]) && (!all_gt[3])
        tke=tke  ||   (!all_lt[2])  &&  (!all_gt[2]) && (!all_gt[1]) && (!all_gt[3])
        tke=tke  ||   (!all_lt[3])  &&  (!all_gt[3]) && (!all_gt[1]) && (!all_gt[2])
    end

    
    faces=[Vector{Tf}(undef,0) for iregion=1:nregions]
    points=[Vector{Tp}(undef,0) for iregion=1:nregions]
    
    for iregion=1:nregions
        for iprime=1:size(primepoints,2)
            @views push!(points[iregion],Tp(primepoints[:,iprime]))
        end
    end
    tet=zeros(Int32,4)
    
    for itet=1:size(cellnodes,2)
        iregion=cellregions[itet]
        for i=1:4
            tet[i]=cellnodes[i,itet]
        end
        if take(coord,tet,xyzcut)
            npts=size(points[iregion],1)
            @views begin
                push!(points[iregion],coord[:,cellnodes[1,itet]])
                push!(points[iregion],coord[:,cellnodes[2,itet]])
                push!(points[iregion],coord[:,cellnodes[3,itet]])
                push!(points[iregion],coord[:,cellnodes[4,itet]])
                push!(faces[iregion],(npts+1,npts+2,npts+3))
                push!(faces[iregion],(npts+1,npts+2,npts+4))
                push!(faces[iregion],(npts+2,npts+3,npts+4))
                push!(faces[iregion],(npts+3,npts+1,npts+4))
            end
        end
    end
    points,faces
end



"""
$(SIGNATURES)

Extract visible boundary faces - those not cut off by the planes
`x=xyzcut[1]` or `y=xyzcut[2]`  or `z=xyzcut[3]`. 

Return corresponding points and facets for each region for drawing as mesh (Makie,MeshCat)
or trisurf (pyplot)
"""
function extract_visible_bfaces3D(coord,bfacenodes,bfaceregions, nbregions, xyzcut;
                                  primepoints=zeros(0,0), Tp=SVector{3,Float32},Tf=SVector{3,Int32})


    nbfaces=size(bfacenodes,2)
    cutcoord=zeros(3)

    function take(coord,simplex,xyzcut)
        for idim=1:3
            all_gt=true
            for inode=1:3
                c=coord[idim,simplex[inode]]-xyzcut[idim]
                all_gt= all_gt && c>0
            end
            if all_gt
                return false
            end
        end
        return true
    end
    

    Tc=SVector{3,eltype(coord)}
    xcoord=reinterpret(Tc,reshape(coord,(length(coord),)))
    
    
    faces=[Vector{Tf}(undef,0) for iregion=1:nbregions]
    points=[Vector{Tp}(undef,0) for iregion=1:nbregions]
    for iregion=1:nbregions
        for iprime=1:size(primepoints,2)
            @views push!(points[iregion],Tp(primepoints[:,iprime]))
        end
    end

    # remove some type instability here
    function collct(points,faces)
        trinodes=[1,2,3]
        for i=1:nbfaces
            iregion=bfaceregions[i]
            trinodes[1]=bfacenodes[1,i]
            trinodes[2]=bfacenodes[2,i]
            trinodes[3]=bfacenodes[3,i]
            if take(coord,trinodes,xyzcut)
                npts=size(points[iregion],1)
                @views push!(points[iregion],xcoord[trinodes[1]])
                @views push!(points[iregion],xcoord[trinodes[2]])
                @views push!(points[iregion],xcoord[trinodes[3]])
                @views push!(faces[iregion],(npts+1,npts+2,npts+3))
            end
        end
    end
    collct(points,faces)
    points,faces
end



"""
  $(SIGNATURES)
  Calculate intersections between tetrahedron with given piecewise linear
  function data and plane 

  Adapted from https://github.com/j-fu/gltools/blob/master/glm-3d.c#L341
 
  A non-empty intersection is either a triangle or a planar quadrilateral,
  define by either 3 or 4 intersection points between tetrahedron edges
  and the plane.

  Input: 
  -       pointlist: 3xN array of grid point coordinates
  -    node_indices: 4 element array of node indices (pointing into pointlist and function_values)
  -   planeq_values: 4 element array of plane equation evaluated at the node coordinates
  - function_values: N element array of function values

  Mutates:
  -  ixcoord: 3x4 array of plane - tetedge intersection coordinates
  - ixvalues: 4 element array of fuction values at plane - tetdedge intersections

  Returns:
  - nxs,ixcoord,ixvalues
  
  This method can be used both for the evaluation of plane sections and for
  the evaluation of function isosurfaces.
"""
function tet_x_plane!(ixcoord,ixvalues,pointlist,node_indices,planeq_values,function_values; tol=0.0)

    # If all nodes lie on one side of the plane, no intersection
    @fastmath if (mapreduce(a->a< -tol,*,planeq_values) || mapreduce(a->a>tol,*,planeq_values))
        return 0
    end
    # Interpolate coordinates and function_values according to
    # evaluation of the plane equation
    nxs=0
    @inbounds @simd for n1=1:3
        N1=node_indices[n1]
        @inbounds @fastmath @simd for n2=n1+1:4
            N2=node_indices[n2]
            if planeq_values[n1]!=planeq_values[n2] && planeq_values[n1]*planeq_values[n2]<tol
                nxs+=1
                t= planeq_values[n1]/(planeq_values[n1]-planeq_values[n2])
                ixcoord[1,nxs]=pointlist[1,N1]+t*(pointlist[1,N2]-pointlist[1,N1])
                ixcoord[2,nxs]=pointlist[2,N1]+t*(pointlist[2,N2]-pointlist[2,N1])
                ixcoord[3,nxs]=pointlist[3,N1]+t*(pointlist[3,N2]-pointlist[3,N1])
                ixvalues[nxs]=function_values[N1]+t*(function_values[N2]-function_values[N1])
            end
        end
    end
    return nxs
end



"""
 We should be able to parametrize this
 with a pushdata function which will remove one copy
 step for GeometryBasics.mesh creation - perhaps a meshcollector struct we
 can dispatch on.
 flevel could be flevels
 xyzcut could be a vector of plane data
 perhaps we can also collect isolines.
 Just an optional collector parameter, defaulting to somethig makie independent.

    Better yet:

 struct TetrahedronMarcher
  ...
 end
 tm=TetrahedronMarcher(planes,levels)

 foreach tet
       collect!(tm, tet_node_coord, node_function_values)
 end
 tm.colors=AbstractPlotting.interpolated_getindex.((cmap,), mcoll.vals, (fminmax,))
 mesh!(collect(mcoll),backlight=1f0) 
""" 



"""
$(SIGNATURES)

Extract isosurfaces and plane interpolation for function on 3D tetrahedral mesh.

The basic observation is that locally on a tetrahedron, cuts with planes and isosurfaces
of P1 functions look the same. This method calculates data for several plane cuts and several
isosurfaces at once. 

Input parameters:
- `coord`: 3 x n_points matrix of point coordinates
- `cellnodes`: 4 x n_cells matrix of point numbers per tetrahedron
- `func`: n_points vector of piecewise linear function values
- `planes`: vector of plane equations `ax+by+cz+d=0`,each  stored as vector [a,b,c,d]
- `flevels`: vector of function isolevels

Keyword arguments:
- `tol`: tolerance for tet x plane intersection
- `primepoints`:  3 x n_prime matrix of "corner points" of domain to be plotted. These are not in the mesh but are used to calculate the axis size e.g. by Makie
- `primevalues`:  n_prime vector of function values in corner points. These can be used to calculate function limits e.g. by Makie
- `Tv`:  type of function values returned
- `Tp`:  type of points returned
- `Tf`:  type of facets returned

Return values: (points, tris, values)
- `points`: vector of points (Tp)
- `tris`: vector of triangles (Tf)
- `values`: vector of function values (Tv)

These can be readily turned into a mesh with function values on it.

Caveat: points with similar coordinates are not identified, e.g. an intersection of a plane and an edge will generate as many edge intersection points as there are tetrahedra adjacent to that edge. As a consequence, normal calculations for visualization alway will end up with facet normals, not point normals, and the visual impression of a rendered isosurface will show its piecewise linear genealogy.

"""
function marching_tetrahedra(coord,cellnodes,func,planes,flevels;
                             tol=1.0e-12,
                             primepoints=zeros(0,0),
                             primevalues=zeros(0),
                             Tv=Float32,
                             Tp=SVector{3,Float32},
                             Tf=SVector{3,Int32})

    # We could rewrite this for Meshing.jl
    # CellNodes::Vector{Ttet}, Coord::Vector{Tpt}
    nplanes=length(planes)
    nlevels=length(flevels)
    nnodes=size(coord,2)
    ntet=size(cellnodes,2)

    all_planeq=Vector{Float32}(undef,nnodes)

    # Create output vectors
    all_ixfaces=Vector{Tf}(undef,0)
    all_ixcoord=Vector{Tp}(undef,0)
    all_ixvalues=Vector{Tv}(undef,0)

    @assert(length(primevalues)==size(primepoints,2))

    for iprime=1:size(primepoints,2)
        @views push!(all_ixcoord,primepoints[:,iprime])
        @views push!(all_ixvalues,primevalues[iprime])
    end
    
    planeq=zeros(4)
    ixcoord=zeros(3,6)
    ixvalues=zeros(6)
    cn=zeros(4)
    node_indices=zeros(Int32,4)

    # Function to evaluate plane equation
    @inbounds @fastmath plane_equation(plane,coord)= coord[1]*plane[1]+coord[2]*plane[2]+coord[3]*plane[3]+plane[4]
    
    function pushtris(ns,ixcoord,ixvalues)
        # number of intersection points can be 3 or 4
        if ns>=3
            last_i=length(all_ixvalues)
            for is=1:ns
                @views push!(all_ixcoord,ixcoord[:,is])
                push!(all_ixvalues,ixvalues[is]) # todo consider nan_replacement here
            end
            push!(all_ixfaces,(last_i+1,last_i+2,last_i+3))
            if ns==4
                push!(all_ixfaces,(last_i+3,last_i+2,last_i+4))
            end
        end
    end

    function calcxs()
        @inbounds for itet=1:ntet
            node_indices[1]=cellnodes[1,itet]
            node_indices[2]=cellnodes[2,itet]
            node_indices[3]=cellnodes[3,itet]
            node_indices[4]=cellnodes[4,itet]
            planeq[1]=all_planeq[node_indices[1]]
            planeq[2]=all_planeq[node_indices[2]]
            planeq[3]=all_planeq[node_indices[3]]
            planeq[4]=all_planeq[node_indices[4]]
            nxs=tet_x_plane!(ixcoord,ixvalues,coord,node_indices,planeq,func,tol=tol)
            pushtris(nxs,ixcoord,ixvalues)
        end
    end
    
    @inbounds for iplane=1:nplanes
        @views @inbounds map!(inode->plane_equation(planes[iplane],coord[:,inode]),all_planeq,1:nnodes)
        calcxs()
    end
    
    # allocation free (besides push!)
    @inbounds for ilevel=1:nlevels
        @views @inbounds @fastmath map!(inode->(func[inode]-flevels[ilevel]),all_planeq,1:nnodes)
        calcxs()
    end
    
    all_ixcoord, all_ixfaces, all_ixvalues
end


"""
    $(SIGNATURES)

Collect isoline snippets on triangles ready for linesegments!
"""
function marching_triangles(coord,cellnodes,func,levels)
    points=Vector{Point2f}(undef,0)
    function isect(nodes)
        (i1,i2,i3)=(1,2,3)

        f=(func[nodes[1]],func[nodes[2]],func[nodes[3]])

        f[1]  <= f[2]  ?  (i1,i2) = (1,2)   : (i1,i2) = (2,1)
        f[i2] <= f[3]  ?  i3=3              : (i2,i3) = (3,i2)
        f[i1] >  f[i2] ?  (i1,i2) = (i2,i1) : nothing

        (n1,n2,n3)=(nodes[i1],nodes[i2],nodes[i3])
        
        dx31=coord[1,n3]-coord[1,n1]
        dx21=coord[1,n2]-coord[1,n1]
        dx32=coord[1,n3]-coord[1,n2]
        
        dy31=coord[2,n3]-coord[2,n1]
        dy21=coord[2,n2]-coord[2,n1]
        dy32=coord[2,n3]-coord[2,n2]

        df31 = f[i3]!=f[i1] ? 1/(f[i3]-f[i1]) : 0.0
        df21 = f[i2]!=f[i1] ? 1/(f[i2]-f[i1]) : 0.0
        df32 = f[i3]!=f[i2] ? 1/(f[i3]-f[i2]) : 0.0

        for level ∈ levels
            if  (f[i1]<=level) && (level<f[i3]) 
	        α=(level-f[i1])*df31
	        x1=coord[1,n1]+α*dx31
	        y1=coord[2,n1]+α*dy31
                
	        if (level<f[i2])
	            α=(level-f[i1])*df21
	            x2=coord[1,n1]+α*dx21
		    y2=coord[2,n1]+α*dy21
                else
	            α=(level-f[i2])*df32
	            x2=coord[1,n2]+α*dx32
	            y2=coord[2,n2]+α*dy32
                end
                push!(points,Point2f(x1,y1))
                push!(points,Point2f(x2,y2))
            end
        end
    end
    
    for itri=1:size(cellnodes,2)
        @views isect(cellnodes[:,itri])
    end

    points
end    




"""
$(SIGNATURES)


Assume that `points` are nodes of a polyline.
Place `nmarkers` equidistant markers  at the polyline, under
the assumption that the points are transformed via the transformation
matrix M vor visualization.
"""
function markerpoints(points,nmarkers,transform)
    dist(p1,p2)=norm(transform*(p1-p2))
    
    llen=0.0
    for i=2:length(points)
       llen+=dist(points[i],points[i-1])
    end

    mdist=llen/(nmarkers-1)
    
    mpoints=[points[1]]
    
    i=2
    l=0.0
    lnext=l+mdist
    while i<length(points)
        d=dist(points[i],points[i-1])
        while l+d <= lnext && i<length(points)
            i=i+1
            l=l+d
            d=dist(points[i],points[i-1])
        end
        
        while lnext <=l+d &&  length(mpoints)<nmarkers-1
            α=(lnext-l)/d
            push!(mpoints,Point2f(α*points[i]+ (1-α)*points[i-1]))
            lnext=lnext+mdist
        end
    end
    push!(mpoints,points[end])
end

function makeplanes(mmin,mmax,n)
    if isa(n,Number)
        if n==0
            return [Inf]
        end
        p=collect(range(mmin,mmax,length=ceil(n)+2))
        p[2:end-1]
    else
        n
    end
end

function makeplanes(xyzmin,xyzmax,x,y,z)
    planes=Vector{Vector{Float64}}(undef,0)
#    ε=1.0e-1*(xyzmax.-xyzmin)
    
    X=makeplanes(xyzmin[1],xyzmax[1],x)
    Y=makeplanes(xyzmin[2],xyzmax[2],y)
    Z=makeplanes(xyzmin[3],xyzmax[3],z)

    for i=1:length(X)
        x=X[i]
        x>xyzmin[1] && x<xyzmax[1]  && push!(planes,[1,0,0,-x])
    end
    
    for i=1:length(Y)
        y=Y[i]
        y>xyzmin[2] && y<xyzmax[2]  && push!(planes,[0,1,0,-y])
    end
    
    for i=1:length(Z)
        z=Z[i]
        z>xyzmin[3] && z<xyzmax[3]  && push!(planes,[0,0,1,-z])
    end
    planes
end



# Calculate isolevel values and function limits
function isolevels(ctx,func)
    
    crange=ctx[:limits]
    if crange==:auto || crange[1]>crange[2] 
        crange=extrema(func)
    end
    
    if isa(ctx[:levels],Number)
        levels=collect(LinRange(crange[1],crange[2],ctx[:levels]+2))
    else
        levels=ctx[:levels]
    end

    if ctx[:colorbarticks] == :default
        colorbarticks = levels
    elseif isa(ctx[:colorbarticks],Number)
        colorbarticks = collect(crange[1]:(crange[2]-crange[1])/(ctx[:colorbarticks]-1):crange[2])
    else
        colorbarticks = ctx[:colorbarticks]
    end
    
    
    #    map(t->round(t,sigdigits=4),levels),crange,map(t->round(t,sigdigits=4),colorbarticks)
    levels,crange,colorbarticks

end

