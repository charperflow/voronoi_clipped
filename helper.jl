#helper.jl
module Helper

using DelaunayTriangulation, Distributions, Random

export calc_areas, voronoi_area_cdf, voronoi_area_pdf, find_homes, get_triangle_areas_pts, scatter_over_hull

distance(x1,x2) = sqrt(sum((x2-x1) .^2))
triangle_area(p, q, r) =  0.5 * (p[1] * q[2] + q[1] * r[2] + r[1] * p[2] - p[1] * r[2] - r[1] * q[2] - q[1] * p[2])

function get_triangle_areas_pts(tri)
    n = length(tri.triangles)
    A = Array{Float64}(undef,n)
    P = []
    idx = 1
    for T in each_triangle(tri)
        i, j, k = indices(T)
        p, q, r = get_point(tri, i, j, k)
        A[idx] = triangle_area(p, q, r)
        push!(P,[p,q,r])
        idx += 1 
    end
    return A, P
end

function calc_areas(vorn)
    n = length(vorn.polygons)
    polygons = vorn.polygons
    pts = vorn.polygon_points
    areas = Array{Float64}(undef,n)
    for i in 1:n
        polygon=polygons[i]
        m = length(polygon)
        area = 0
        numerator = 0
        count = 1
        xs = [pts[polygon[j]][1] for j in 1:m]
        ys = [pts[polygon[j]][2] for j in 1:m] 
        while count < m 
            numerator += (xs[count]*ys[count+1] - ys[count]*xs[count+1])
            count += 1
        end
        numerator += (xs[m]*ys[1] - ys[m]*xs[1])
        area = numerator/2
        areas[i] = area
    end #for
    return areas
end

function scatter_over_hull(tri,n)
    pts = Array{Tuple{Float64,Float64}}(undef,n)
    tri_As, tri_pts = get_triangle_areas_pts(tri)
    tri_As_ps = tri_As ./ sum(tri_As)
    dist_tri_A = Categorical(tri_As_ps)
    sample_T = rand(dist_tri_A,n)
    Ts = tri_pts[sample_T]

    P1s = [x[1] for x in Ts]
    P2s = [x[2] for x in Ts]
    P3s = [x[3] for x in Ts]

    as = [P3s[i] .- P1s[i] for i in 1:n]
    bs = [P2s[i] .- P1s[i] for i in 1:n]

    u1s = rand(Uniform(0,1),n)
    u2s = rand(Uniform(0,1),n)
    u1s = [(u1s[i]+u2s[i]) > 1 ? 1-u1s[i] : u1s[i] for i in 1:n]
    u2s = [(u1s[i]+u2s[i]) > 1 ? 1-u2s[i] : u2s[i] for i in 1:n]

    Ws = [P1s[i] .+ (u1s[i] .* as[i]) .+ (u2s[i] .* bs[i]) for i in 1:n]

    return Ws
end

function voronoi_area_cdf(areas;res=.001)
    N = size(areas)[1]
    steps = Int(floor(1/res))
    X = [res*i for i in 1:steps]
    Y = zeros(steps)
    max_area = maximum(areas)
    total_area = sum(areas)
    for i in 1:steps
        count = 0
        for ii in 1:N
            if (areas[ii]/max_area)<=X[i]
                count = count+areas[ii]/total_area
            end
        end
        Y[i] = count
    end
    return X,Y
end #function

function voronoi_area_pdf(areas;res=.001)
    N = size(areas)[1]
    steps = Int(floor(1/res))
    X = [res*i for i in 1:steps]
    Y = zeros(steps)
    max_area = maximum(areas)
    total_area = sum(areas)
    for i in 1:steps
        count = 0
        lb = (i == 1 ? 0 : X[i-1])
        ub = X[i]
        for ii in 1:N
            if (areas[ii]/max_area >= lb) && (areas[ii]/max_area <ub)
                count = count + areas[ii]/total_area
            end
        end
        Y[i] = count/res
    end
    return X,Y
end #function

function find_homes(pts,gens)
    M = size(pts)[1]
    N = size(gens)[1]
    G = zeros(N)
    ties = []
    for i in 1:M
        min_dist = distance(pts[i],gens[1])
        home = 1
        ties = []
        for j in 2:N
            dist = distance(pts[i],gens[j])
            if dist < min_dist
                home = j
                min_dist = dist
            elseif dist == min_dist
                push!(ties,j)
            end #conditional
        end #j

        #Now we keep track of which generator is the home
        G[home] += 1
        #splitting up any ties that may have occured 
        #should never be more than 3 
        n = size(ties)[1]
        if !isempty(ties)
            for k in 1:n
                G[ties[k]]+=(1/n)
            end
        end
    end# i loop 
    return G
end #function

end