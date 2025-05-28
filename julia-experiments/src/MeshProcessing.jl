module MeshProcessing

export get_index_of_closest_node, create_edge_points_for_quadratic_elements!, 
get_maximum_freq_allowed, read_bdf

using Statistics
using Printf

function get_index_of_closest_node(point, nodes)
    min_dist_sq = Inf
    index = 1
    for i in axes(nodes, 2)
        dx = nodes[1,i] - point[1]
        dy = nodes[2,i] - point[2]
        dz = nodes[3,i] - point[3]
        dist_sq = dx^2 + dy^2 + dz^2
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq
            index = i
        end
    end
    return index
end

function create_edge_points_for_quadratic_elements!(mesh)

    existing_edge_nodes = Set{Tuple{Int,Int}}()
    edge_nodes_indices = Dict{Tuple{Int,Int}, Int}()
    new_elem_vol = Matrix{Int}(undef, 10, mesh.num_elem_vol)

    for e in axes(mesh.elem_vol, 2)

        # copy the 4 vertex nodes to the new matrix
        new_elem_vol[1,e] = mesh.elem_vol[1,e]
        new_elem_vol[2,e] = mesh.elem_vol[2,e]
        new_elem_vol[3,e] = mesh.elem_vol[3,e]
        new_elem_vol[4,e] = mesh.elem_vol[4,e]
        
        # add the 6 edge nodes
        vertex_pairs = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]
        for i = axes(vertex_pairs,1)
            cur_vertexes_idxs = Tuple(sort([
                mesh.elem_vol[vertex_pairs[i][1],e],
                mesh.elem_vol[vertex_pairs[i][2],e]
            ]))
            if !in(cur_vertexes_idxs, existing_edge_nodes) # create new node

                # calculate coordinates
                x = mean([mesh.nodes[1,cur_vertexes_idxs[1]], mesh.nodes[1,cur_vertexes_idxs[2]]])
                y = mean([mesh.nodes[2,cur_vertexes_idxs[1]], mesh.nodes[2,cur_vertexes_idxs[2]]])
                z = mean([mesh.nodes[3,cur_vertexes_idxs[1]], mesh.nodes[3,cur_vertexes_idxs[2]]])
                mesh.nodes = hcat(mesh.nodes, [x,y,z]);

                # put new node in the set of existing edge nodes
                push!(existing_edge_nodes, cur_vertexes_idxs)

                # increase node count
                mesh.num_nodes += 1

                # put new node in the indices dictionary
                edge_nodes_indices[cur_vertexes_idxs] = mesh.num_nodes
                
                # put the new node in the matrix of volume elements
                new_elem_vol[4+i,e] = mesh.num_nodes

            else # node already created

                # put the already existing node in the matrix of volume elements
                new_elem_vol[4+i,e] = edge_nodes_indices[cur_vertexes_idxs]

            end
        end
    end
    mesh.elem_vol = new_elem_vol;
end

function get_phase_vel(complex_speed_of_sound)
    return 1 ./ (real(1 ./ complex_speed_of_sound));
end

function get_distance(a,b)
    return sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2);
end

function get_maximum_freq_allowed(mesh, c, nodes_per_lambda)
    freq_allowed = Vector{Float64}(undef, 6*mesh.num_elem_vol);
    i = 0;
    for ne = 1:mesh.num_elem_vol
        for nv = [[4,1], [4,2], [4,3], [1,2], [1,3], [2,3]]
            i += 1
            dist = get_distance(
                mesh.nodes[:,mesh.elem_vol[nv[1],ne]],
                mesh.nodes[:,mesh.elem_vol[nv[2],ne]]
            )
            freq_allowed[i] = get_phase_vel(c) / (nodes_per_lambda*dist)
        end
    end
    sort!(freq_allowed);
    outlier_partition = 0.001;
    last_outlier_index = Int(round(outlier_partition*6*mesh.num_elem_vol));
    outliers = freq_allowed[1:last_outlier_index];
    freq_max_allowed = outliers[last_outlier_index];
    return freq_max_allowed
end

mutable struct Mesh
    nodes::Matrix{Float64}
    elem_surf::Matrix{Int}
    elem_vol::Matrix{Int}
    id_surf::Vector{Int}
    id_vol::Vector{Int}
    num_nodes::Int
    num_elem_surf::Int
    num_elem_vol::Int
end

function read_bdf(infilename)

    function replace_pm(S_in)
        S = S_in
        positions = vcat(findall(==('+'), S), findall(==('-'), S))
        sort!(positions, rev=true)
        
        for pos in positions
            if pos != 1 && S[prevind(S, pos)] in '0':'9'
                S = S[1:prevind(S, pos)] * "e" * S[pos:end]
            end
        end
        return S
    end
    
    lines = readlines(infilename)
    
    # First pass: Counting
    num_nodes = 0
    num_ele_surf = 0
    num_ele_vol = 0
    for line in lines
        if(startswith(line, "GRID")) num_nodes+=1 end
        if(startswith(line, "CTRIA3")) num_ele_surf+=1 end
        if(startswith(line, "CTETRA")) num_ele_vol+=1 end
    end
    
    # Initialize arrays
    nodes = Matrix{Float64}(undef, 3, num_nodes)
    elements_surf = Matrix{Int}(undef, 3, num_ele_surf)
    elements_vol = Matrix{Int}(undef, 4, num_ele_vol)
    id_surf = Vector{Int}(undef, num_ele_surf)
    id_vol = Vector{Int}(undef, num_ele_vol)
    
    # Second pass: Parse data
    idx_cur_node = 0
    idx_cur_ele_surf = 0
    idx_cur_ele_vol = 0

    for line in lines
        if startswith(line, "GRID")
            idx_cur_node += 1
            if startswith(line, "GRID*")
                # Long format parsing
                onn = parse(Int, line[9:24])
                sx = replace_pm(line[41:56])
                sy = replace_pm(line[57:72])
                sz = replace_pm(line[73:88])
            else
                # Short format parsing
                onn = parse(Int, line[9:16])
                sx = replace_pm(line[25:32])
                sy = replace_pm(line[33:40])
                sz = replace_pm(line[41:48])
            end
            
            nodes[:,idx_cur_node] .= [parse(Float64,sx), parse(Float64,sy), parse(Float64,sz)]
        
        elseif startswith(line, "CTRIA3")
            idx_cur_ele_surf += 1
            data = parse.(Int, split(line[17:end]))
            id_surf[idx_cur_ele_surf] = data[1]
            elements_surf[:,idx_cur_ele_surf] = [data[2], data[3], data[4]]
        
        elseif startswith(line, "CTETRA")
            idx_cur_ele_vol += 1
            data = parse.(Int, split(line[17:end]))
            id_vol[idx_cur_ele_vol] = data[1]
            elements_vol[:,idx_cur_ele_vol] = [data[2], data[3], data[4], data[5]]
        
        end
    end
    
    return Mesh(
        nodes,
        elements_surf,
        elements_vol,
        id_surf,
        id_vol,
        num_nodes,
        num_ele_surf,
        num_ele_vol
    )
end

end