using Statistics

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