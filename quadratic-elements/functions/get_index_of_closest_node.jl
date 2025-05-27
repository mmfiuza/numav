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