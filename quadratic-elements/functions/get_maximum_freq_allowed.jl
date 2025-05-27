
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