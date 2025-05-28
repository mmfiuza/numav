module  Solver

export fem_solve

using SparseArrays
using LinearAlgebra

function fem_solve(
mesh, idx_source, element_order, rho_air, c_air, freq, volume_velocity_amplitude::Function)

    # define shape functions
    function shape_function(ξ1,ξ2,ξ3)
        ξ4 = 1-ξ1-ξ2-ξ3

        if element_order == 1
            return [
                ξ1,
                ξ2,
                ξ3,
                ξ4,
            ];

        elseif element_order == 2
            return [
                ξ1*(2ξ1-1),
                ξ2*(2ξ2-1),
                ξ3*(2ξ3-1),
                ξ4*(2ξ4-1),
                4ξ1*ξ2,
                4ξ1*ξ3,
                4ξ1*ξ4,
                4ξ2*ξ3,
                4ξ2*ξ4,
                4ξ3*ξ4
            ];
        end
    end
        
    function shape_function_gradient(ξ1,ξ2,ξ3)

        if element_order == 1
            return [
                1 0 0 -1;
                0 1 0 -1;
                0 0 1 -1;
            ];

        elseif element_order == 2
            ξ4 = 1-ξ1-ξ2-ξ3
            return [
                4ξ1-1 0     0     1-4ξ4 4ξ2 4ξ3 4(ξ4-ξ1) 0   -4ξ2     -4ξ3     ;
                0     4ξ2-1 0     1-4ξ4 4ξ1 0   -4ξ1     4ξ3 4(ξ4-ξ2) -4ξ3     ;
                0     0     4ξ3-1 1-4ξ4 0   4ξ1 -4ξ1     4ξ2 -4ξ2     4(ξ4-ξ3) ;
            ];
        end
    end

    if element_order == 1
        nodes_per_element = 4
    elseif element_order == 2
        nodes_per_element = 10
    end

    stiffness_elemental_matrix = 
        [zeros(Float64, nodes_per_element, nodes_per_element) for _ in 1:mesh.num_elem_vol]
    mass_elemental_matrix = 
        [zeros(Float64, nodes_per_element, nodes_per_element) for _ in 1:mesh.num_elem_vol]

    for ne in 1:mesh.num_elem_vol

        # matrix of xyz coordinates of each vertex in the element
        elem_coord_matrix = hcat(
            mesh.nodes[1, mesh.elem_vol[:,ne]], # x
            mesh.nodes[2, mesh.elem_vol[:,ne]], # y
            mesh.nodes[3, mesh.elem_vol[:,ne]]  # z
        );

        # volume of a tetrahedron in the reference domain
        reference_element_volume = 1/6

        # stiffness matrices
        if element_order == 1
            num_points_stiffness = 1
            gauss_points = [[1/4, 1/4, 1/4]]
            weights = [1] * reference_element_volume
        elseif element_order == 2
            num_points_stiffness = 4
            a = (5 - sqrt(5)) / 20;
            b = (5 + 3*sqrt(5)) / 20;
            gauss_points = [[a,a,a], [b,a,a], [a,b,a], [a,a,b]]
            weights = [1/4, 1/4, 1/4, 1/4] * reference_element_volume
        end
        for ng in 1:num_points_stiffness
            jacobian = shape_function_gradient(gauss_points[ng]...)*elem_coord_matrix;
            det_jacobian = det(jacobian)
            inv_jacobian = inv(jacobian)
            b_matrix = inv_jacobian * shape_function_gradient(gauss_points[ng]...)
            btb = transpose(b_matrix) * b_matrix
            stiffness_elemental_matrix[ne] += (1/rho_air) * btb * det_jacobian * weights[ng]
        end
        
        # mass matrices
        if element_order == 1
            num_points_mass = 4
            a = (5 - sqrt(5)) / 20;
            b = (5 + 3*sqrt(5)) / 20;
            gauss_points = [[a,a,a], [b,a,a], [a,b,a], [a,a,b]]
            weights = [1/4, 1/4, 1/4, 1/4] * reference_element_volume
        elseif element_order == 2
            num_points_mass = 15
            a = 1/4
            b = (7 + sqrt(15)) / 34
            c = (7 - sqrt(15)) / 34
            d = (13 - 3*sqrt(15)) / 34
            e = (13 + 3*sqrt(15)) / 34
            f = (5 - sqrt(15)) / 20
            g = (5 + sqrt(15)) / 20
            gauss_points = [
                [a,a,a], [b,b,b], [b,b,d], [b,d,b], [d,b,b], 
                [c,c,c], [c,c,e], [c,e,c], [e,c,c], [f,f,g],
                [f,g,f], [g,f,f], [f,g,g], [g,f,g], [g,g,f]
            ]
            w1 = 8/405
            w2 = (2665 - 14*sqrt(15)) / 226_800
            w3 = (2665 + 14*sqrt(15)) / 226_800
            w4 = 5/567
            weights = [w1,w2,w2,w2,w2,w3,w3,w3,w3,w4,w4,w4,w4,w4,w4]
        end
        for ng in 1:num_points_mass
            jacobian = shape_function_gradient(gauss_points[ng]...)*elem_coord_matrix;
            det_jacobian = det(jacobian)
            nnt = shape_function(gauss_points[ng]...) * 
                transpose(shape_function(gauss_points[ng]...))
            mass_elemental_matrix[ne] += 1/(rho_air*c_air^2) * nnt * det_jacobian * weights[ng]
        end
    end

    # flatten matrices
    stiffness_vals = vcat(vec.(stiffness_elemental_matrix)...)
    mass_vals = vcat(vec.(mass_elemental_matrix)...)

    # convert to complex
    stiffness_vals = Vector{ComplexF64}(stiffness_vals)
    mass_vals = Vector{ComplexF64}(mass_vals)

    i = repeat(1:nodes_per_element, nodes_per_element)
    j = repeat(1:nodes_per_element, inner=nodes_per_element)
    row_sparse_indices = vec(mesh.elem_vol[i,:])
    col_sparse_indices = vec(mesh.elem_vol[j,:])

    stiffness_global_matrix = sparse(
        row_sparse_indices, col_sparse_indices, stiffness_vals, mesh.num_nodes, mesh.num_nodes)
        
    mass_global_matrix = sparse(
        row_sparse_indices, col_sparse_indices, mass_vals, mesh.num_nodes, mesh.num_nodes)

    force_vector = zeros(ComplexF64, mesh.num_nodes,1)

    complex_pressure_amplitude = Matrix{ComplexF64}(undef, mesh.num_nodes, length(freq))

    for nf in axes(freq,1)
        
        global_matrix = stiffness_global_matrix - (2pi*freq[nf])^2*mass_global_matrix
        
        force_vector[idx_source] = -im*2pi*freq[nf]*volume_velocity_amplitude(freq[nf])

        complex_pressure_amplitude[:,nf] = global_matrix \ force_vector

        print(nf, "\n")
    end

    return complex_pressure_amplitude
end

end