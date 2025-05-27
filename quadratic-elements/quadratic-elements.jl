using SparseArrays
using LinearAlgebra

# air properties
temperature_celsius = 20;
atmospheric_pressure = 101325; # [Pa]
rho_air = atmospheric_pressure / (287*(temperature_celsius + 273.15)); # [kg/m^3]
c_air = sqrt(1.400 * 287.0 * (temperature_celsius + 273.15)); # [m/s]

# room dimensions
lx = 5; # [m]
ly = 4; # [m]
lz = 3; # [m]
volume = lx*ly*lz; # [m3]

# define positions
source_position = [0.0, 0.0, 0.0]; # x,y,z [m]
receiver_position = [3.07629, 1.71063, 1.83652]; # x,y,z [m]

volume_velocity_amplitude = 1 ./ omega;

#%% Analytic solution

# calculate modes
max_mode_idx = 100;
include("./functions/calculate_modes.jl")
psi_source, freq_modes = calculate_modes(source_position, lx, ly, lz, c_air, max_mode_idx);
psi_receiver, _ = calculate_modes(receiver_position, lx, ly, lz, c_air, max_mode_idx);
omega_modes = 2*pi*freq_modes;
k_modes = omega_modes/c_air;

freq_ana = 1 : 0.1 : 100;
omega_ana = 2*pi*freq_ana;
k = omega_ana/c_air;

pressure_amplitude = im*omega*rho_air.*volume_velocity_amplitude/volume .* sum(
    transpose(psi_receiver).*transpose(psi_source) ./ (k.^2 .- transpose(k_modes).^2),
    dims=2
);

spl_ana = 20*log10.(abs.(pressure_amplitude)/sqrt(2)/20e-6);

#%% FEM solution

@enum ElementOrder Linear Quadratic
element_order = Linear

# read mesh
include("./functions/read_bdf.jl")
mesh = read_bdf("./mesh.bdf");

# frequency array for the simulation
freq = 1 : 0.1 : 100;

# eveluate maximum frequency allowed
include("./functions/get_maximum_freq_allowed.jl")
freq_max_allowed = get_maximum_freq_allowed(mesh, c_air, 6)
println("Maximum allowed frequency: $(round(Int, freq_max_allowed)) Hz")

# create new node on each edge middle point for quadratic elements
if element_order == Quadratic
    include("./functions/create_edge_points_for_quadratic_elements!.jl")
    create_edge_points_for_quadratic_elements!(mesh)
end

# discover the node index for sources and receivers
include("./functions/get_index_of_closest_node.jl")
idx_source = get_index_of_closest_node(source_position, mesh.nodes)
idx_receiver = get_index_of_closest_node(receiver_position, mesh.nodes)

# define shape functions
if element_order == Linear
    
    function shape_function(ξ1,ξ2,ξ3)
        ξ4 = 1-ξ1-ξ2-ξ3
        return [
            ξ1,
            ξ2,
            ξ3,
            ξ4,
        ];
    end
    
    function shape_function_gradient(ξ1,ξ2,ξ3)
        return [
            1 0 0 -1;
            0 1 0 -1;
            0 0 1 -1;
        ];
    end
    Quadratic
elseif element_order == Quadratic
    
    function shape_function(ξ1,ξ2,ξ3)
        ξ4 = 1-ξ1-ξ2-ξ3
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
                
    function shape_function_gradient(ξ1,ξ2,ξ3)  
        ξ4 = 1-ξ1-ξ2-ξ3
        return [
            4ξ1-1 0     0     1-4ξ4 4ξ2 4ξ3 4(ξ4-ξ1) 0   -4ξ2     -4ξ3     ;
            0     4ξ2-1 0     1-4ξ4 4ξ1 0   -4ξ1     4ξ3 4(ξ4-ξ2) -4ξ3     ;
            0     0     4ξ3-1 1-4ξ4 0   4ξ1 -4ξ1     4ξ2 -4ξ2     4(ξ4-ξ3) ;
        ];
    end
    
end

if element_order == Linear
    nodes_per_element = 4
elseif element_order == Quadratic
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
    if element_order == Linear
        num_points_stiffness = 1
        gauss_points = [[1/4, 1/4, 1/4]]
        weights = [1] * reference_element_volume
    elseif element_order == Quadratic
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
    if element_order == Linear
        num_points_mass = 4
        a = (5 - sqrt(5)) / 20;
        b = (5 + 3*sqrt(5)) / 20;
        gauss_points = [[a,a,a], [b,a,a], [a,b,a], [a,a,b]]
        weights = [1/4, 1/4, 1/4, 1/4] * reference_element_volume
    elseif element_order == Quadratic
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
        nnt = shape_function(gauss_points[ng]...) * transpose(shape_function(gauss_points[ng]...))
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
    
    local global_matrix = stiffness_global_matrix - (2pi*freq[nf])^2*mass_global_matrix
    
    force_vector[idx_source] = -im*2pi*freq[nf]*volume_velocity_amplitude[nf]

    complex_pressure_amplitude[:,nf] = global_matrix \ force_vector

    print(nf, "\n")
end

spl_fem = 20*log10.(abs.(complex_pressure_amplitude)/sqrt(2)/20e-6);


#%% plot

# using CairoMakie

# fig = CairoMakie.Figure(size = (600, 400));
# ax = CairoMakie.Axis(fig[1,1],
#     xscale = log10,
#     xticks = (10: 10 : 100, string.(10: 10 : 100)), 
#     xlabel = "Frequency (Hz)",
#     ylabel = "SPL (dB)",
#     title = "Frequency response",
#     titlesize = 14,
#     xticksize = 10,
#     yticksize = 10,
# );
# CairoMakie.xlims!(ax, 25, 100);
# CairoMakie.ylims!(ax, 0, 200);
# CairoMakie.lines!(ax, vec(freq_ana), vec(spl));
# CairoMakie.lines!(ax, vec(freq), vec(spl_fem[idx_receiver,:]));
# fig
# CairoMakie.save("plot.pdf", fig)

using Plots

Plots.plot(
    xaxis = :log10,
    xlim = [25, 100],
    ylim = [0, 200],
    xticks = (10: 10 : 100, string.(10: 10 : 100)),
    xlabel = "Frequency (Hz)",
    ylabel = "SPL (dB)",
    title = "Frequency response",
    grid = true,
    gridalpha = 0.25,
    legendposition = (0.11, 0.9),
    # annotation = (39.5, 127, Plots.text("text here", 10, :left)),
    framestyle = :box,
    # label = "something",
    linestyle = :solid,
    titlefontsize = 14,
    legendfontsize = 10,
    xtickfontsize = 10,
    ytickfontsize = 10,
    xguidefontsize = 12,
    yguidefontsize = 12,
    size = (600, 400),
)
Plots.plot!(
    vec(freq_ana), vec(spl_ana),
    label = "analytic"
)
Plots.plot!(
    vec(freq), vec(spl_fem[idx_receiver,:]),
    label = "numeric"
)

# savefig("plot.pdf");

#%% plot mesh

using Makie, GLMakie

segments_set = Set{Tuple{Int64, Int64}}()

for i in 1:size(mesh.elem_vol, 2)
    elem = mesh.elem_vol[:, i]
    corners = elem[1:4]   # Vertex indices
    midpts = elem[5:10]   # Mid-edge indices

    # Define edges and their midpoints (order: [v1-v2, v1-v3, v1-v4, v2-v3, v2-v4, v3-v4])
    edges = [
        (corners[1], corners[2], midpts[1]),
        (corners[1], corners[3], midpts[2]),
        (corners[1], corners[4], midpts[3]),
        (corners[2], corners[3], midpts[4]),
        (corners[2], corners[4], midpts[5]),
        (corners[3], corners[4], midpts[6])
    ]

    for (a, b, m) in edges
        push!(segments_set, (min(a, m), max(a, m)))  # Segment: corner → midpoint
        push!(segments_set, (min(m, b), max(m, b)))  # Segment: midpoint → corner
    end
end

segments = collect(segments_set)

segments_points = Vector{Point{3,Float64}}(undef, 2*length(segments))
count = 1
for (i, j) in segments
    segments_points[count] = Point3(mesh.nodes[1, i], mesh.nodes[2, i], mesh.nodes[3, i])
    count += 1
    segments_points[count] = Point3(mesh.nodes[1, j], mesh.nodes[2, j], mesh.nodes[3, j])
    count += 1
end

# segments_points = [Point3(mesh.nodes[1, i], mesh.nodes[2, i], mesh.nodes[3, i]) 
#     for (i, j) in segments for _ in 1:2]

fig = Figure()
ax = Axis3(fig[1, 1], title="Second-Order Tetrahedral Mesh", aspect=:data)

# Plot all line segments
linesegments!(ax, segments_points, color=:black, linewidth=1)

# add midpoints as markers
midpoint_indices = unique(mesh.elem_vol[5:10, :])
midpoint_coords = [Point3(mesh.nodes[:, i]) for i in midpoint_indices]
scatter!(ax, midpoint_coords, color=:red, markersize=5)

fig 