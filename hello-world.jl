



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
receiver_position = [1.4782, 1.81717, 1.30383]; # x,y,z [m]

# calculate modes
max_mode_idx = 10;
include("./functions/calculate_modes.jl")
psi_source, freq_modes = calculate_modes(source_position, lx, ly, lz, c_air, max_mode_idx);
psi_receiver, _ = calculate_modes(receiver_position, lx, ly, lz, c_air, max_mode_idx);
omega_modes = 2*pi*freq_modes;
k_modes = omega_modes/c_air;

freq_ana = 1 : 0.1 : 100;
omega = 2*pi*freq_ana;
k = omega/c_air;

volume_velocity_amplitude = 1 ./ omega;

pressure_amplitude = im*omega*rho_air.*volume_velocity_amplitude/volume .* sum(
    transpose(psi_receiver).*transpose(psi_source) ./ (k.^2 .- transpose(k_modes).^2),
    dims=2
);

spl = 20*log10.(abs.(pressure_amplitude)/sqrt(2)/20e-6);



nodes_per_element = 4

# read mesh
include("./functions/read_bdf.jl")
mesh = read_bdf("./mesh.bdf");

# frequency array for the simulation
freq = 1 : 0.1 : 100;

# eveluate maximum frequency allowed
include("./functions/get_maximum_freq_allowed.jl")
freq_max_allowed = get_maximum_freq_allowed(mesh, c_air, 6)
println("Maximum allowed frequency: $(round(Int, freq_max_allowed)) Hz")

# # create new node on each edge middle point
# include("./functions/create_edge_points_for_quadratic_elements!.jl")
# create_edge_points_for_quadratic_elements!(mesh)

# discover the node index for sources and receivers
include("./functions/get_index_of_closest_node.jl")
idx_source = get_index_of_closest_node(source_position, mesh.nodes)
idx_receiver = get_index_of_closest_node(receiver_position, mesh.nodes)

function shape_function(xi1,xi2,xi3)
    xi4 = 1-xi1-xi2-xi3
    return [
        xi4,
        xi1,
        xi2,
        xi3,
    ];
end

function shape_function_gradient(xi1,xi2,xi3)
    return [
        -1 1 0 0 ;
        -1 0 1 0 ;
        -1 0 0 1 ;
    ];
end

# function shape_function(xi1,xi2,xi3)
#     xi4 = 1-xi1-xi2-xi3
#     return [
#         xi4*(2xi4-1),
#         xi1*(2xi1-1),
#         xi2*(2xi2-1),
#         xi3*(2xi3-1),
#         4xi4*xi1,
#         4xi4*xi2,
#         4xi4*xi3,
#         4xi1*xi2,
#         4xi1*xi3,
#         4xi2*xi3
#     ];
# end

# function shape_function_gradient(xi1,xi2,xi3)
#     xi4 = 1-xi1-xi2-xi3
#     return [
#         1-4xi4 4xi1-1 0      0      4(xi4-xi1) -4xi2      -4xi3      4xi2 4xi3 0    ;
#         1-4xi4 0      4xi2-1 0      -4xi1      4(xi4-xi2) -4xi3      4xi1 0    4xi3 ;
#         1-4xi4 0      0      4xi3-1 -4xi1      -4xi2      4(xi4-xi3) 0    4xi1 4xi2 ;
#     ];
# end

using SparseArrays
using LinearAlgebra
mass_elemental_matrix = 
    [Matrix{Float64}(undef, nodes_per_element, nodes_per_element) for _ in 1:mesh.num_elem_vol]
stiffness_elemental_matrix = [
    Matrix{Float64}(undef, nodes_per_element, nodes_per_element) for _ in 1:mesh.num_elem_vol]

for ne in 1:mesh.num_elem_vol

    elem_coord_matrix = hcat(
        mesh.nodes[1, mesh.elem_vol[:,ne]],
        mesh.nodes[2, mesh.elem_vol[:,ne]],
        mesh.nodes[3, mesh.elem_vol[:,ne]]
    );

    jacobian_matrix = shape_function_gradient(0.25,0.25,0.25)*elem_coord_matrix;
    det_jacobian_matrix = det(jacobian_matrix)
    inv_jacobian_matrix = inv(jacobian_matrix)

    b_matrix = inv_jacobian_matrix * shape_function_gradient(0.25,0.25,0.25)
    btb_matrix = transpose(b_matrix) * b_matrix
    btb_det_jacobian = btb_matrix * det_jacobian_matrix / 6
    stiffness_elemental_matrix[ne] = btb_det_jacobian ./ rho_air

    nnt_matrix = shape_function(0.25,0.25,0.25) * transpose(shape_function(0.25,0.25,0.25))
    nnt_det_jacobian = nnt_matrix * det_jacobian_matrix / 6
    mass_elemental_matrix[ne] = nnt_det_jacobian ./ (rho_air * c_air^2)

end

# flatten
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

using CairoMakie

fig = CairoMakie.Figure(size = (600, 400));
ax = CairoMakie.Axis(fig[1,1],
    xscale = log10,
    xticks = (10: 10 : 100, string.(10: 10 : 100)), 
    xlabel = "Frequência (Hz)",
    ylabel = "Nível de pressão sonora (dB)",
    title = "Função de resposta em frequência",
    titlesize = 14,
    xticksize = 10,
    yticksize = 10,
);
CairoMakie.xlims!(ax, 25, 100);
CairoMakie.ylims!(ax, 0, 600);
CairoMakie.lines!(ax, vec(freq_ana), vec(spl));
CairoMakie.lines!(ax, vec(freq), vec(spl_fem[idx_receiver,:]));
fig
# CairoMakie.save("plot.pdf", fig)

# using Plots

# Plots.plot(
#     vec(freq), vec(spl),
#     xaxis = :log10,
#     xlim = [25, 100],
#     ylim = [0, 200],
#     xticks = (10: 10 : 100, string.(10: 10 : 100)), 
#     xlabel = "Frequência (Hz)",
#     ylabel = "Nível de pressão sonora (dB)",
#     title = "Função de resposta em frequência",
#     grid = true,
#     gridalpha = 0.25,
#     legendposition = (0.11, 0.9),
#     annotation = (39.5, 127, Plots.text("Segundo\n modo", 10, :left)),
#     framestyle = :box,
#     label = "Alguma coisa",
#     linestyle = :solid,
#     titlefontsize = 14,
#     legendfontsize = 10,
#     xtickfontsize = 10,
#     ytickfontsize = 10,
#     xguidefontsize = 12,
#     yguidefontsize = 12,
#     size = (600, 400),
# )
# savefig("plot.pdf");