
# air properties
temperature_celsius = 20;
atmospheric_pressure = 101325; # [Pa]

# room dimensions
lx = 5; # [m]
ly = 4; # [m]
lz = 3; # [m]

# define positions
source_position = [0.0, 0.0, 0.0]; # x,y,z [m]
receiver_position = [3.07629, 1.71063, 1.83652]; # x,y,z [m]

# define source volume velocity
volume_velocity_amplitude(freq) = 1 ./ freq

# %% Analytic solution

include("../src/Analytic.jl")
using .Analytic

freq_ana = 1 : 0.1 : 100

complex_pressure_amplitude_ana = simulate_analytic_rectangular_room(
    [lx,ly,lz], source_position, receiver_position, temperature_celsius, 
    atmospheric_pressure, freq_ana, volume_velocity_amplitude);

spl_ana = 20*log10.(abs.(complex_pressure_amplitude_ana)/sqrt(2)/20e-6);

# %% FEM solution

include("../src/FemFreq3d.jl")
using .FemFreq3d

element_order = 1

# frequency array for the simulation
freq_fem = 1 : 0.1 : 100;

# read mesh
mesh = read_bdf("./julia-experiments/examples/mesh.bdf");

# discover the node index for sources and receivers
idx_receiver = get_index_of_closest_node(receiver_position, mesh.nodes)
idx_source = get_index_of_closest_node(source_position, mesh.nodes)

complex_pressure_amplitude_fem = simulate(mesh, idx_source, element_order, temperature_celsius, 
    atmospheric_pressure, freq_fem, volume_velocity_amplitude)

spl_fem = 20*log10.(abs.(complex_pressure_amplitude_fem)/sqrt(2)/20e-6);

# %% Plot frf

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
    framestyle = :box,
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
    vec(freq_fem), vec(spl_fem[idx_receiver,:]),
    label = "numeric"
)

# %% Plot mesh

using GLMakie

function plot_mesh(mesh)

    segments_set = Set{Tuple{Int64, Int64}}()

    for i in 1:size(mesh.elem_vol, 2)
        elem = mesh.elem_vol[:, i]
        corners = elem[1:4] # vertex indices
        midpts = elem[5:10] # mid-edge indices

        # Define edges and their midpoints
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

    fig = GLMakie.Figure()
    ax = GLMakie.Axis3(fig[1, 1], title="Mesh", aspect=:data)

    # Plot all line segments
    GLMakie.linesegments!(ax, segments_points, color=:black, linewidth=1)

    # add midpoints as markers
    midpoint_indices = unique(mesh.elem_vol[5:10, :])
    midpoint_coords = [Point3(mesh.nodes[:, i]) for i in midpoint_indices]
    GLMakie.scatter!(ax, midpoint_coords, color=:red, markersize=5)

    fig
end

plot_mesh(mesh)

