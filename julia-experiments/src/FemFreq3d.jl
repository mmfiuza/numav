module FemFreq3d

export simulate, get_index_of_closest_node, read_bdf

include("MeshProcessing.jl")
include("Solver.jl")

using .MeshProcessing
using .Solver

function simulate(mesh, idx_source, element_order, temperature_celsius, atmospheric_pressure, 
freq, volume_velocity_amplitude::Function)

    rho_air = atmospheric_pressure / (287*(temperature_celsius + 273.15)); # [kg/m^3]
    c_air = sqrt(1.400 * 287.0 * (temperature_celsius + 273.15)); # [m/s]

    # eveluate maximum frequency allowed
    freq_max_allowed = get_maximum_freq_allowed(mesh, c_air, 6)
    println("Maximum allowed frequency: $(round(Int, freq_max_allowed)) Hz")

    # create new node on each edge middle point for quadratic elements
    if element_order == 2
        create_edge_points_for_quadratic_elements!(mesh)
    end

    complex_pressure_amplitude = fem_solve(
        mesh, idx_source, element_order, rho_air, c_air, freq, volume_velocity_amplitude)

    return complex_pressure_amplitude
end

end