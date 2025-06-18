# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

using Numav

# create the simulation object with some numerical method
sml = Numav.Simulation{
    Numav.Phenomenon.acoustic,
    Numav.NumericalMethod.fem,
    Numav.Domain.frequency,
    Numav.Dimension.d3
}()

# set element order
Numav.set_element_order(sml, 2)

# determine simulation frequency range
freq_min = 0
freq_max = 100
Numav.set_freq_limits(sml, freq_min, freq_max)

# load the mesh
Numav.load_mesh(sml, "example1.bdf")

# add a volume material (air in this case)
id_air = 1
rho_air = 1.20
c_air = 343
Numav.add_volume_material(sml, id_air, rho_air, c_air)

# add volume velocity sources
source_coordinates_1 = [1.0, 1.5, 2.0]
complex_volume_velocity_amplitude(f) = 1/f
Numav.add_source(sml,
    Numav.TypeOfSource.point, source_coordinates_1,
    Numav.PhysicalQuantity.volume_velocity, complex_volume_velocity_amplitude
)

# id_surface_source_1 = 2
# Numav.add_source(sml,
#     TypeOfSource.surface, id_surface_source_1,
#     PhysicalQuantity.volume_velocity, get_complex_volume_velocity_amplitude
# )

# # add pressure sources
# source_coordinates_2 = [2.0, 2.5, 1.0]
# Numav.add_source(sml,
#     TypeOfSource.point, source_coordinates_2,
#     PhysicalQuantity.pressure, get_complex_pressure_amplitude
# )
# id_surface_source_2 = 3
# Numav.add_source(sml,
#     TypeOfSource.surface, id_surface_source_2,
#     PhysicalQuantity.pressure, get_complex_pressure_amplitude
# )

# # add specific surface impedance
# id_surface_impedance = 4
# Numav.add_specific_surface_acoustic_impedance(sml,
#     id_surface_impedance, get_specific_surface_acoustic_impedance
# )

# # run the simulation
# complex_pressure_amplitude = Numav.run(sml)
