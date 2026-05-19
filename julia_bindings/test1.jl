# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Numav

# create the simulation object with some numerical method
s = Simulation{
    Phenomenon.acoustic,
    NumericalMethod.fem,
    Domain.frequency,
    Dimension.d3,
    ElementOrder.o1
}()

# determine simulation frequency
freq_max = 100
set_maximum_frequency(s, freq_max)

# load the mesh
load_mesh(s, "tests/test1.bdf")

# add a volume materials
rho_air_1(f) = 1.20
c_air_1(f) = 343
add_volume_material(s, 1, rho_air_1, c_air_1)
rho_air_2(f) = 1.20 + 2im
c_air_2(f) = 343 + 3im
add_volume_material(s, 2, rho_air_2, c_air_2)

# add volume velocity sources
Q(f) = 10/f
add_sound_source(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
U(f) = 15/f
add_sound_source(s, physical_group=2, particle_velocity=Q)

# add pressure sources
P(f) = 2f
add_sound_source(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_sound_source(s, physical_group=3, pressure=P)

# add specific surface acoustic impedance
Z(f) = 1f + 2im
add_specific_surface_acoustic_impedance(s, physical_group=4, impedance=Z)

set_result_export_path(s, "result.nmvr")

# run the simulation
result = simulate(s)
