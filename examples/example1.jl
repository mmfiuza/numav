# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

using Numav

# create the simulation object with some numerical method
s = Simulation{
    Phenomenon.acoustic,
    NumericalMethod.fem,
    Domain.frequency,
    Dimension.d3
}()

# set element order
set_element_order(s, 2)

# determine simulation frequency range
freq_min = 0
freq_max = 100
set_freq_limits(s, freq_min, freq_max)

# load the mesh
load_mesh(s, "example1.bdf")

# add a volume material (air in this case)
id_air = 1
rho_air = 1.20
c_air = 343
add_volume_material(s, id_air, rho_air, c_air)

# add volume velocity sources
Q(f) = 10/f
add_source(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
add_source(s, surface_id=2, volume_velocity=Q)

# add pressure sources
P(f) = 2/f
add_source(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_source(s, surface_id=3, pressure=P)

# add specific surface impedance
Z(f) = 1f + 2im
add_specific_surface_acoustic_impedance(s, surface_id=4, impedance=Z)

# run the simulation
result = Numav.run(s)
