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
set_frequency_range(s, 40, 100)
set_frequency_sampling_density(s, FrequencySamplingDensity.quadratic)

# load the mesh
load_mesh(s, "tests/test1.bdf")

# add a volume materials
rho1(f) = 1.20
c1(f) = 343
add_volume_material(s, physical_group=1, 
    density="julia-bindings/tests/pqv.txt", sound_speed=c1)
rho2(f) = 1.20 + 2im
c2(f) = 343 + 3im
add_volume_material(s, physical_group=2, density=rho2, sound_speed=2)

# add volume velocity sources
Q(f) = 10/f
add_sound_source(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
U(f) = 15/f
add_sound_source(s, physical_group=2, particle_velocity=2)

# add pressure sources
P(f) = 2f
add_sound_source(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_sound_source(s, physical_group=3, pressure=P)

# add specific surface acoustic impedance
Z(f) = 1f + 2im
add_surface_material(s, physical_group=4, 
    impedance="julia-bindings/tests/pqv.txt")

set_result_export_path(s, "result.h5")

# run the simulation
simulate(s)
