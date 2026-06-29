# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

@test begin

using Numav

s = Simulation{
    Phenomenon.acoustic,
    NumericalMethod.fem,
    Domain.frequency,
    Dimension.d3,
    ElementShape.tetrahedron,
    ElementOrder.o1
}()

set_maximum_frequency(s, 100)
set_frequency_sampling_density(s, FrequencySamplingDensity.quadratic)
load_mesh(s, "test1.bdf")

pqv = "pqv.txt"

add_volume_material(
    s, physical_group=1, density=f->f, sound_speed=f->f
)
add_volume_material(
    s, physical_group=2, density=1, sound_speed=1
)
add_volume_material(
    s, physical_group=3, density=pqv, sound_speed=pqv
)

add_sound_source(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=f->f)
add_sound_source(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=1)
add_sound_source(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=pqv)
add_sound_source(s, physical_group=1, particle_velocity=f->f)
add_sound_source(s, physical_group=2, particle_velocity=1)
add_sound_source(s, physical_group=3, particle_velocity=pqv)
add_sound_source(s, coordinates=[1.0, 1.0, 1.0], pressure=f->f)
add_sound_source(s, coordinates=[1.0, 1.0, 1.0], pressure=1)
add_sound_source(s, coordinates=[1.0, 1.0, 1.0], pressure=pqv)
add_sound_source(s, physical_group=4, pressure=f->f)
add_sound_source(s, physical_group=5, pressure=1)
add_sound_source(s, physical_group=6, pressure=pqv)

add_surface_material(s, physical_group=7, impedance=f->f)
add_surface_material(s, physical_group=8, impedance=1)
add_surface_material(s, physical_group=9, impedance=pqv)

set_result_export_path(s, "test1.h5")

run(s)

return true # Passed if code reaches here

end # @test

