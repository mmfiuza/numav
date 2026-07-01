# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

@test begin

using Numav

pqv = "pqv.txt"
func(x) = x

s1 = create_simulation(
    numerical_method = :fem,
    equation = :helmholtz,
    element_shape = :tetrahedron,
    element_order = :linear
)
set_maximum_frequency!(s1, 100)
set_frequency_sampling_density!(s1, FrequencySamplingDensity.quadratic)
load_mesh!(s1, "test1.bdf")
add_volume_material!(s1, physical_group=1, density=f->f, speed_of_sound=f->f)
add_volume_material!(s1, physical_group=2, density=1, speed_of_sound=1)
add_volume_material!(s1, physical_group=3, density=pqv, speed_of_sound=pqv)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], volume_velocity=f->f)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], volume_velocity=1)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], volume_velocity=pqv)
add_sound_source!(s1, physical_group=1, particle_velocity=f->f)
add_sound_source!(s1, physical_group=2, particle_velocity=1)
add_sound_source!(s1, physical_group=3, particle_velocity=pqv)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], pressure=f->f)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], pressure=1)
add_sound_source!(s1, coordinates=[1.0, 1.0, 1.0], pressure=pqv)
add_sound_source!(s1, physical_group=4, pressure=f->f)
add_sound_source!(s1, physical_group=5, pressure=1)
add_sound_source!(s1, physical_group=6, pressure=pqv)
add_surface_material!(s1, physical_group=7, impedance=f->f)
add_surface_material!(s1, physical_group=8, impedance=1)
add_surface_material!(s1, physical_group=9, impedance=pqv)
set_result_export_path!(s1, "test1.h5")
run!(s1)

s2 = create_simulation(
    numerical_method = :fem,
    equation = :helmholtz,
    element_shape = :tetrahedron,
    element_order = :linear
)
set_maximum_frequency!(s2, 100)
set_frequency_sampling_density!(s2, FrequencySamplingDensity.quadratic)
load_mesh!(s2, "test1.bdf")
add_volume_material!(s2, physical_group=1, density=func, speed_of_sound=func)
add_volume_material!(s2, physical_group=2, density=func, speed_of_sound=1)
add_volume_material!(s2, physical_group=3, density=func, speed_of_sound=pqv)
add_sound_source!(s2, coordinates=[1, 1, 1], volume_velocity=func)
set_result_export_path!(s2, "test1.h5")
run!(s2)

return true # Passed if code reaches here

end # @test

