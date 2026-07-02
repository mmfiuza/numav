# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

@test begin

using Numav

redirect_stdout(devnull) do

result_path = "test1.h5"
mesh_path = "test1.bdf"
pqv = "pqv.txt"
func(x) = x

s = create_simulation(
    numerical_method = :fem,
    equation = :helmholtz,
    element_shape = :tetrahedron,
    element_order = :linear
)
@test_throws ErrorException run!(s)
set_frequency!(s, steps=[20, 30, 50, 80, 90, 100])
@test_throws ArgumentError set_frequency!(s, steps=[20, 30], max=100)
load_mesh!(s, mesh_path)
@test_throws ErrorException load_mesh!(s, mesh_path)
@test_throws ErrorException run!(s)
add_volume_material!(s, physical_group=1, density=func, speed_of_sound=f->f)
add_volume_material!(s, physical_group=2, density=func, speed_of_sound=1)
add_volume_material!(s, physical_group=3, density=func, speed_of_sound=pqv)
@test_throws ErrorException add_volume_material!(
    s, physical_group=9, density=1, speed_of_sound=1)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=f->f)
add_sound_source!(s, coordinates=[1, 1, 1], volume_velocity=1)
add_sound_source!(s, coordinates=1:3, volume_velocity=pqv)
add_sound_source!(s, physical_group=1, particle_velocity=f->f)
add_sound_source!(s, physical_group=2, particle_velocity=1)
add_sound_source!(s, physical_group=3, particle_velocity=pqv)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], pressure=f->f)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], pressure=1)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], pressure=pqv)
add_sound_source!(s, physical_group=4, pressure=f->f)
add_sound_source!(s, physical_group=5, pressure=1)
add_sound_source!(s, physical_group=6, pressure=pqv)
@test_throws ErrorException add_sound_source!(
    s, physical_group=100, pressure=1)
@test_throws ArgumentError add_sound_source!(
    s, coordinates=[1.0, 1.0], pressure=1)
add_surface_material!(s, physical_group=7, specific_acoustic_impedance=f->f)
add_surface_material!(s, physical_group=8, specific_acoustic_impedance=1)
add_surface_material!(s, physical_group=9, specific_acoustic_impedance=pqv)
@test_throws ErrorException add_surface_material!(
    s, physical_group=100, specific_acoustic_impedance=1)
@test_throws ErrorException run!(s)
set_result_export_path!(s, result_path)
run!(s)
@test_throws ErrorException run!(s)

s = create_simulation(
    numerical_method = :fem,
    equation = :helmholtz,
    element_shape = :tetrahedron,
    element_order = :quadratic
)
set_frequency!(s, max=100)
load_mesh!(s, mesh_path)
add_volume_material!(s, physical_group=1, density=1, speed_of_sound=f->f)
add_volume_material!(s, physical_group=2, density=1, speed_of_sound=1)
add_volume_material!(s, physical_group=3, density=1, speed_of_sound=pqv)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=func)
set_result_export_path!(s, result_path)
run!(s)

s = create_simulation(
    numerical_method = :fem,
    equation = :helmholtz,
    element_shape = :tetrahedron,
    element_order = :linear
)
set_frequency!(s, min=30, max=100, step_count=1000, sampling_density=:constant)
load_mesh!(s, mesh_path)
add_volume_material!(s, physical_group=1, density=pqv, speed_of_sound=f->f)
add_volume_material!(s, physical_group=2, density=pqv, speed_of_sound=1)
add_volume_material!(s, physical_group=3, density=pqv, speed_of_sound=pqv)
add_sound_source!(s, coordinates=[1.0, 1.0, 1.0], volume_velocity=func)
set_result_export_path!(s, result_path)
run!(s)

end # redirect_stdout

return true # Passed if code reaches here

end # @test

