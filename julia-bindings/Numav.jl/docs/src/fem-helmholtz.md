# FEM-Helmholtz

This section explains how Numav can perform acoustic simulations based on the Helmholtz equation.

## Quick start

Here is a complete example of a simulation code:

```julia
using Numav

# Create the simulation object
s = create_simulation(
    numerical_method = Fem,
    equation = Helmholtz,
    element_shape = Tetrahedron,
    element_order = Linear
)

# Set simulated frequencies
set_frequency!(s, max=200) # (Hz)

# Load mesh
load_mesh!(s, "path/to/mesh_file.bdf")

# Define material (frequency-dependent)
rho(f) = 1.20
c(f) = 343
add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)

# Add a surface impedance boundary condition
Z(f) = 1f + 2im
add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)

# Add a monopole volume velocity source at a point
Q(f) = 10/f
add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)

# Add a particle velocity source on a surface
U(f) = 15/f
add_sound_source!(s, physical_group=2, particle_velocity=U)

# Add pressure sources (to point or surface)
P(f) = 2f
add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_sound_source!(s, physical_group=3, pressure=P)

# Set output file and run
set_result_export_path!(s, "result.h5")
run!(s)
```

For more details, look at the [Functions](@ref) section.

## Functions

```@docs; canonical=false
create_simulation
set_frequency!
load_mesh!
add_volume_material!
add_surface_material!
add_sound_source!
set_result_export_path!
run!
```

```@meta
```
