# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# export functions

export load_mesh!

function _check_if_simulation_has_mesh(s::Simulation)
    if ( 
        !hasproperty(s, :numerical_method) ||
        !(s.numerical_method isa FemNumericalMethod)
    )
        _throw_simulation_not_applicable()
    end
    return
end

"""
Loads the mesh file. It defines the geometry of the domain and must contain the physical groups that are referenced when assigning boundary conditions. To generate meshes, a good open-source option is [Gmsh](https://gmsh.info/).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| `path_to_mesh` | `String` | Path to the mesh file (`.bdf`) |

---
# Examples

> ```julia
> load_mesh!(s, "some_mesh_file.bdf")
> ```

!!! warning
    Currently, only meshes in the **small field** BDF format are supported.

!!! warning
    The coordinates in the mesh file must be given in **meters**.
"""
function load_mesh!(
    simulation::Simulation,
    path_to_mesh::AbstractString
)
    _check_if_simulation_has_mesh(simulation)
    _cpp_load_mesh!(simulation._cpp_simulation, String(path_to_mesh))
end
