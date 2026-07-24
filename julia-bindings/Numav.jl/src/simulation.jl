# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    Simulation,
    create_simulation,
    set_result_export_path!,
    run!

abstract type Simulation end

abstract type NumericalMethod end
struct FemNumericalMethod <: NumericalMethod end

abstract type Equation end
struct HelmholtzEquation <: Equation end

abstract type ElementShape end
struct TetrahedronElementShape <: ElementShape end

abstract type ElementOrder end
struct LinearElementOrder <: ElementOrder end
struct QuadraticElementOrder <: ElementOrder end

"""
Creates and returns an instance of the `Simulation` type.

| Keyword arguments | Type | Supported options | Description |
|:--|:--|:--|:--|
| `numerical_method` | `Numav.Option` | `Fem` | Numerical method used |
| `equation` | `Numav.Option`  | `Helmholtz` | Differential equation to be solved |
| `element_shape` | `Numav.Option`  | `Tetrahedron` | Geometrical shape of elements |
| `element_order` | `Numav.Option`  | `Linear`, `Quadratic` | Polynomial order of the finite elements |

---
# Examples

> ```julia
> s = create_simulation(
>     numerical_method = Fem,
>     equation = Helmholtz,
>     element_shape = Tetrahedron,
>     element_order = Linear
> )
> ```
"""
function create_simulation(;
    numerical_method::Option,
    equation::Option,
    element_shape::Option,
    element_order::Option
)
    if numerical_method === Fem && equation === Helmholtz
        cpp_args = ()

        if numerical_method === Fem
            cpp_args = (cpp_args..., _cpp_NumericalMethod_fem)
        else
            throw(ArgumentError("Invalid `numerical_method` option"))
        end

        if equation === Helmholtz 
            cpp_args = (cpp_args..., _cpp_Equation_helmholtz)
        else
            throw(ArgumentError("Invalid `equation` option"))
        end

        if element_shape === Tetrahedron
            element_shape_type = TetrahedronElementShape
            cpp_args = (cpp_args..., _cpp_ElementShape_tetrahedron)
        else
            throw(ArgumentError("Invalid `element_shape` option"))
        end

        if element_order === Linear
            element_order_type = LinearElementOrder
            cpp_args = (cpp_args..., _cpp_ElementOrder_linear)
        elseif element_order === Quadratic
            element_order_type = QuadraticElementOrder
            cpp_args = (cpp_args..., _cpp_ElementOrder_quadratic)
        else
            throw(ArgumentError("Invalid `element_order` option"))
        end

        return SimulationFemHelmholtz{
                element_shape_type,
                element_order_type,
            }( _cpp_simulation = _cpp_Simulation{cpp_args...}() )
    end
    throw(ArgumentError("Invalid options"))
end

"""
Specifies the file path where simulation results will be written when calling [`run!`](@ref).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| `file_path` | `String` | Output file path (`.h5`) |

---
# Examples

> ```julia
> set_result_export_path!(s, "result.h5")
> ```

"""
function set_result_export_path!(
    simulation::Simulation,
    path_to_hdf5_file::AbstractString
)
    _cpp_set_result_export_path!(
        simulation._cpp_simulation,
        String(path_to_hdf5_file)
    )
end

"""
Assembles and solves the system of equations for all frequencies in the defined range while writing the results to the specified export path during the solution.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |

---
# Examples

> ```julia
> run!(s)
> ```

---
# Output format

Results are exported as files in the [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/). It contains all the results and passed inputs to setup the simulation. To read the results, it is recommended to use [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) or [HDFView](https://www.hdfgroup.org/download-hdfview/).

With [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/), you can post-process results in Julia like:
```julia
using HDF5
r = h5open("result.h5", "r")
# inspect contents
```
"""
function run!(
    simulation::Simulation
)
    _cpp_run!(simulation._cpp_simulation)
end
