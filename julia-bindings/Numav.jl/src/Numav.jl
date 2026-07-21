# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

# export structs
export Simulation
    
# export options
export
    Fem,
    Helmholtz,
    Tetrahedron,
    Constant,
    Linear,
    Quadratic

# export functions
export
    create_simulation,
    set_frequency!,
    load_mesh!,
    add_volume_material!,
    add_surface_material!,
    add_sound_source!,
    set_result_export_path!,
    run!

# include dependencies
using numav_julia_jll
using CxxWrap
using DelimitedFiles
using Interpolations
using Base.Threads

# wrap the C++ part
@wrapmodule(() -> libnumav_julia)
function __init__()
    @initcxx
    # Set env variables to enable MKL in SDL mode to configure at runtime
    ENV["MKL_INTERFACE_LAYER"] = "ILP64"
    ENV["MKL_THREADING_LAYER"] = "INTEL"
end

# global storage for user-defined functions
const _user_functions::Vector{Ref{Function}} = [ ]
global _next_index::UInt = 1
_lock = ReentrantLock()
function _cmplx_split_and_store(f::Function)
    idx = lock(_lock) do
        global _next_index
        idx = _next_index
        _next_index += 1
        push!(_user_functions, f)
        return idx
    end
    name_real = Symbol("_f", idx, "_r")
    name_imag = Symbol("_f", idx, "_i")
    @eval begin
        ($name_real)(x::Float64) = Float64(real(_user_functions[$idx][](x)))
        ($name_imag)(x::Float64) = Float64(imag(_user_functions[$idx][](x)))
    end
    f_real = @eval CxxWrap.@safe_cfunction($(name_real), Float64, (Float64,))
    f_imag = @eval CxxWrap.@safe_cfunction($(name_imag), Float64, (Float64,))
    return (f_real, f_imag)
end

abstract type Option end

struct FemType <: Option end
const Fem = FemType()

struct HelmholtzType <: Option end
const Helmholtz = HelmholtzType()

struct TetrahedronType <: Option end
const Tetrahedron = TetrahedronType()

struct ConstantType <: Option end
const Constant = ConstantType()

struct LinearType <: Option end
const Linear = LinearType()

struct QuadraticType <: Option end
const Quadratic = QuadraticType()

abstract type Members end

struct Simulation{option_1, option_2, option_3, option_4}
    _m::Members
end

struct MembersFemHelmholtz{S, O} <: Members
    _cpp_simulation::_cpp_Simulation
end

# type aliases
const Fdpq = Union{Function, Number, AbstractString}

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
        cpp_args = (cpp_args..., _cpp_ElementShape_tetrahedron)
    else
        throw(ArgumentError("Invalid `element_shape` option"))
    end

    if element_order === Linear
        cpp_args = (cpp_args..., _cpp_ElementOrder_linear)
    elseif element_order === Quadratic
        cpp_args = (cpp_args..., _cpp_ElementOrder_quadratic)
    else
        throw(ArgumentError("Invalid `element_order` option"))
    end

    if numerical_method === Fem && equation === Helmholtz
        return Simulation{
                numerical_method,
                equation,
                element_shape,
                element_order
            }(
                MembersFemHelmholtz{element_shape, element_order}(
                    _cpp_Simulation{cpp_args...}()
                )
            )
    end
    throw(ErrorException("Invalid options"))
end

function _cubic_range(start::Real, stop::Real, length::Integer)
    @assert length != 0 && length != 1
    x = range(0, 1, length=length)
    return @. cbrt( start^3 + (stop^3 - start^3) * x );
end

"""
Sets the frequecies to perform the simulations with two mutually exclusive ways to use it:
- **automatic mode**: pass `max` (and optionally `min`, `length`, `sampling_density` and `step`) and Numav will automatically determine the frequency vector;
- **manual mode**: pass `vector` directly, giving full control over which frequencies are evaluated.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `max` | `Real` | Upper frequency limit (Hz). Required unless `steps` is provided. Mutually exclusive with `steps`. |
| `min` | `Real` | Lower bound of the frequency range (Hz). Only usable together with `max`. Defaults to `0` if omitted. |
| `length` | `Integer` | Number of frequency steps to compute within the `min`/`max` range. Only usable together with `max`. Defaults to `4096` if omitted. |
| `sampling_density` | `Numav.Option` | Sampling strategy to use within the range. Only usable together with `max`. Defaults to `Quadratic` if omitted. |
| `step` | `Real` | Difference in Hertz between each frequency step, supposing equally spaced steps. Only usable together with `max`. Cannot be used together with `length` and `sampling_density` |
| `vector` | `Vector{Real}` | List of frequencies in Hertz to solve at. Cannot be used together with `max`, `min`, `length`, `sampling_density` and `step`. |

---
# `sampling_density` options

| Mode | Description |
|:--|:--|
| `Quadratic` (_default_) | Frequency discretization density grows quadratically with frequency. It tends to follow the modal density of rooms, being the recommended option for good accuracy and computation time balance. |
| `Constant` | Frequency steps are evenly spaced (uniform spacing). Suitable for broadband analyses where equal resolution at all frequencies is desired. |
 

---
# Examples

> Solve up to a maximum frequency, letting Numav automatically pick the steps from 0 Hz:
> ```julia
> set_frequency!(s, max=200) # solve from 0 Hz to 200 Hz
> ```
 
> Solve within an explicit range, by adding `min`:
> ```julia
> set_frequency!(s, min=20, max=200) # solve from 20 Hz to 200 Hz
> ```
 
> Control the number of frequency steps with `length` (defaults to `4096` when not specified):
> ```julia
> set_frequency!(s, max=200, length=500) # evaluate at 500 frequency points
> ```
> A higher step count gives finer frequency resolution at the cost of longer computation time.

> Control how steps are distributed with `sampling_density` (defaults to `Quadratic` when not specified):
> ```julia
> set_frequency!(s, max=200, sampling_density=Constant)
> ```

> Manually specify the exact frequencies with `vector`, for example to match measurement points or concentrate resolution around a resonance:
> ```julia
> set_frequency!(s, vector=[10, 40, 60, 100]) # Solve only at these frequencies
> ```
"""
function set_frequency!( 
    simulation::Simulation{Fem, Helmholtz};
    max::Union{Real, Nothing} = nothing,
    min::Union{Real, Nothing} = nothing,
    length::Union{Integer, Nothing} = nothing,
    sampling_density::Union{Option, Nothing} = nothing,
    step::Union{Real, Nothing} = nothing,
    vector::Union{AbstractVector{<:Real}, Nothing} = nothing
)
    if isnothing(max) && isnothing(vector)
        throw(ArgumentError("Neither `max` nor `vector` passed"))
    end
    if !isnothing(vector)
        if (
            !isnothing(max) ||
            !isnothing(min) ||
            !isnothing(length) ||
            !isnothing(sampling_density) ||
            !isnothing(step)
        )
            throw(ArgumentError(
                "`vector` and other parameter(s) passed simultaneously"
            ))
        end
    end
    if !isnothing(step) && !isnothing(length)
        throw(ArgumentError(
            "`step` and `length` passed simultaneously"
        ))
    end
    if !isnothing(step) && !isnothing(sampling_density)
        throw(ArgumentError(
            "`step` and `sampling_density` passed simultaneously"
        ))
    end
    if !isnothing(max)
        min = something(min, 0)
        if !isnothing(step)
            vector = range(min, max, step=step)
        else
            length = something(length, 4096)
            sampling_density = something(sampling_density, Quadratic)
            if sampling_density === Constant
                vector = range(min, max, length=length)
            elseif sampling_density === Quadratic
                vector = _cubic_range(min, max, length)
            else
                throw(ArgumentError("Invalid `sampling_density` option"))
            end
        end
    end
    _cpp_set_frequency_vector!(
        simulation._m._cpp_simulation,
        Float64.(Vector(vector))
    )
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
    simulation::Simulation{Fem},
    path_to_mesh::AbstractString
)
    _cpp_load_mesh!(simulation._m._cpp_simulation, String(path_to_mesh))
end

function _read_pqv_table(filename::AbstractString)
    data = readdlm(filename, ',')
    freq = Float64.(data[:, 1])
    pqv = Complex.(Float64.(data[:, 2]), Float64.(data[:, 3]))
    return (freq, pqv)
end

function _pqv_to_function(pqv::Fdpq)::Function
    if pqv isa Number
        return (_ -> pqv)
    elseif pqv isa AbstractString
        interp = linear_interpolation(_read_pqv_table(pqv)...) 
        return (x -> interp(x))
    elseif pqv isa Function
        return pqv
    end
    throw(ArgumentError("Invalid type of `pqv`"))
end

"""
Assigns acoustic material properties to a volumetric region of the mesh, identified by its physical group tag.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID (or vector of IDs) from the mesh |
| `density` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Density in kg/m³ |
| `speed_of_sound` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Speed of sound in m/s |

---
# Examples

> ```julia
> rho(f) = 1.20 # air density in kg/m³
> c(f) = 343 # speed of sound in m/s
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Physical quantities can by complex to model sound absorption:
> ```julia
> rho(f) = 1.20 + 0.001*f
> c(f) = 340 + 0.1*sqrt(f)
> add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_volume_material!(s, physical_group=[4,6], density=1.2, speed_of_sound=343)
> ```
"""
function add_volume_material!( 
    simulation::Simulation{Fem, Helmholtz};
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    density::Fdpq,
    speed_of_sound::Fdpq
)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_volume_material!(
                simulation,
                physical_group = pg,
                density = density,
                speed_of_sound = speed_of_sound
            )
        end
        return
    end
    density = _pqv_to_function(density)
    speed_of_sound = _pqv_to_function(speed_of_sound)
    _cpp_add_volume_material!(
        simulation._m._cpp_simulation,
        UInt64(physical_group),
        _cmplx_split_and_store(density)...,
        _cmplx_split_and_store(speed_of_sound)...
    )
end

"""
Assigns a specific acoustic impedance boundary condition to a surface in the mesh. This is used to model absorbers, reflecting surfaces or other boundary treatments.

The specific acoustic impedance is the ratio of complex amplitude of acoustic pressure to complex amplitude of normal particle velocity (in Pa·s/m).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of the boundary surface |
| `specific_acoustic_impedance` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | specific surface acoustic impedance (Pa·s/m) |

---
# Examples

> ```julia
> Z(f) = 1f + 2im
> add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)
> ```

> Multiple physical groups can be assigned at once:
> ```julia
> add_surface_material!(s, physical_group=[3,1], specific_acoustic_impedance=1.0)
> ```
"""
function add_surface_material!( 
    simulation::Simulation{Fem, Helmholtz};
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Fdpq
)
    if physical_group isa AbstractVector
        for pg in physical_group
            add_surface_material!(
                simulation,
                physical_group = pg,
                specific_acoustic_impedance = specific_acoustic_impedance
            )
        end
        return
    end
    specific_acoustic_impedance = _pqv_to_function(specific_acoustic_impedance)
    _cpp_add_surface_material!(
        simulation._m._cpp_simulation,
        UInt64(physical_group),
        _cpp_PhysicalQuantity_impedance,
        _cmplx_split_and_store(specific_acoustic_impedance)...
    )
end

"""
Sources can be applied either at a specific point in space (via `coordinates`) or over a surface region (via `physical_group`). Three excitation types are supported: `volume_velocity`, `particle velocity`, and `pressure`.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `coordinates` | `Vector{Real}`, `Vector{Vector{Real}}` | `[x, y, z]` location of a point source in m |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of a surface or volume region |
| `volume_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Volume velocity in m³/s |
| `particle_velocity` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Normal particle velocity in m/s |
| `pressure` | [frequency-dependent physical quantity](@ref "Frequency-dependent physical quantities") | Acoustic pressure in Pa |

---
# Examples

> Volume velocity source (monopole):
> ```julia
> Q(f) = 10/f # Volume velocity in m³/s as a function of frequency
> add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
> ```
> Suitable for modeling punctual omnidirectional sources.

> Particle velocity source (vibrating surface):
> ```julia
> U(f) = 15/f # Particle velocity in m/s
> add_sound_source!(s, physical_group=2, particle_velocity=U)
> ```
> Prescribes the normal component of particle velocity on all surfaces of a physical group. Useful for modeling vibrating panels or pistons.

> Pressure source:
> ```julia
> P(f) = 2f # Pressure in Pa as a function of frequency
> 
> # At a point in space
> add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
> 
> # On a mesh surface
> add_sound_source!(s, physical_group=3, pressure=P)
> ```
> Prescribes acoustic pressure, either at a point or on a surface.

> Multiple points or physical groups can be assigned at once:
> ```julia
> p1 = [1.0, 3.0, 2.0]
> p2 = [3.0, 1.0, 1.0]
> add_sound_source!(s, coordinates=[p1,p2], volume_velocity=0.01)
> add_sound_source!(s, physical_group=[3,5,9,2], particle_velocity=0.01)
> ```

!!! note
    Each call to `add_sound_source!` should specify either `coordinates` or `physical_group` (not both), and exactly **one** excitation type keyword.
"""
function add_sound_source!( 
    simulation::Simulation{Fem, Helmholtz};
    coordinates::Union{
        AbstractVector{<:Real},
        AbstractVector{<:AbstractVector{<:Real}},
        Nothing
    } = nothing,
    physical_group::Union{
        Integer,
        AbstractVector{<:Integer},
        Nothing
    } = nothing,
    volume_velocity::Union{Fdpq, Nothing} = nothing,
    particle_velocity::Union{Fdpq, Nothing} = nothing,
    pressure::Union{Fdpq, Nothing} = nothing
)
    if isnothing(coordinates) && isnothing(physical_group)
        throw(ArgumentError(
            "`coordinates` and `physical_group` not defined"
        ))
    end
    if !isnothing(coordinates) && !isnothing(physical_group)
        throw(ArgumentError(
            "`coordinates` and `physical_group` defined simultaneously"
        ))
    end
    pq_count = count(
        !isnothing, (volume_velocity, particle_velocity, pressure)
    )
    if pq_count == 0
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` and `pressure` not defined"
        ))
    end
    if pq_count > 1
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` or `pressure` defined " *
            "simultaneously"
        ))
    end
    if coordinates isa AbstractVector{<:Real} && length(coordinates) != 3
        throw(ArgumentError(
            "x,y,z coordinates does not have 3 components"
        ))
    end
    if coordinates isa AbstractVector{<:AbstractVector}
        for c in coordinates
            add_sound_source!(
                simulation,
                coordinates = c,
                physical_group = physical_group,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end
    if physical_group isa AbstractVector
        for pg in physical_group
            add_sound_source!(
                simulation,
                coordinates = coordinates,
                physical_group = pg,
                volume_velocity = volume_velocity,
                particle_velocity = particle_velocity,
                pressure = pressure
            )
        end
        return
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    pqv = Ref{Fdpq}()
    if !isnothing(volume_velocity)
        pq_type = _cpp_PhysicalQuantity_volume_velocity
        pqv[] = volume_velocity
    elseif !isnothing(particle_velocity)
        pq_type = _cpp_PhysicalQuantity_particle_velocity
        pqv[] = particle_velocity
    elseif !isnothing(pressure)
        pq_type = _cpp_PhysicalQuantity_pressure
        pqv[] = pressure
    end

    pqv[] = _pqv_to_function(pqv[])
    source_args =
    if !isnothing(coordinates)
        (_cpp_SourceType_point, Float64.(coordinates))
    elseif !isnothing(physical_group)
        (_cpp_SourceType_surface, UInt64(physical_group))
    end

    _cpp_add_sound_source!(
        simulation._m._cpp_simulation,
        source_args...,
        pq_type,
        _cmplx_split_and_store(pqv[])...
    )
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
    simulation::Simulation{Fem, Helmholtz},
    path_to_hdf5_file::AbstractString
)
    _cpp_set_result_export_path!(
        simulation._m._cpp_simulation,
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
function run!(simulation::Simulation{Fem, Helmholtz})
    _cpp_run!(simulation._m._cpp_simulation)
end

end # module Numav
