# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    add_volume_material!,
    add_surface_material!,
    add_sound_source!

@kwdef mutable struct SimulationFemHelmholtz{
    S<:ElementShape,
    O<:ElementOrder,
} <: Simulation

    numerical_method::FemNumericalMethod = FemNumericalMethod()
    equation::HelmholtzEquation = HelmholtzEquation()
    element_shape::S = S()
    element_order::O = O()

    _cpp_simulation::_cpp_Simulation
    _fi_to_freq::Vector{Float64} = []
end

function _check_if_simulation_has_fem_helmholtz(s::Simulation)
    if (
        !hasproperty(s, :numerical_method) ||
        !hasproperty(s, :equation) ||
        !(
            (s.numerical_method isa FemNumericalMethod) &&
            (s.equation isa HelmholtzEquation)
        )
    )
        _throw_simulation_not_applicable()
    end
    return
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
    simulation::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    density::Fdpq,
    speed_of_sound::Fdpq
)
    _check_if_simulation_has_fem_helmholtz(simulation)
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
    density = _fdpq_to_function(density)
    speed_of_sound = _fdpq_to_function(speed_of_sound)
    _cpp_add_volume_material!(
        simulation._cpp_simulation,
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
    simulation::SimulationFemHelmholtz;
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Fdpq
)
    _check_if_simulation_has_fem_helmholtz(simulation)
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
    specific_acoustic_impedance = _fdpq_to_function(specific_acoustic_impedance)
    _cpp_add_surface_material!(
        simulation._cpp_simulation,
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
    simulation::SimulationFemHelmholtz;
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
    _check_if_simulation_has_fem_helmholtz(simulation)
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
    fdpq_count = count(
        !isnothing, (volume_velocity, particle_velocity, pressure)
    )
    if fdpq_count == 0
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` and `pressure` not defined"
        ))
    end
    if fdpq_count > 1
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
    fdpqv = Ref{Fdpq}()
    if !isnothing(volume_velocity)
        pq_type = _cpp_PhysicalQuantity_volume_velocity
        fdpqv[] = volume_velocity
    elseif !isnothing(particle_velocity)
        pq_type = _cpp_PhysicalQuantity_particle_velocity
        fdpqv[] = particle_velocity
    elseif !isnothing(pressure)
        pq_type = _cpp_PhysicalQuantity_pressure
        fdpqv[] = pressure
    end

    fdpqv[] = _fdpq_to_function(fdpqv[])
    source_args =
    if !isnothing(coordinates)
        (_cpp_SourceType_point, Float64.(coordinates))
    elseif !isnothing(physical_group)
        (_cpp_SourceType_surface, UInt64(physical_group))
    end

    _cpp_add_sound_source!(
        simulation._cpp_simulation,
        source_args...,
        pq_type,
        _cmplx_split_and_store(fdpqv[])...
    )
end
