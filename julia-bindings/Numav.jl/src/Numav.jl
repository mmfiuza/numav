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
    add_sound_source!,
    add_surface_material!,
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
const Pq = Union{Function, Number, AbstractString}

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

function _pqv_to_function(pqv::Pq)::Function
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

function add_volume_material!( 
    simulation::Simulation{Fem, Helmholtz};
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    density::Pq,
    speed_of_sound::Pq
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

function add_surface_material!( 
    simulation::Simulation{Fem, Helmholtz};
    physical_group::Union{Integer, AbstractVector{<:Integer}},
    specific_acoustic_impedance::Pq
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
    volume_velocity::Union{Pq, Nothing} = nothing,
    particle_velocity::Union{Pq, Nothing} = nothing,
    pressure::Union{Pq, Nothing} = nothing
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
    pqv = Ref{Pq}()
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

function set_result_export_path!(
    simulation::Simulation{Fem, Helmholtz},
    path_to_hdf5_file::AbstractString
)
    _cpp_set_result_export_path!(
        simulation._m._cpp_simulation,
        String(path_to_hdf5_file)
    )
end

function run!(simulation::Simulation{Fem, Helmholtz})
    _cpp_run!(simulation._m._cpp_simulation)
end

end # module Numav
