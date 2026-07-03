# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

# export structs
export Simulation

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

# type aliases
const Pq = Union{Function, Number, AbstractString}
const SimulationFemHelmTet{O} = Simulation{
    NumericalMethod_fem,
    Equation_helmholtz,
    ElementShape_tetrahedron,
    O
}

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

function create_simulation(;
    numerical_method::Symbol,
    equation::Symbol,
    element_shape::Symbol,
    element_order::Symbol
)
    args = ()
    if numerical_method == :fem
        args = (args..., NumericalMethod_fem)
    end
    if equation == :helmholtz 
        args = (args..., Equation_helmholtz)
    end
    if element_shape == :tetrahedron
        args = (args..., ElementShape_tetrahedron)
    end
    if element_order == :linear
        args = (args..., ElementOrder_linear)
    end
    if element_order == :quadratic
        args = (args..., ElementOrder_quadratic)
    end
    return Simulation{args...}()
end

function cubic_range(start::Real, stop::Real, length::Integer)
    @assert length != 0 && length != 1
    x = range(0, 1, length=length)
    return @. cbrt( start^3 + (stop^3 - start^3) * x );
end

function set_frequency!( 
    simulation::SimulationFemHelmTet{O};
    max::Union{Real, Nothing} = nothing,
    vector::Union{AbstractVector{<:Real}, Nothing} = nothing,
    min::Union{Real, Nothing} = nothing,
    length::Union{Integer, Nothing} = nothing,
    sampling_density::Union{Symbol, Nothing} = nothing,
    step::Union{Real, Nothing} = nothing,
) where O
    if isnothing(max) && isnothing(vector)
        throw(ArgumentError("`max` or `vector` must be passed"))
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
                "If `vector` is passed, no other parameter is allowed"
            ))
        end
    end
    if !isnothing(step) && !isnothing(length)
        throw(ArgumentError(
            "`step` and `length` cannot be passed simultaneously"
        ))
    end
    if !isnothing(step) && !isnothing(sampling_density)
        throw(ArgumentError(
            "`step` and `sampling_density` cannot be passed simultaneously"
        ))
    end
    if !isnothing(max)
        min = something(min, 0)
        if !isnothing(step)
            vector = range(min, max, step=step)
        else
            length = something(length, 4096)
            sampling_density = something(sampling_density, :quadratic)
            if sampling_density == :constant
                vector = range(min, max, length=length)
            elseif sampling_density == :quadratic
                vector = cubic_range(min, max, length)
            else
                throw(ArgumentError(
                    "`sampling_density` must be `:constant` or `:quadratic`"
                ))
            end
        end
    end
    _cpp_set_frequency_vector!(simulation, Float64.(Vector(vector)))
end

function load_mesh!(
    simulation::SimulationFemHelmTet{O},
    path_to_mesh::AbstractString
) where O
    _cpp_load_mesh!(simulation, String(path_to_mesh))
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
    throw(ArgumentError("type of `pqv` could not be handled"))
end

function add_volume_material!( 
    simulation::SimulationFemHelmTet{O};
    physical_group::Integer,
    density::Pq,
    speed_of_sound::Pq
) where O
    density = _pqv_to_function(density)
    speed_of_sound = _pqv_to_function(speed_of_sound)
    _cpp_add_volume_material!(
        simulation,
        UInt64(physical_group),
        _cmplx_split_and_store(density)...,
        _cmplx_split_and_store(speed_of_sound)...
    )
end

function add_surface_material!( 
    simulation::SimulationFemHelmTet{O};
    physical_group::Integer,
    specific_acoustic_impedance::Pq
) where O
    specific_acoustic_impedance = _pqv_to_function(specific_acoustic_impedance)
    _cpp_add_surface_material!(
        simulation,
        UInt64(physical_group),
        PhysicalQuantity_impedance,
        _cmplx_split_and_store(specific_acoustic_impedance)...
    )
end

function add_sound_source!( 
    simulation::SimulationFemHelmTet{O};
    coordinates::Union{AbstractVector{<:Real}, Nothing} = nothing,
    physical_group::Union{Integer, Nothing} = nothing,
    volume_velocity::Union{Pq, Nothing} = nothing,
    particle_velocity::Union{Pq, Nothing} = nothing,
    pressure::Union{Pq, Nothing} = nothing
) where O
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
    if !isnothing(coordinates) && length(coordinates) != 3
        throw(ArgumentError(
            "`coordinates` vector does not have 3 elements"
        ))
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    pqv = Ref{Pq}()
    if !isnothing(volume_velocity)
        pq_type = PhysicalQuantity_volume_velocity
        pqv[] = volume_velocity
    elseif !isnothing(particle_velocity)
        pq_type = PhysicalQuantity_particle_velocity
        pqv[] = particle_velocity
    elseif !isnothing(pressure)
        pq_type = PhysicalQuantity_pressure
        pqv[] = pressure
    end

    pqv[] = _pqv_to_function(pqv[])
    source_args =
    if !isnothing(coordinates)
        (SourceType_point, Float64.(coordinates))
    elseif !isnothing(physical_group)
        (SourceType_surface, UInt64(physical_group))
    end

    _cpp_add_sound_source!(
        simulation,
        source_args...,
        pq_type,
        _cmplx_split_and_store(pqv[])...
    )
end

function set_result_export_path!(
    simulation::SimulationFemHelmTet{O},
    path_to_hdf5_file::AbstractString
) where O
    _cpp_set_result_export_path!(simulation, String(path_to_hdf5_file))
end

function run!(simulation::SimulationFemHelmTet{O}) where O
    _cpp_run!(simulation)
end

end # module Numav
