# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

# export enums
export 
    NumericalMethod,
    Equation,
    ElementShape,
    ElementOrder,
    SourceType,
    PhysicalQuantity,
    FrequencySamplingDensity

# export structs
export Simulation

# export functions
export
    create_simulation,
    set_maximum_frequency!,
    set_frequency_range!,
    set_frequency_steps_count!,
    set_frequency_sampling_density!,
    set_frequency_steps!,
    load_mesh!,
    add_volume_material!,
    add_sound_source!,
    add_surface_material!,
    set_result_export_path!,
    run!

# wrap the NumericalMethod enum class
module NumericalMethod
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :NumericalMethod_module)
    function __init__()
        @initcxx
    end
end

# wrap the Equation enum class
module Equation
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :Equation_module)
    function __init__()
        @initcxx
    end
end

# wrap the ElementShape enum class
module ElementShape
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :ElementShape_module)
    function __init__()
        @initcxx
    end
end

# wrap the ElementOrder enum class
module ElementOrder
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :ElementOrder_module)
    function __init__()
        @initcxx
    end
end

# wrap the SourceType enum class
module SourceType
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :SourceType_module)
    function __init__()
        @initcxx
    end
end

# wrap the PhysicalQuantity enum class
module PhysicalQuantity
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :PhysicalQuantity_module)
    function __init__()
        @initcxx
    end
end

# wrap the FrequencySamplingDensity enum class
module FrequencySamplingDensity
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :FrequencySamplingDensity_module)
    function __init__()
        @initcxx
    end
end

# wrap the main Numav module
using CxxWrap
using numav_julia_jll
@wrapmodule(() -> libnumav_julia)
function __init__()
    @initcxx
    # Set env variables to enable MKL in SDL mode to configure at runtime
    ENV["MKL_INTERFACE_LAYER"] = "ILP64"
    ENV["MKL_THREADING_LAYER"] = "INTEL"
end

using Base.Threads

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

struct _Empty end
const _empty = _Empty()

function create_simulation(;
    numerical_method::Symbol,
    equation::Symbol,
    element_shape::Symbol,
    element_order::Symbol
)
    args = ()
    if numerical_method == :fem
        args = (args..., NumericalMethod.fem)
    end
    if equation == :helmholtz 
        args = (args..., Equation.helmholtz)
    end
    if element_shape == :tetrahedron
        args = (args..., ElementShape.tetrahedron)
    end
    if element_order == :linear
        args = (args..., ElementOrder.linear)
    end
    if element_order == :quadratic
        args = (args..., ElementOrder.quadratic)
    end
    return Simulation{args...}()
end

function add_volume_material!( 
    simulation::Simulation{
        NumericalMethod.fem,
        Equation.helmholtz,
        ElementShape.tetrahedron,
        O
    };
    physical_group::Integer,
    density::Union{Function, Number, String},
    speed_of_sound::Union{Function, Number, String}
) where O
    density_args =
    if density isa Function
        _cmplx_split_and_store(density)
    elseif density isa Number
        (ComplexF64(density),)
    else
        (density,)
    end

    speed_of_sound_args =
    if speed_of_sound isa Function
        _cmplx_split_and_store(speed_of_sound)
    elseif speed_of_sound isa Number
        (ComplexF64(speed_of_sound),)
    else
        (speed_of_sound,)
    end

    _add_volume_material!(
        simulation,
        UInt64(physical_group),
        density_args...,
        speed_of_sound_args...
    )
end

function add_surface_material!( 
    simulation::Simulation{
        NumericalMethod.fem,
        Equation.helmholtz,
        ElementShape.tetrahedron,
        O
    };
    physical_group::Integer,
    impedance::Union{Function, Number, String, _Empty},
) where O
    impedance_args =
    if impedance isa Function
        _cmplx_split_and_store(impedance)
    elseif impedance isa Number
        (ComplexF64(impedance),)
    elseif impedance isa String
        (impedance,)
    end

    _add_surface_material!(
        simulation,
        UInt64(physical_group),
        PhysicalQuantity.impedance,
        impedance_args...
    )
end

function add_sound_source!( 
    simulation::Simulation{
        NumericalMethod.fem,
        Equation.helmholtz,
        ElementShape.tetrahedron,
        O
    };
    coordinates::Union{Vector, _Empty} = _empty,
    physical_group::Union{Integer, _Empty} = _empty,
    volume_velocity::Union{Function, Number, String, _Empty} = _empty,
    particle_velocity::Union{Function, Number, String, _Empty} = _empty,
    pressure::Union{Function, Number, String, _Empty} = _empty
) where O
    if coordinates == _empty && physical_group == _empty
        throw(ArgumentError(
            "`coordinates` and `physical_group` not defined"
        ))
    end
    if coordinates != _empty && physical_group != _empty
        throw(ArgumentError(
            "`coordinates` and `physical_group` defined simultaneously"
        ))
    end
    non_empty_count = 0
    if volume_velocity != _empty
        non_empty_count += 1
    end
    if particle_velocity != _empty
        non_empty_count += 1
    end
    if pressure != _empty
        non_empty_count += 1
    end
    if non_empty_count == 0
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` and `pressure` not defined"
        ))
    end
    if non_empty_count > 1
        throw(ArgumentError(
            "`volume_velocity`, `particle_velocity` or `pressure` defined " *
            "simultaneously"
        ))
    end
    if coordinates != _empty && length(coordinates) != 3
        throw(ArgumentError(
            "`coordinates` vector does not have 3 elements"
        ))
    end

    # Check if volume_velocity, particle_velocity or pressure was given
    pqv = Ref{Any}()
    if volume_velocity != _empty
        pq_type = PhysicalQuantity.volume_velocity
        pqv[] = volume_velocity
    elseif particle_velocity != _empty
        pq_type = PhysicalQuantity.particle_velocity
        pqv[] = particle_velocity
    elseif pressure != _empty
        pq_type = PhysicalQuantity.pressure
        pqv[] = pressure
    end

    pq_args =
    if pqv[] isa Function
        _cmplx_split_and_store(pqv[])
    elseif pqv[] isa Number
        (ComplexF64(pqv[]),)
    elseif pqv[] isa String
        (pqv[],)
    end

    source_args =
    if coordinates != _empty
        (SourceType.point, Float64.(coordinates))
    elseif physical_group != _empty
        (SourceType.surface, UInt64(physical_group))
    end

    _add_sound_source!(simulation, source_args..., pq_type, pq_args...)
end

end # module Numav
