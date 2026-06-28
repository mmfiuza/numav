# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

# export enums
export 
    Phenomenon,
    NumericalMethod,
    Domain,
    Dimension,
    SourceType,
    PhysicalQuantity,
    ElementShape,
    ElementOrder,
    FrequencySamplingDensity

# export structs
export Simulation

# export functions
export
    set_maximum_frequency,
    set_frequency_range,
    set_frequency_steps_count,
    set_frequency_sampling_density,
    set_frequency_steps,
    load_mesh,
    add_volume_material,
    add_sound_source,
    add_surface_material,
    set_result_export_path,
    simulate

# wrap the Phenomenon enum class
module Phenomenon
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_Phenomenon)
    function __init__()
        @initcxx
    end
end

# wrap the NumericalMethod enum class
module NumericalMethod
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_NumericalMethod)
    function __init__()
        @initcxx
    end
end

# wrap the Domain enum class
module Domain
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_Domain)
    function __init__()
        @initcxx
    end
end

# wrap the Dimension enum class
module Dimension
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_Dimension)
    function __init__()
        @initcxx
    end
end

# wrap the SourceType enum class
module SourceType
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_SourceType)
    function __init__()
        @initcxx
    end
end

# wrap the PhysicalQuantity enum class
module PhysicalQuantity
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_PhysicalQuantity)
    function __init__()
        @initcxx
    end
end

# wrap the ElementShape enum class
module ElementShape
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_ElementShape)
    function __init__()
        @initcxx
    end
end

# wrap the ElementOrder enum class
module ElementOrder
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_ElementOrder)
    function __init__()
        @initcxx
    end
end

# wrap the FrequencySamplingDensity enum class
module FrequencySamplingDensity
    using CxxWrap
    using numav_julia_jll
    @wrapmodule(() -> libnumav_julia, :define_module_FrequencySamplingDensity)
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
global _next_index::Int = 1
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

function add_volume_material( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementShape.tetrahedron,
        ElementOrder.o1
    };
    physical_group::Int,
    density::Union{Function, Number, String},
    sound_speed::Union{Function, Number, String}
)
    density_args =
    if density isa Function
        _cmplx_split_and_store(density)
    elseif density isa Number
        (ComplexF64(density),)
    else
        (density,)
    end

    sound_speed_args =
    if sound_speed isa Function
        _cmplx_split_and_store(sound_speed)
    elseif sound_speed isa Number
        (ComplexF64(sound_speed),)
    else
        (sound_speed,)
    end

    _add_volume_material(
        simulation, physical_group, density_args..., sound_speed_args...
    )
end

function add_surface_material( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementShape.tetrahedron,
        ElementOrder.o1
    };
    physical_group::Int,
    impedance::Union{Function, Number, String, _Empty},
)
    impedance_args =
    if impedance isa Function
        _cmplx_split_and_store(impedance)
    elseif impedance isa Number
        (ComplexF64(impedance),)
    elseif impedance isa String
        (impedance,)
    end

    _add_surface_material(
        simulation,
        physical_group,
        PhysicalQuantity.impedance,
        impedance_args...
    )
end

function add_sound_source( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementShape.tetrahedron,
        ElementOrder.o1
    };
    coordinates::Union{Vector, _Empty} = _empty,
    physical_group::Union{Int, _Empty}  = _empty,
    volume_velocity::Union{Function, Number, String, _Empty} = _empty,
    particle_velocity::Union{Function, Number, String, _Empty} = _empty,
    pressure::Union{Function, Number, String, _Empty} = _empty
)
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
        (SourceType.point, coordinates)
    elseif physical_group != _empty
        (SourceType.surface, physical_group)
    end

    _add_sound_source(simulation, source_args..., pq_type, pq_args...)
end

end # module Numav
