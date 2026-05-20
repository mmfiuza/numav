# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

# export enums
export 
    Phenomenon,
    NumericalMethod,
    Domain,
    Dimension,
    TypeOfSource,
    PhysicalQuantity,
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
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_Phenomenon)
    function __init__()
        @initcxx
    end
end

# wrap the NumericalMethod enum class
module NumericalMethod
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_NumericalMethod)
    function __init__()
        @initcxx
    end
end

# wrap the Domain enum class
module Domain
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_Domain)
    function __init__()
        @initcxx
    end
end

# wrap the Dimension enum class
module Dimension
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_Dimension)
    function __init__()
        @initcxx
    end
end

# wrap the TypeOfSource enum class
module TypeOfSource
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_TypeOfSource)
    function __init__()
        @initcxx
    end
end

# wrap the PhysicalQuantity enum class
module PhysicalQuantity
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_PhysicalQuantity)
    function __init__()
        @initcxx
    end
end

# wrap the ElementOrder enum class
module ElementOrder
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_ElementOrder)
    function __init__()
        @initcxx
    end
end

# wrap the FrequencySamplingDensity enum class
module FrequencySamplingDensity
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl, :define_module_FrequencySamplingDensity)
    function __init__()
        @initcxx
    end
end

# wrap the main Numav module
using CxxWrap
using libnumav_jl_jll
@wrapmodule(() -> libnumav_jl)
function __init__()
    @initcxx
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
        ElementOrder.o1
    };
    physical_group::Int,
    density::Function,
    sound_speed::Function
)
    density_real, density_imag = _cmplx_split_and_store(density)
    sound_speed_real, sound_speed_imag = _cmplx_split_and_store(sound_speed)

    # call the C++ method
    _add_volume_material(
        simulation,
        physical_group,
        density_real,
        density_imag,
        sound_speed_real,
        sound_speed_imag
    )
end

function add_sound_source( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementOrder.o1
    };
    coordinates::Union{Vector, _Empty} = _empty,
    physical_group::Union{Int, _Empty}  = _empty,
    volume_velocity::Union{Function, _Empty} = _empty,
    particle_velocity::Union{Function, _Empty} = _empty,
    pressure::Union{Function, _Empty} = _empty
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

    # Check if velocity or pressure was given
    pq_func = Ref{Function}()
    if volume_velocity != _empty
        pq_type = PhysicalQuantity.volume_velocity
        pq_func[] = volume_velocity
    elseif particle_velocity != _empty
        pq_type = PhysicalQuantity.particle_velocity
        pq_func[] = particle_velocity
    elseif pressure != _empty
        pq_type = PhysicalQuantity.pressure
        pq_func[] = pressure
    end

    pq_func_real, pq_func_imag = _cmplx_split_and_store(pq_func[])

    # call the C++ method
    if coordinates != _empty
        _add_sound_source(
            simulation,
            TypeOfSource.point,
            coordinates,
            pq_type, 
            pq_func_real,
            pq_func_imag
        )
    else
        _add_sound_source(
            simulation,
            TypeOfSource.surface,
            physical_group,
            pq_type,
            pq_func_real,
            pq_func_imag
        )
    end
end

function add_surface_material( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementOrder.o1
    };
    physical_group::Int,
    impedance::Function,
)
    pq_func_real, pq_funq_imag = _cmplx_split_and_store(impedance)

    # call the C++ method
    _add_surface_material(
        simulation,
        physical_group,
        PhysicalQuantity.impedance,
        pq_func_real,
        pq_funq_imag
    )
end

end # module Numav
