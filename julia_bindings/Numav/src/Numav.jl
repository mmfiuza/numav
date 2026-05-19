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
    add_specific_surface_acoustic_impedance,
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

struct _Empty end
const _empty = _Empty()

function _split_cmplx_func(func::Function)
    real_cfunc = CxxWrap.@safe_cfunction(
        x -> Float64(real(func(x))),
        Float64,
        (Float64,)
    )
    imag_cfunc = CxxWrap.@safe_cfunction(
        x -> Float64(imag(func(x))),
        Float64,
        (Float64,)
    )
    return (real_cfunc, imag_cfunc)
end

function add_volume_material( 
    simulation::Simulation{
        Phenomenon.acoustic,
        NumericalMethod.fem,
        Domain.frequency,
        Dimension.d3,
        ElementOrder.o1
    },
    physical_group::Int,
    rho::Function,
    c::Function
)
    rho_real, rho_imag = _split_cmplx_func(rho)
    c_real, c_imag = _split_cmplx_func(c)

    # call the C++ function
    _add_volume_material(
        simulation, physical_group, rho_real, rho_imag, c_real, c_imag
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
    physical_quantity_function = Ref{Function}()
    if volume_velocity != _empty
        physical_quantity_type = PhysicalQuantity.volume_velocity
        physical_quantity_function[] = volume_velocity
    elseif particle_velocity != _empty
        physical_quantity_type = PhysicalQuantity.particle_velocity
        physical_quantity_function[] = particle_velocity
    elseif pressure != _empty
        physical_quantity_type = PhysicalQuantity.pressure
        physical_quantity_function[] = pressure
    end

    physical_quantity_real, physical_quantity_imag = 
        _split_cmplx_func(physical_quantity_function[])

    # call the C++ function
    if coordinates != _empty
        _add_sound_source(
            simulation,
            TypeOfSource.point,
            coordinates,
            physical_quantity_type, 
            physical_quantity_real,
            physical_quantity_imag
        )
    else
        _add_sound_source(
            simulation,
            TypeOfSource.surface,
            physical_group,
            physical_quantity_type,
            physical_quantity_real,
            physical_quantity_imag
        )
    end
end

function add_specific_surface_acoustic_impedance( 
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
    impedance_real, impedance_imag = _split_cmplx_func(impedance)
    # call the C++ function
    _add_surface_specific_acoustic_impedance(
        simulation,
        physical_group,
        impedance_real,
        impedance_imag
    )
end

end # module Numav
