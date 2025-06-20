# Copyright (c) 2025 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

    # export structs
    export Simulation, Result
    
    # export enums
    export Phenomenon, NumericalMethod, Domain, Dimension, TypeOfSource, PhysicalQuantity
    
    # export functions
    export set_element_order, set_freq_limits, load_mesh, add_volume_material, add_source,
        add_specific_surface_acoustic_impedance, run

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

    # wrap the main Numav module
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl)
    function __init__()
        @initcxx
    end

    struct _Empty end
    const _empty = _Empty()
    
    # declare the expression to crate safa_cfuntions - type: Float64(Float64)
    function_ref_Float64_Float64 = Ref{Function}()
    safe_cfunction_expr_Float64_Float64 = quote
        CxxWrap.@safe_cfunction(function_ref_Float64_Float64[], Float64, (Float64,))
    end

    function add_source( 
        simulation::Simulation{
            Phenomenon.acoustic,
            NumericalMethod.fem,
            Domain.frequency,
            Dimension.d3
        };
        coordinates::Union{Vector,_Empty} = _empty,
        surface_id::Union{Int,_Empty}  = _empty,
        volume_velocity::Union{Function,_Empty} = _empty,
        pressure::Union{Function,_Empty} = _empty
    )
        if coordinates==_empty && surface_id==_empty
            throw(ArgumentError("`coordinates` and `surface_id` not defined"))
        elseif coordinates!=_empty && surface_id!=_empty
            throw(ArgumentError("`coordinates` and `surface_id` defined simultaneously"))
        elseif volume_velocity==_empty && pressure==_empty
            throw(ArgumentError("`volume_velocity` and `pressure` not defined"))
        elseif volume_velocity!=_empty && pressure!=_empty
            throw(ArgumentError("`volume_velocity` and `pressure` defined simultaneously"))
        end

        # Check if velocity or pressure was given
        physical_quantity_function = Ref{Function}()
        if volume_velocity!=_empty
            physical_quantity = PhysicalQuantity.volume_velocity
            physical_quantity_function[] = volume_velocity
        else
            physical_quantity = PhysicalQuantity.pressure
            physical_quantity_function[] = pressure
        end

        # real part function split 
        function_ref_Float64_Float64[] = 
            function real_part(n::Float64)
                return real(physical_quantity_function[](n)) |> Float64
            end
        real_physical_quantity_cfunction = eval(safe_cfunction_expr_Float64_Float64)

        # imaginary part function split
        function_ref_Float64_Float64[] =
            function imag_part(n::Float64)
                return imag(physical_quantity_function[](n)) |> Float64
            end
        imag_physical_quantity_cfunction = eval(safe_cfunction_expr_Float64_Float64)

        # call the C++ function
        if coordinates!=_empty
            _add_source(
                simulation, TypeOfSource.point, coordinates, physical_quantity, 
                real_physical_quantity_cfunction, imag_physical_quantity_cfunction
            )
        else
            _add_source(
                simulation, TypeOfSource.surface, surface_id, physical_quantity,
                real_physical_quantity_cfunction, imag_physical_quantity_cfunction
            )
        end
    end

    function add_specific_surface_acoustic_impedance( 
        simulation::Simulation{
            Phenomenon.acoustic,
            NumericalMethod.fem,
            Domain.frequency,
            Dimension.d3
        };
        surface_id::Int,
        impedance::Function,
    )
        # real part function split 
        function_ref_Float64_Float64[] = 
            function real_part(n::Float64)
                return real(impedance[](n)) |> Float64
            end
        real_impedance_cfunction = eval(safe_cfunction_expr_Float64_Float64)

        # imaginary part function split
        function_ref_Float64_Float64[] =
            function imag_part(n::Float64)
                return imag(impedance[](n)) |> Float64
            end
        imag_impedance_cfunction = eval(safe_cfunction_expr_Float64_Float64)

        # call the C++ function
        _add_surface_specific_acoustic_impedance(
            simulation, surface_id, real_impedance_cfunction, imag_impedance_cfunction
        )
    end

end