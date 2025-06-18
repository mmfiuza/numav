# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

module Numav

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
    
    quantity_value_ref = Ref{Function}()
    safe_cfunction_expr = quote
        quantity_value_cfunc = CxxWrap.@safe_cfunction(quantity_value_ref, Float64, (Float64,))
    end
    function add_source(
        simulation::Simulation{
            Phenomenon.acoustic,
            NumericalMethod.fem,
            Domain.frequency,
            Dimension.d3
        },
        source_type::TypeOfSource.TypeOfSource_type,
        source_coordinates::Vector{Float64},
        physical_quantity::PhysicalQuantity.PhysicalQuantity_type,
        quantity_value::Function
    )
        quantity_value_ref[] = quantity_value
        eval(safe_cfunction_expr)
        _add_source(
            simulation, source_type, source_coordinates, physical_quantity, quantity_value_cfunc)
    end

end