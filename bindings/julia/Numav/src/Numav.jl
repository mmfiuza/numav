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
    
    # declare the expression to crate safa_cfuntions - type: Float64(Float64)
    function_ref_Float64_Float64 = Ref{Function}()
    safe_cfunction_expr_Float64_Float64 = quote
        CxxWrap.@safe_cfunction(function_ref_Float64_Float64[], Float64, (Float64,))
    end

    function add_source(
        simulation::Simulation{
            Phenomenon.acoustic, NumericalMethod.fem,
            Domain.frequency, Dimension.d3},
        source_type::TypeOfSource.type,
        source_coordinates::Vector{Float64},
        physical_quantity::PhysicalQuantity.type,
        quantity_value::Function
    )
        # real part function split
        function_ref_Float64_Float64[] = 
            function real_part(n::Float64) return Float64(real(quantity_value(n))) end
        real_safe_cfunc = eval(safe_cfunction_expr_Float64_Float64)

        # imaginary part function split
        function_ref_Float64_Float64[] =
            function imag_part(n::Float64) return Float64(imag(quantity_value(n))) end
        imag_safe_cfunc = eval(safe_cfunction_expr_Float64_Float64)

        # call C++ function
        _add_source(
            simulation, source_type, source_coordinates, physical_quantity,
            real_safe_cfunc, imag_safe_cfunc
        )
    end

end