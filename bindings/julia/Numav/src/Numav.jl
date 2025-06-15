# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

module Numav

    # wrap the main Numav module
    using CxxWrap
    using libnumav_jl_jll
    @wrapmodule(() -> libnumav_jl)
    function __init__()
        @initcxx
    end

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

end