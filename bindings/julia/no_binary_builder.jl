# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

module Numav

    using CxxWrap

    @wrapmodule(() -> libnumav_jl)

    function __init__()
        @initcxx
    end
    
end