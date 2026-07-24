# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

module Numav

    include("cpp.jl") # wrapped C++ part
    include("utils.jl") # various things
    include("options.jl") # define every Option singleton
    include("simulation.jl") # general definition of the Simulation type

    include("mesh.jl")
    include("frequency.jl")
    include("fem-helmholtz.jl")
    
end
