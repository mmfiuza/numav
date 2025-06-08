# Copyright (c) 2025 Matheus Machado Fiuza <matheus.fiuza@eac.ufsm.br>

lines = readlines("/root/.julia/dev/libnumav_jl_jll/Project.toml")

for line in lines
    if(startswith(line, "uuid"))
        global uuid = strip(line[5:end], ['\"', ' ', '='])
        break
    end
end

open("/root/.julia/artifacts/Overrides.toml", "w") do file
    write(file, "[" * uuid * "]\n")
    write(file, "libnumav_jl = \"/workspaces/numav/install\"\n")
end
