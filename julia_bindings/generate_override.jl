# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# get libnumav_jl_jll uuid
lines = readlines("/usr/local/share/julia/dev/libnumav_jl_jll/Project.toml")
for line in lines
    if(startswith(line, "uuid"))
        global uuid = strip(line[5:end], ['\"', ' ', '='])
        break
    end
end

# generate Overrides.toml
open("/usr/local/share/julia/artifacts/Overrides.toml", "w") do file
    write(file, "[" * uuid * "]\n")
    write(file, "libnumav_jl = \"/workspace/julia_bindings/override\"\n")
end
