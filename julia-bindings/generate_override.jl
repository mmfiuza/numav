# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# get numav_julia_jll uuid
lines = readlines("/usr/local/share/julia/dev/numav_julia_jll/Project.toml")
for line in lines
    if(startswith(line, "uuid"))
        global uuid = strip(line[5:end], ['\"', ' ', '='])
        break
    end
end

# generate Overrides.toml
open("/usr/local/share/julia/artifacts/Overrides.toml", "w") do file
    write(file, "[" * uuid * "]\n")
    write(file, "numav_julia = \"/workspace/julia-bindings/override\"\n")
end
