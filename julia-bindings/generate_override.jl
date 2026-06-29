# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

lib_name = "numav_julia_jll"
override_path = "/workspace/julia-bindings/override"
depot_path = Base.DEPOT_PATH[1]
jll_path = "$(depot_path)/dev/$(lib_name)"


# get uuid
lines = readlines("$(jll_path)/Project.toml")
for line in lines
    if(startswith(line, "uuid"))
        global uuid = strip(line[5:end], ['\"', ' ', '='])
        break
    end
end

# generate Overrides.toml
override_toml_file_path = "$(depot_path)/artifacts/Overrides.toml"
lib_name_without_jll = replace(lib_name, "_jll" => "")
open(override_toml_file_path, "w") do file
    write(file, "[" * uuid * "]\n")
    write(file, "$(lib_name_without_jll) = \"$(override_path)\"\n")
end

println("")
println("Overridden the package in: ", jll_path)
println("")
println("The overrides file is now: ", override_toml_file_path)
println("")
println("The text in this overrides file is:")
println(read(override_toml_file_path, String))
println("")
