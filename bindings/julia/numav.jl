module Numav
    using CxxWrap
    using Libdl

    function get_lib_path()
        libname = Sys.iswindows() ? "numav_julia" : "libnumav_julia"
        return joinpath(@__DIR__, "..", "..", "build", "lib", libname)
    end

    @wrapmodule(get_lib_path)

    function __init__()
        @initcxx
    end
end