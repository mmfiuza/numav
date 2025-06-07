module Numav

    using CxxWrap
    using libnumav_jll

    function _get_lib_path()
        return joinpath(libnumav_jll.find_artifact_dir(), "lib", "libnumav_jl")
    end

    @wrapmodule libnumav_jll.get_libnumav_jl_path()

    function __init__()
        @initcxx
    end

end