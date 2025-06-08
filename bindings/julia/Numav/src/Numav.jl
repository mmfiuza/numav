module Numav

    using CxxWrap
    using libnumav_jl_jll

    @wrapmodule(() -> libnumav_jl)

    function __init__()
        @initcxx
    end

end