module Numav
    using CxxWrap
    using Libdl
    using libnumav_jll

    @wrapmodule libnumav

    function __init__()
        @initcxx
    end
end