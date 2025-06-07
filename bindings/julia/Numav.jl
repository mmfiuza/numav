module Numav
    using CxxWrap
    using Libdl

    using CxxWrap
    using Numav_jll

    @wrapmodule Numav

    function __init__()
        @initcxx
    end
end