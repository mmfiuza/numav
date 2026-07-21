# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# include dependencies
using numav_julia_jll
using CxxWrap
using Base.Threads

# wrap the C++ part
@wrapmodule(() -> libnumav_julia)
function __init__()
    @initcxx
    # Set env variables to enable MKL in SDL mode to configure at runtime
    ENV["MKL_INTERFACE_LAYER"] = "ILP64"
    ENV["MKL_THREADING_LAYER"] = "INTEL"
end

# global storage for user-defined functions
const _user_functions::Vector{Ref{Function}} = [ ]
global _next_index::UInt = 1
_lock = ReentrantLock()

function _cmplx_split_and_store(f::Function)
    idx = lock(_lock) do
        global _next_index
        idx = _next_index
        _next_index += 1
        push!(_user_functions, f)
        return idx
    end
    name_real = Symbol("_f", idx, "_r")
    name_imag = Symbol("_f", idx, "_i")
    @eval begin
        ($name_real)(x::Float64) = Float64(real(_user_functions[$idx][](x)))
        ($name_imag)(x::Float64) = Float64(imag(_user_functions[$idx][](x)))
    end
    f_real = @eval CxxWrap.@safe_cfunction($(name_real), Float64, (Float64,))
    f_imag = @eval CxxWrap.@safe_cfunction($(name_imag), Float64, (Float64,))
    return (f_real, f_imag)
end
