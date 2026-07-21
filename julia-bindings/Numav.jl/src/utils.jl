# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# include dependencies
using Interpolations
using DelimitedFiles

# type aliases
const Fdpq = Union{Function, Number, AbstractString}

function _read_fdpq_table(filename::AbstractString)
    data = readdlm(filename, ',')
    freq = Float64.(data[:, 1])
    fdpqv = Complex.(Float64.(data[:, 2]), Float64.(data[:, 3]))
    return (freq, fdpqv)
end

function _fdpq_to_function(fdpqv::Fdpq)::Function
    if fdpqv isa Number
        return (_ -> fdpqv)
    elseif fdpqv isa AbstractString
        interp = linear_interpolation(_read_fdpq_table(fdpqv)...) 
        return (x -> interp(x))
    elseif fdpqv isa Function
        return fdpqv
    end
    throw(ArgumentError("Invalid type of `pqv`"))
end

function _cubic_range(start::Real, stop::Real, length::Integer)
    @assert length != 0 && length != 1
    x = range(0, 1, length=length)
    return @. cbrt( start^3 + (stop^3 - start^3) * x );
end
