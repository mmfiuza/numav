# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    set_frequency!,
    get_frequency_vector

function _check_if_simulation_has_frequency(s::Simulation)
    if (
        !hasproperty(s, :equation) ||
        !(s.equation isa HelmholtzEquation)
    )
        _throw_simulation_not_applicable()
    end
    return
end

"""
Sets the frequecies to perform the simulations with two mutually exclusive ways to use it:
- **automatic mode**: pass `max` (and optionally `min`, `length`, `sampling_density` and `step`) and Numav will automatically determine the frequency vector;
- **manual mode**: pass `vector` directly, giving full control over which frequencies are evaluated.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |
| **Keyword arguments** | | |
| `max` | `Real` | Upper frequency limit (Hz). Required unless `steps` is provided. Mutually exclusive with `steps`. |
| `min` | `Real` | Lower bound of the frequency range (Hz). Only usable together with `max`. Defaults to `0` if omitted. |
| `length` | `Integer` | Number of frequency steps to compute within the `min`/`max` range. Only usable together with `max`. Defaults to `4096` if omitted. |
| `sampling_density` | `Numav.Option` | Sampling strategy to use within the range. Only usable together with `max`. Defaults to `Quadratic` if omitted. |
| `step` | `Real` | Difference in Hertz between each frequency step, supposing equally spaced steps. Only usable together with `max`. Cannot be used together with `length` and `sampling_density` |
| `vector` | `Vector{Real}` | List of frequencies in Hertz to solve at. Cannot be used together with `max`, `min`, `length`, `sampling_density` and `step`. |

---
# `sampling_density` options

| Mode | Description |
|:--|:--|
| `Quadratic` (_default_) | Frequency discretization density grows quadratically with frequency. It tends to follow the modal density of rooms, being the recommended option for good accuracy and computation time balance. |
| `Constant` | Frequency steps are evenly spaced (uniform spacing). Suitable for broadband analyses where equal resolution at all frequencies is desired. |
 

---
# Examples

> Solve up to a maximum frequency, letting Numav automatically pick the steps from 0 Hz:
> ```julia
> set_frequency!(s, max=200) # solve from 0 Hz to 200 Hz
> ```
 
> Solve within an explicit range, by adding `min`:
> ```julia
> set_frequency!(s, min=20, max=200) # solve from 20 Hz to 200 Hz
> ```
 
> Control the number of frequency steps with `length` (defaults to `4096` when not specified):
> ```julia
> set_frequency!(s, max=200, length=500) # evaluate at 500 frequency points
> ```
> A higher step count gives finer frequency resolution at the cost of longer computation time.

> Control how steps are distributed with `sampling_density` (defaults to `Quadratic` when not specified):
> ```julia
> set_frequency!(s, max=200, sampling_density=Constant)
> ```

> Manually specify the exact frequencies with `vector`, for example to match measurement points or concentrate resolution around a resonance:
> ```julia
> set_frequency!(s, vector=[10, 40, 60, 100]) # Solve only at these frequencies
> ```

!!! tip
    Use [`get_frequency_vector`](@ref) to visually check the defined frequency vector.
"""
function set_frequency!( 
    simulation::Simulation;
    max::Union{Real, Nothing} = nothing,
    min::Union{Real, Nothing} = nothing,
    length::Union{Integer, Nothing} = nothing,
    sampling_density::Union{Option, Nothing} = nothing,
    step::Union{Real, Nothing} = nothing,
    vector::Union{AbstractVector{<:Real}, Nothing} = nothing
)
    _check_if_simulation_has_frequency(simulation)
    if isnothing(max) && isnothing(vector)
        throw(ArgumentError("Neither `max` nor `vector` passed"))
    end
    if !isnothing(vector)
        if (
            !isnothing(max) ||
            !isnothing(min) ||
            !isnothing(length) ||
            !isnothing(sampling_density) ||
            !isnothing(step)
        )
            throw(ArgumentError(
                "`vector` and other parameter(s) passed simultaneously"
            ))
        end
    end
    if !isnothing(step) && !isnothing(length)
        throw(ArgumentError(
            "`step` and `length` passed simultaneously"
        ))
    end
    if !isnothing(step) && !isnothing(sampling_density)
        throw(ArgumentError(
            "`step` and `sampling_density` passed simultaneously"
        ))
    end
    if !isnothing(max)
        min = something(min, 0)
        if !isnothing(step)
            vector = range(min, max, step=step)
        else
            length = something(length, 4096)
            sampling_density = something(sampling_density, Quadratic)
            if sampling_density === Constant
                vector = range(min, max, length=length)
            elseif sampling_density === Quadratic
                vector = _cubic_range(min, max, length)
            else
                throw(ArgumentError("Invalid `sampling_density` option"))
            end
        end
    end
    simulation._fi_to_freq = Float64.(Vector(vector))
    _cpp_set_frequency_vector!(
        simulation._cpp_simulation,
        simulation._fi_to_freq
    )
end

"""
Returns a copy of the frequency vector set by [`set_frequency!`](@ref).

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | the simulation instance |

---
# Examples

> ```julia
> get_frequency_vector(s)
> ```
"""
function get_frequency_vector(
    simulation::Simulation
)
    _check_if_simulation_has_frequency(simulation)
    return copy(simulation._fi_to_freq)
end
