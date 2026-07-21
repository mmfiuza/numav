# Frequency-dependent physical quantities

All material properties and source amplitudes are **frequency-dependent physical quantities**, such as in the functions:
- [`add_volume_material!`](@ref)
- [`add_surface_material!`](@ref)
- [`add_sound_source!`](@ref)

There are three ways to specify them:

## 1. Function

It takes a real frequency value in Hz as its argument and returns some value that can be complex or not:

```julia
# traditional function
function c(f)
    return 343.0 + 0.5im * f^2/1000
end

# function in assignment form
rho(f) = 1.2 + 0.01im*f

# anonymous function
f -> 343.0 + 0.05 * sqrt(f)

add_volume_material!(s, physical_group=1, density=c, speed_of_sound=rho)
```

## 2. Constant

A plain scalar value (complex or not). Numav will treat it as uniform across all frequencies:

```julia
add_sound_source!(s, coordinates=[0.0, 0.0, 0.0], pressure=0.1)
```

## 3. Tabulated file (path string)

A `String` containing the path to a plain-text CSV-like file with three columns: `frequency`, `real`, and `imag`. Numav reads the table and interpolates between the provided data points at each required frequency step. There is no header row.

```
  0, 1, 2
 25, 1, 2
 50, 5, 1
 75, 3, 2
100, 1, 0
```

```julia
add_sound_source!(s, physical_group=1, pressure="data/pressure.txt")
```

The file format rules are:
- Three columns separated by commas: `frequency` (Hz), `real` part, `imaginary` part.
- No header row.
- For purely real quantities (e.g. density), set the `imag` column to `0`.

!!! warning
    Any table must contain data covering all the frequency range that you plan to simulate.
