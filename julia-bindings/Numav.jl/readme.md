# Numav.jl

Numav.jl is a Julia wrapper for the C++ library Numav, to perform acoustics and vibrations. Currently, it only supports frequency domain FEM simulations for acoustics. However, it is a plan to cover more numerical methods soon.

---

## Installation
Just type into the Julia REPL:
```julia
] add Numav
```

---

## Quick start

Here is a complete example on how to run your first simulation in Numav:

```julia
using Numav

# Create the simulation object
s = Simulation{
    Phenomenon.acoustic,
    NumericalMethod.fem,
    Domain.frequency,
    Dimension.d3,
    ElementOrder.o1 # can be o2
}()

# Set simulated frequencies (Hz)
set_maximum_frequency(s, 200)

# Optional: set frequency range
# set_frequency_range(s, 20, 200)

# Optional: manually specify frequencies
# set_frequency_steps(s, [20, 50, 100, 200])

# Load mesh
load_mesh(s, "path/to/mesh_file.bdf")

# Define material (frequency-dependent)
rho(f) = 1.20
c(f) = 343
add_volume_material(s, physical_group=1, density=rho, sound_speed=c)

# Add a monopole volume velocity source at a point
Q(f) = 10/f
add_sound_source(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)

# Add a particle velocity source on a surface
U(f) = 15/f
add_sound_source(s, physical_group=2, particle_velocity=U)

# Add pressure sources (to point or surface)
P(f) = 2f
add_sound_source(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_sound_source(s, physical_group=3, pressure=P)

# Add a surface impedance boundary condition
Z(f) = 1f + 2im
add_surface_material(s, physical_group=4, impedance=Z)

# Set output file
set_result_export_path(s, "result.h5")

# Run the simulation
simulate(s)
```

---

## Frequency-dependent quantities

All material properties and source amplitudes are **frequency-dependent quantities**. Numav evaluates them at each frequency step automatically. There are three ways to specify them:

### 1. Function

It takes a real frequency value in Hz as its argument and returns some value that can be complex or not:

```julia
# Constant function
rho(f) = 1.20

# Analytically frequency-dependent
Q(f) = 10 / f

# Complex-valued
Z(f) = 1.0f + 2.0im

# Anonymous function
c = f -> 343.0 + 0.05 * sqrt(f)

# Regular function
function U(f)
    return 343.0 + 0.5im * f^2/1000
end
```

### 2. Constant

A plain scalar value (complex or not). Numav will treat it as uniform across all frequencies:

```julia
add_volume_material(s, physical_group=1, density=1.20, sound_speed=343)
```

### 3. Tabulated file (path string)

A `String` containing the path to a plain-text CSV-like file with three columns: `frequency`, `real`, and `imag`. Numav reads the table and interpolates between the provided data points at each required frequency step. There is no header row.

```
0,   1, 2
25,  1, 2
50,  5, 1
75,  3, 2
100, 1, 0
```

```julia
add_surface_material(s, physical_group=4, impedance="data/wall_impedance.txt")
```

The file format rules are:
- Three columns separated by commas: `frequency` (Hz), `real` part, `imaginary` part.
- No header row.
- For purely real quantities (e.g. density), set the `imag` column to `0`.

> **Note:** Any table should contain data covering all the frequency range that you plan to simulate.
---

## API reference

### `Simulation{...}()`

Creates the top-level simulation object. All simulation parameters are encoded as type parameters.

```julia
s = Simulation{
    Phenomenon.acoustic, # Physics type    
    NumericalMethod.fem, # Numerical method
    Domain.frequency,    # Analysis domain
    Dimension.d3,        # Spatial dimension
    ElementOrder.o1      # Finite element polynomial order
}()
```

**Type parameters:**

| Parameter | Supported options | Description |
|---|---|---|
| `Phenomenon` | `acoustic` | The physical phenomenon being simulated |
| `NumericalMethod` | `fem` | Numerical method used |
| `Domain` | `frequency` | Whether to solve in the frequency or time domain |
| `Dimension` | `d3` | Spatial dimensionality of the problem |
| `ElementOrder` | `o1`, `o2` | Polynomial order of the finite elements (1st or 2nd order) |

---

### Frequency specification

There are three functions to specify which frequencies the solver should evaluate:
- `set_maximum_frequency` (recommended);
- `set_frequency_range`;
- `set_frequency_steps`.

They are mutually exclusive. Use whichever best fits your workflow.

#### `set_maximum_frequency(s, freq_max)`

The simplest way to define the frequency range. Numav will automatically determine a suitable set of frequency steps from 0 Hz up to `freq_max` (see [`set_frequency_sampling_density`](#set_frequency_sampling_densitys-mode)).

```julia
set_maximum_frequency(s, 200)  # Solve from 0 Hz to 200 Hz
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `freq_max` | `Real` | Upper frequency limit (Hz) |

---

#### `set_frequency_range(s, freq_min, freq_max)`

Defines both a lower and upper bound for the frequency sweep. Numav will automatically determine the frequency steps within this range, based on the sampling density.

```julia
set_frequency_range(s, 20, 200)  # Solve from 20 Hz to 200 Hz
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `freq_min` | `Real` | Lower bound of the frequency range (Hz) |
| `freq_max` | `Real` | Upper bound of the frequency range (Hz) |

---

#### `set_frequency_steps(s, frequencies)`

Manually specifies the exact set of frequencies to solve at, as a vector of values in Hz. Use this when you need full control over which frequencies are evaluated, for example to match measurement points or concentrate resolution around a resonance.

```julia
set_frequency_steps(s, [10, 40, 60, 100])  # Solve only at these four frequencies
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `frequencies` | `Vector{Real}` | List of frequencies to solve at (Hz), in any order |

> **Note:** When using `set_frequency_steps`, `set_frequency_steps_count` and `set_frequency_sampling_density` cannot be used.

---

#### `set_frequency_steps_count(s, freq_count)`

Sets the number of frequency steps the solver will evaluate within the range defined by `set_maximum_frequency` or `set_frequency_range`. If this function is never called, the default step count is **4096**. This function cannot be used if `set_frequency_steps` was used before.

```julia
set_frequency_steps_count(s, 500)  # Evaluate at 500 frequency points
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `freq_count` | `Int` | Number of frequency steps to compute |

A higher step count gives finer frequency resolution at the cost of longer computation time. The distribution of those steps within the range is controlled by [`set_frequency_sampling_density`](#set_frequency_sampling_densitys-mode).

---

#### `set_frequency_sampling_density(s, mode)`

An optional function to control how automatically-generated frequency steps are distributed within the range defined by `set_maximum_frequency` or `set_frequency_range`. If this function is never called, the **quadratic mode will be picked by default**. This function cannot be used if `set_frequency_steps` was used before.

```julia
set_frequency_sampling_density(s, FrequencySamplingDensity.constant)
set_frequency_sampling_density(s, FrequencySamplingDensity.quadratic)
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `mode` | `FrequencySamplingDensity` | Sampling strategy to use |

**Sampling modes:**

| Mode | Description |
|---|---|
| `FrequencySamplingDensity.constant` | Frequency steps are evenly spaced (uniform spacing). Suitable for broadband analyses where equal resolution at all frequencies is desired. |
| `FrequencySamplingDensity.quadratic`<br>(default) | Frequency discretization density grows quadratically with frequency. This tends to follow the modal density of rooms, being the recommended option for good accuracy and computation time balance. |

---

### `load_mesh(s, file_path)`

Loads the mesh file. It defines the geometry of the domain and must contain the physical groups that are referenced when assigning boundary conditions. Currently, only meshes in the **small field BDF** format are supported. To generate meshes, a good open-source option is [Gmsh](https://gmsh.info/). The coordinates present in the mesh file must be given in meters.

```julia
load_mesh(s, "some_mesh_file.bdf")
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `file_path` | `String` | Path to the mesh file (`.bdf`) |

---

### `add_volume_material(s; physical_group, density, sound_speed)`

Assigns acoustic material properties to a volumetric region of the mesh, identified by its physical group tag. Both `density` and `sound_speed` are [frequency-dependent quantities](#frequency-dependent-quantities).

```julia
rho(f) = 1.20 # Air density in kg/m³
c(f) = 343    # Speed of sound in m/s

add_volume_material(s, physical_group=1, density=rho, sound_speed=c)
```

| Keyword | Type | Description |
|---|---|---|
| `physical_group` | `Int` | Physical group ID from the mesh |
| `density` | frequency-dependent | Density (kg/m³) |
| `sound_speed` | frequency-dependent | Speed of sound (m/s) |

Frequency-dependent properties can model media such as absorptive foams or frequency-dependent fluids:

```julia
rho(f) = 1.20 + 0.001*f
c(f) = 340 + 0.1*sqrt(f)
add_volume_material(s, physical_group=1, density=rho, sound_speed=c)
```

---

### `add_sound_source(s; ...)`

Sources can be applied either at a specific point in space (via `coordinates`) or over a surface/volume region (via `physical_group`). Three excitation types are supported: `volume_velocity`, `particle velocity`, and `pressure`. All excitation amplitudes are [frequency-dependent quantities](#frequency-dependent-quantities).

#### Volume velocity source (monopole)

Suitable for modeling compact omnidirectional sources.

```julia
Q(f) = 10/f # Volume velocity in m³/s as a function of frequency
add_sound_source(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
```

#### Particle velocity source (vibrating surface)

Prescribes the normal component of particle velocity on all surfaces of a physical group. Useful for modeling vibrating panels or pistons.

```julia
U(f) = 15/f # Particle velocity in m/s
add_sound_source(s, physical_group=2, particle_velocity=U)
```

#### Pressure source

Prescribes acoustic pressure, either at a point or on a surface.

```julia
P(f) = 2f # Pressure in Pa as a function of frequency

# At a point in space
add_sound_source(s, coordinates=[2.0, 2.5, 1.0], pressure=P)

# On a mesh surface
add_sound_source(s, physical_group=3, pressure=P)
```

**Summary of keyword arguments:**

| Keyword | Type | Description |
|---|---|---|
| `coordinates` | `Vector{Float64}` | `[x, y, z]` location of a point source in m |
| `physical_group` | `Int` | Physical group ID of a surface or volume region |
| `volume_velocity` | frequency-dependent | Volume velocity `Q(f)` in m³/s |
| `particle_velocity` | frequency-dependent | Normal particle velocity `U(f)` in m/s |
| `pressure` | frequency-dependent | Acoustic pressure `P(f)` in Pa |

> **Note:** Each call to `add_sound_source` should specify either `coordinates` or `physical_group` (not both), and exactly one excitation type keyword.

---

### `add_surface_material(s; physical_group, impedance)`

Assigns a specific acoustic impedance boundary condition to a surface in the mesh. This is used to model absorbing or reflecting walls, liners, and other boundary treatments.

The impedance is the ratio of acoustic pressure to normal particle velocity (in Pa·s/m), and is a [frequency-dependent quantity](#frequency-dependent-quantities). It can be complex-valued to represent both resistive and reactive components.

```julia
# As a function
Z(f) = 1f + 2im
add_surface_material(s, physical_group=4, impedance=Z)

# As a tabulated file
add_surface_material(s, physical_group=4, impedance="data/wall_impedance.txt")
```

| Keyword | Type | Description |
|---|---|---|
| `physical_group` | `Int` | Physical group ID of the boundary surface |
| `impedance` | frequency-dependent | specific surface acoustic impedance (Pa·s/m) |

---

### `set_result_export_path(s, file_path)`

Specifies the file path where simulation results will be written. Results are stored in [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/) (`.h5`).

```julia
set_result_export_path(s, "result.h5")
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `file_path` | `String` | Output file path (`.h5`) |

---

### `simulate(s)`

Assembles and solves the system of equations for all frequencies in the defined range, and writes the results to the specified export path during the solution.

```julia
simulate(s)
```

#### Output format

Results are exported as HDF5 files. The file typically contains all the results as well as the inputs provided by you. To read the results it is recommended to use [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) or [HDFView](https://www.hdfgroup.org/download-hdfview/). In Julia, you can post-process results with:

```julia
using HDF5
fid = h5open("result.h5", "r")
# inspect contents
```

---
