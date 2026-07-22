# Numav.jl

Numav.jl is a Julia wrapper for the C++ library Numav, to perform acoustics and vibrations simulations. Currently, it only supports frequency domain FEM simulations for acoustics. However, it is a plan to cover more numerical methods soon.

---

## Installation
Just type into the Julia REPL:
```julia
]add Numav
```

---

## Quick start

Here is a complete example on how to run your first simulation in Numav:

```julia
using Numav

# Create the simulation object
s = create_simulation(
    numerical_method = Fem,
    equation = Helmholtz,
    element_shape = Tetrahedron,
    element_order = Linear
)

# Set simulated frequencies
set_frequency!(s, max=200) # (Hz)

# Load mesh
load_mesh!(s, "path/to/mesh_file.bdf")

# Define material (frequency-dependent)
rho(f) = 1.20
c(f) = 343
add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)

# Add a monopole volume velocity source at a point
Q(f) = 10/f
add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)

# Add a particle velocity source on a surface
U(f) = 15/f
add_sound_source!(s, physical_group=2, particle_velocity=U)

# Add pressure sources (to point or surface)
P(f) = 2f
add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)
add_sound_source!(s, physical_group=3, pressure=P)

# Add a surface impedance boundary condition
Z(f) = 1f + 2im
add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)

# Set output file
set_result_export_path!(s, "result.h5")

# Run the simulation
run!(s)
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
add_volume_material!(s, physical_group=1, density=1.20, speed_of_sound=343)
```

### 3. Tabulated file (path string)

A `String` containing the path to a plain-text CSV-like file with three columns: `frequency`, `real`, and `imag`. Numav reads the table and interpolates between the provided data points at each required frequency step. There is no header row.

```
  0, 1, 2
 25, 1, 2
 50, 5, 1
 75, 3, 2
100, 1, 0
```

```julia
add_surface_material!(s, physical_group=4, specific_acoustic_impedance="data/wall_impedance.txt")
```

The file format rules are:
- Three columns separated by commas: `frequency` (Hz), `real` part, `imaginary` part.
- No header row.
- For purely real quantities (e.g. density), set the `imag` column to `0`.

> **Note:** Any table should contain data covering all the frequency range that you plan to simulate.
---

## API reference

### `create_simulation(...)`

Creates an instance of type `Simulation`.

```julia
s = create_simulation(
    numerical_method = Fem,
    equation = Helmholtz,
    element_shape = Tetrahedron,
    element_order = Linear
)
```

**Type parameters:**

| Parameter | Supported options | Description |
|---|---|---|
| `numerical_method` | `Fem` | Numerical method used |
| `equation` | `Helmholtz` | Differential equation to be solved |
| `element_shape` | `Tetrahedron` | Geometrical shape of elements |
| `element_order` | `Linear`, `Quadratic` | Polynomial order of the finite elements (1st or 2nd order) |

---
 
### `set_frequency!(s; max, min, length, sampling_density, step, vector)`

Sets the frequecies to perform the simulations with two mutually exclusive ways to use it:
- **automatic mode**: pass `max` (and optionally `min`, `length`, `sampling_density` and `step`), and Numav will automatically determine a suitable set of frequency steps;
- **manual mode**: pass `vector` directly, giving full control over which frequencies are evaluated.
---

```julia
set_frequency!(s, max=200)  # Solve from 0 Hz to 200 Hz
```
 
| Keyword | Type | Description |
|---|---|---|
| `max` | `Real` | Upper frequency limit (Hz). Required unless `steps` is provided; mutually exclusive with `steps`. |
| `min` | `Real` | Lower bound of the frequency range (Hz). Only usable together with `max`. Defaults to `0` if omitted. |
| `length` | `Integer` | Number of frequency steps to compute within the `min`/`max` range. Only usable together with `max`. Defaults to `4096` if omitted. |
| `sampling_density` | `Numav.Option` | Sampling strategy to use within the the range: `Constant` or `Quadratic`. Only usable together with `max`. Defaults to `Quadratic` if omitted. |
| `step` | `Real` | Defference in Hertz between each frequency step, supposing equally spaced steps. Only usable together with `max`. Cannot be used together with `length` and `sampling_density` |
| `vector` | `Vector{Real}` | List of frequencies in Hertz to solve at. Cannot be used together with `max`, `min`, `length`, `sampling_density` and `step`. |
 
**Solve up to a maximum frequency**, letting Numav automatically pick the steps from 0 Hz:
 
```julia
set_frequency!(s, max=200) # Solve from 0 Hz to 200 Hz
```
 
**Solve within an explicit range**, by adding `min`:
 
```julia
set_frequency!(s, min=20, max=200) # Solve from 20 Hz to 200 Hz
```
 
**Control the number of frequency steps** with `length` (defaults to `4096` when not specified):
 
```julia
set_frequency!(s, max=200, length=500) # Evaluate at 500 frequency points
```
 
A higher step count gives finer frequency resolution at the cost of longer computation time.
 
**Control how steps are distributed** with `sampling_density` (defaults to `Quadratic` when not specified):
 
```julia
set_frequency!(s, max=200, sampling_density=Constant)
set_frequency!(s, max=200, sampling_density=Quadratic)
```
 
**Sampling modes:**
 
| Mode | Description |
|---|---|
| `Constant` | Frequency steps are evenly spaced (uniform spacing). Suitable for broadband analyses where equal resolution at all frequencies is desired. |
| `Quadratic`<br>(default) | Frequency discretization density grows quadratically with frequency. This tends to follow the modal density of rooms, being the recommended option for good accuracy and computation time balance. |
 
**Manually specify the exact frequencies** with `vector`, for example to match measurement points or concentrate resolution around a resonance:
 
```julia
set_frequency!(s, vector=[10, 40, 60, 100]) # Solve only at these frequencies
```
 
---

### `load_mesh!(s, file_path)`

Loads the mesh file. It defines the geometry of the domain and must contain the physical groups that are referenced when assigning boundary conditions. Currently, only meshes in the **small field BDF** format are supported. To generate meshes, a good open-source option is [Gmsh](https://gmsh.info/). The coordinates present in the mesh file must be given in meters.

```julia
load_mesh!(s, "some_mesh_file.bdf")
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `file_path` | `String` | Path to the mesh file (`.bdf`) |

---

### `add_volume_material!(s; physical_group, density, speed_of_sound)`

Assigns acoustic material properties to a volumetric region of the mesh, identified by its physical group tag. Both `density` and `speed_of_sound` are [frequency-dependent quantities](#frequency-dependent-quantities).

```julia
rho(f) = 1.20 # air density in kg/m³
c(f) = 343 # speed of sound in m/s
add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
```

| Keyword | Type | Description |
|---|---|---|
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID (or vector of IDs) from the mesh |
| `density` | frequency-dependent | Density in kg/m³ |
| `speed_of_sound` | frequency-dependent | Speed of sound in m/s |

Frequency-dependent properties can model any media, including frequency-dependent complex quantities to add absorption:

```julia
rho(f) = 1.20 + 0.001*f
c(f) = 340 + 0.1*sqrt(f)
add_volume_material!(s, physical_group=1, density=rho, speed_of_sound=c)
```

Additionally, multiple physical groups can be assigned at once:
```julia
add_volume_material!(s, physical_group=[4,6], density=1.2, speed_of_sound=343)
```

---

### `add_sound_source!(s; ...)`

Sources can be applied either at a specific point in space (via `coordinates`) or over a surface/volume region (via `physical_group`). Three excitation types are supported: `volume_velocity`, `particle velocity`, and `pressure`. All excitation amplitudes are [frequency-dependent quantities](#frequency-dependent-quantities).

#### Volume velocity source (monopole)

Suitable for modeling punctual omnidirectional sources.

```julia
Q(f) = 10/f # Volume velocity in m³/s as a function of frequency
add_sound_source!(s, coordinates=[1.0, 1.5, 2.0], volume_velocity=Q)
```

#### Particle velocity source (vibrating surface)

Prescribes the normal component of particle velocity on all surfaces of a physical group. Useful for modeling vibrating panels or pistons.

```julia
U(f) = 15/f # Particle velocity in m/s
add_sound_source!(s, physical_group=2, particle_velocity=U)
```

#### Pressure source

Prescribes acoustic pressure, either at a point or on a surface.

```julia
P(f) = 2f # Pressure in Pa as a function of frequency

# At a point in space
add_sound_source!(s, coordinates=[2.0, 2.5, 1.0], pressure=P)

# On a mesh surface
add_sound_source!(s, physical_group=3, pressure=P)
```

| Keyword | Type | Description |
|---|---|---|
| `coordinates` | `Vector{Real}`, `Vector{Vector{Real}}` | `[x, y, z]` location of a point source in m |
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of a surface or volume region |
| `volume_velocity` | frequency-dependent | Volume velocity in m³/s |
| `particle_velocity` | frequency-dependent | Normal particle velocity in m/s |
| `pressure` | frequency-dependent | Acoustic pressure in Pa |

> **Note:** Each call to `add_sound_source!` should specify either `coordinates` or `physical_group` (not both), and exactly one excitation type keyword.

Additionally, multiple physical groups or points can be assigned at once:
```julia
add_sound_source!(s, coordinates=[[1.0, 3.0, 2.0], [1.0, 1,0, 1.0]], volume_velocity=0.03)
add_sound_source!(s, physical_group=[3,5,9,2], density=1.2, particle_velocity=0.01)
```

---

### `add_surface_material!(s; physical_group, specific_acoustic_impedance)`

Assigns a specific acoustic impedance boundary condition to a surface in the mesh. This is used to model absorbing or reflecting walls, liners, and other boundary treatments.

The impedance is the ratio of acoustic pressure to normal particle velocity (in Pa·s/m), and is a [frequency-dependent quantity](#frequency-dependent-quantities). It can be complex-valued.

```julia
Z(f) = 1f + 2im
add_surface_material!(s, physical_group=4, specific_acoustic_impedance=Z)
```

| Keyword | Type | Description |
|---|---|---|
| `physical_group` | `Integer`, `Vector{Integer}` | Physical group ID of the boundary surface |
| `specific_acoustic_impedance` | frequency-dependent | specific surface acoustic impedance (Pa·s/m) |

Additionally, multiple physical groups can be assigned at once:
```julia
add_surface_material!(s, physical_group=[3,1], specific_acoustic_impedance=1.0)
```

---

### `set_result_export_path!(s, file_path)`

Specifies the file path where simulation results will be written. Results are stored in [HDF5 format](https://www.hdfgroup.org/solutions/hdf5/) (`.h5`).

```julia
set_result_export_path!(s, "result.h5")
```

| Argument | Type | Description |
|---|---|---|
| `s` | `Simulation` | The simulation object |
| `file_path` | `String` | Output file path (`.h5`) |

---

### `run!(s)`

Assembles and solves the system of equations for all frequencies in the defined range, and writes the results to the specified export path during the solution.

```julia
run!(s)
```

#### Output format

Results are exported as HDF5 files. The file typically contains all the results as well as the inputs provided by you. To read the results it is recommended to use [HDF5.jl](https://juliaio.github.io/HDF5.jl/stable/) or [HDFView](https://www.hdfgroup.org/download-hdfview/). In Julia, you can post-process results with:

```julia
using HDF5
fid = h5open("result.h5", "r")
# inspect contents
```

---
