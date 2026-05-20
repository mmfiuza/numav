# Glossary

- entity: Some geometrical body (points, lines, surfaces, volumes, etc).

# Acronyms:

- IDX (index): An unsigned integer that uniquely identifies an entity.

- EPG (external physical group): An unsigned integer present in an external mesh file (.bdf, etc) or API that classifies an entity into some group, usually a material. An EPG sequence is NOT necessarily consecutive, since mesh formats can have arbitrary conventions.
    * ESPG: External Surface Physical Group
    * EVPG: External Volume Physical Group

- IPG (internal physical group): An unsigned integer that classifies an entity into some group in Numav's internal implementation. An IPG sequence is ALWAYS consecutive starting at 0, because that is conveninent to manipulate in C++.
    - ISPG: Internal Surface Physical Group
        - ISPGI: Internal Surface Physical Group with Impedance
        - ISPGV: Internal Surface Physical Group with Velocity
        - ISPGP: Internal Surface Physical Group with Pressure
    - IVPG: Internal Volume Physical Group

- FI (frequency independent): Some quantity that will be the constant relative to sound frequency.

- FD (frequency dependent): Some quantity that will change relative to sound frequency.

- NNZ (non-zero): Some quantity that is not null. Usually refers to the non-zero terms of a sparse matrix.

- SEI (Surface Element Index):
    - ISEI (Impedance Surface Element Index)
    - VSEI (Velocity Surface Element Index)
    - PSEI (Pressure Surface Element Index)
- VEI (Volume Element Index)

- FIPI (Frequency Independent Part Index)

- GPI (Gauss Point Index)

- DI (Dimension Index)

- NI (Node Index)

- NCI (Node Combination Index)

- ENI (Elementary Node Index)

- PVNI (Point Velocity Node Index)

- PPNI (Point Pressure Node Index)

- PVI (Pressure Value Index)

- RI (Receiver Index)

- NZI (Non-Zero Index)

- PQ (Physical Quantity)
