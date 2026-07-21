# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    Fem,
    Helmholtz,
    Tetrahedron,
    Constant,
    Linear,
    Quadratic

abstract type Option end

struct FemType <: Option end
const Fem = FemType()

struct HelmholtzType <: Option end
const Helmholtz = HelmholtzType()

struct TetrahedronType <: Option end
const Tetrahedron = TetrahedronType()

struct ConstantType <: Option end
const Constant = ConstantType()

struct LinearType <: Option end
const Linear = LinearType()

struct QuadraticType <: Option end
const Quadratic = QuadraticType()
