# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

export
    Fem,
    Helmholtz,
    Tetrahedron,
    Constant,
    Linear,
    Quadratic

abstract type Option end

struct FemOptionType <: Option end
const Fem = FemOptionType()

struct HelmholtzOptionType <: Option end
const Helmholtz = HelmholtzOptionType()

struct TetrahedronOptionType <: Option end
const Tetrahedron = TetrahedronOptionType()

struct ConstantOptionType <: Option end
const Constant = ConstantOptionType()

struct LinearOptionType <: Option end
const Linear = LinearOptionType()

struct QuadraticOptionType <: Option end
const Quadratic = QuadraticOptionType()
