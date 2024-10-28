# abstract type for all boundary conditions
abstract type AbstractBoundaryConditions end

struct PeriodicBoundaryConditions <: AbstractBoundaryConditions end