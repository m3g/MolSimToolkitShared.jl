module MolSimToolkitShared

using Compat: @compat
using TestItems: @testitem

using LinearAlgebra: eigvecs
using StaticArrays: SVector, MMatrix, SMatrix

# Currently only shared function names
function distance end
function distances end
function coordination_number end
function bulk_coordination end

# Utility functions
include("./wrap.jl")
include("./structural_alignment.jl")

# Public API
@compat public distance, distances
@compat public coordination_number, bulk_coordination
@compat public center_of_mass
@compat public wrap, wrap_to_first
@compat public align, align!, rmsd

end
