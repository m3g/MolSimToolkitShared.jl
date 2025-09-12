module MolSimToolkitShared

using Compat: @compat
using TestItems: @testitem

using LinearAlgebra: eigvecs, norm, cross, dot
using StaticArrays: SVector, MMatrix, SMatrix

# Currently only shared function names
function distance end
function distances end
function bulk_coordination end

function coordination_number(::String, args...; kargs...)
    throw(ArgumentError("""\n
        Invalid arguments for the `coordination_number` function.
        Plese check the documentation for the correct call signature, by typing:

        julia> ? coordination_number

    """))
end
@testitem "coordination_number" begin
    import MolSimToolkitShared: coordination_number
    # The coordination_number(::String, args...; kargs...) is a placeholder for the docs only
    @test_throws ArgumentError coordination_number("string.txt", 1.0)
end

# Utility functions
include("./wrap.jl")
include("./structural_alignment.jl")
include("./dihedral.jl")

# Public API
@compat public distance, distances
@compat public coordination_number, bulk_coordination
@compat public center_of_mass
@compat public wrap, wrap_to_first
@compat public align, align!, rmsd
@compat public dihedral, dihedrals

end
