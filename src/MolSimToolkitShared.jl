module MolSimToolkitShared

import PrecompileTools
using Compat: @compat
using TestItems: @testitem

using LinearAlgebra: eigvecs, norm, cross, dot
using StaticArrays: SVector, MMatrix, SMatrix

# Currently only shared function names
function distance end
function distances end
function bulk_coordination end
function positions end
function get_atoms end

@testitem "shared names" begin
    if VERSION >= v"1.11"
        @test Base.ispublic(MolSimToolkitShared, :distance)
        @test Base.ispublic(MolSimToolkitShared, :distances)
        @test Base.ispublic(MolSimToolkitShared, :bulk_coordination)
        @test Base.ispublic(MolSimToolkitShared, :positions)
        @test Base.ispublic(MolSimToolkitShared, :get_atoms)
    end
end

function coordination_number(args...; kargs...)
    throw(ArgumentError("""\n
        Invalid arguments for the `coordination_number` function.
        Please check the documentation for the correct call signature, by typing:

        julia> ? coordination_number

    """))
end
@testitem "coordination_number" begin
    import MolSimToolkitShared: coordination_number
    # The coordination_number(args...; kargs...) is a placeholder for the docs only
    @test_throws ArgumentError coordination_number("string.txt", 1.0)
end

# Utility functions
include("./wrap.jl")
include("./structural_alignment.jl")
include("./dihedral.jl")

# Public API
@compat public distance, distances
@compat public positions
@compat public coordination_number, bulk_coordination
@compat public center_of_mass
@compat public wrap, wrap_to_first
@compat public align, align!, rmsd, alignment_movements, apply_alignment_transformation!
@compat public dihedral, dihedrals
@compat public get_atoms

# Precompilation tools
include("precompile.jl")

end
