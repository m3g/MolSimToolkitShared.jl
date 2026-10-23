"""
    MolSimToolkitShared

The shared function names and utilities of the ecosystem (MolSimToolkit.jl, PDBTools.jl,
ComplexMixtures.jl, Packmol.jl) are now defined in PDBToolsBase.jl. This package imports them
from there, such that packages that depend on MolSimToolkitShared and packages that depend on
PDBToolsBase extend the same functions. New code should depend on PDBToolsBase directly.
"""
module MolSimToolkitShared

using TestItems: @testitem

# The shared functions are those of PDBToolsBase
import PDBToolsBase: distance, distances
import PDBToolsBase: positions
import PDBToolsBase: coordination_number
import PDBToolsBase: center_of_mass
import PDBToolsBase: wrap, wrap_to_first
import PDBToolsBase: align, align!, rmsd, alignment_movements, apply_alignment_transformation!
import PDBToolsBase: dihedral, dihedrals
import PDBToolsBase: get_atoms
import PDBToolsBase: save, load

# Not part of PDBToolsBase: only MolSimToolkit.jl defines methods for it
function bulk_coordination end

# Public API
public distance, distances
public positions
public coordination_number, bulk_coordination
public center_of_mass
public wrap, wrap_to_first
public align, align!, rmsd, alignment_movements, apply_alignment_transformation!
public dihedral, dihedrals
public get_atoms
public save
public load

@testitem "shared names" begin
    import PDBToolsBase
    for name in (
        :distance, :distances, :positions, :coordination_number, :center_of_mass,
        :wrap, :wrap_to_first, :align, :align!, :rmsd, :alignment_movements,
        :apply_alignment_transformation!, :dihedral, :dihedrals, :get_atoms, :save, :load,
    )
        @test Base.ispublic(MolSimToolkitShared, name)
        # the same function object as in PDBToolsBase
        @test getfield(MolSimToolkitShared, name) === getfield(PDBToolsBase, name)
    end
    @test Base.ispublic(MolSimToolkitShared, :bulk_coordination)
end

@testitem "coordination_number" begin
    import MolSimToolkitShared: coordination_number
    # The coordination_number(args...; kargs...) is a placeholder for the docs only
    @test_throws ArgumentError coordination_number("string.txt", 1.0)
end

end
