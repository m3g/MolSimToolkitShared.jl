```@meta
CollapsedDocStrings = true
```

# MolSimToolkitShared.jl

[MolSimToolkitShared.jl](https://github.com/m3g/MolSimToolkitShared.jl) is a 
low-level package with the function names that are shared among other packages, such as
[MolSimToolkit.jl](https://github.com/m3g/MolSimToolkit.jl), 
[PDBTools.jl](https://github.com/m3g/PDBTools.jl),
[ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl) and
[Packmol.jl](https://github.com/m3g/Packmol.jl).

!!! note "Superseded by PDBToolsBase.jl"
    The shared functions are now defined in
    [PDBToolsBase.jl](https://github.com/m3g/PDBToolsBase.jl), and MolSimToolkitShared imports
    them from there: packages that depend on MolSimToolkitShared and packages that depend on
    PDBToolsBase extend the same functions. New code should depend on PDBToolsBase directly.

Normally this package won't be used by the end-user, but as a dependency of other packages.

```@docs
MolSimToolkitShared
```

## API

These function names are public:

```
distance, distances
positions
coordination_number, bulk_coordination
center_of_mass
wrap, wrap_to_first
align, align!, rmsd, alignment_movements, apply_alignment_transformation!
dihedral, dihedrals
get_atoms
save, load
```

All of them, except `bulk_coordination`, are the functions of PDBToolsBase, and are documented
in the [PDBToolsBase.jl documentation](https://m3g.github.io/PDBToolsBase.jl/dev/api/).
