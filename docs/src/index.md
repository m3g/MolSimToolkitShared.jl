```@meta
CollapsedDocStrings = true
```

# MolSimToolkitShared.jl

[MolSimToolkitShared.jl](https://github.com/m3g/MolSimToolkitShared.jl) is a 
low-level package defining a series of functions and function names that 
are shared among other packages, such as [MolSimToolkit.jl](https://github.com/m3g/MolSimToolkit.jl), 
[PDBTools.jl](https://github.com/m3g/PDBTools.jl),
and [ComplexMixtures.jl](https://github.com/m3g/ComplexMixtures.jl).

Normally this package won't be used by the end-user, but as a dependency of other packages.


## API

### Function name placeholders

These function names are considered public:

```
distance, distances
coordination_number, bulk_coordination
center_of_mass
wrap, wrap_to_first
align, align!, rmsd
```

### Methods

```@autodocs
Modules = [ MolSimToolkitShared ]
```
