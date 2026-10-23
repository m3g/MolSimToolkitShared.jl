MolSimToolkitShared.jl Changelog
===========================
  
[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-experimental]: https://img.shields.io/badge/Experimental-yellow.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-fix]: https://img.shields.io/badge/Fix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 1.7.0-DEV
-------------
- ![INFO][badge-info] The shared functions (`distance`, `distances`, `positions`, `coordination_number`, `center_of_mass`, `wrap`, `wrap_to_first`, `align`, `align!`, `rmsd`, `alignment_movements`, `apply_alignment_transformation!`, `dihedral`, `dihedrals`, `get_atoms`, `save`, `load`) are now imported from [PDBToolsBase.jl](https://github.com/m3g/PDBToolsBase.jl), such that packages that depend on MolSimToolkitShared and packages that depend on PDBToolsBase extend the same functions. `bulk_coordination` remains defined here. New code should depend on PDBToolsBase directly.
- ![INFO][badge-info] Requires Julia 1.11 or later (as PDBToolsBase).

Version 1.6.0
-------------
- ![FEATURE][badge-feature] add `save` and `load` to public shared function names.

Version 1.5.1
-------------
- ![FIX][badge-fix] `coordination_number` with string first argument does not throw anymore, instead the fully non-specialized call throws the argument error.

Version 1.5.0
-------------
- ![FEATURE][badge-feature] add `get_atoms` as a public shared function.

Version 1.4.0
-------------
- ![FEATURE][badge-feature] add `positions` as a public shared function.

Version 1.3.0
-------------
- ![FEATURE][badge-feature] split alignment functions to expose `alignment_movements` and `apply_alignment_transformations!` functions.

Version 1.2.1
-------------
- ![ENHANCEMENT][badge-enhancement] add precompilation statements.
- ![INFO][badge-info] Update CI action files.

Version 1.2.0
-------------
- ![FEATURE][badge-feature] add `dihedral` and `dihedrals` functions.
- ![INFO][badge-info] add CHANGELOG.md file and CI run.
