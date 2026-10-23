[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://m3g.github.io/MolSimToolkitShared.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://m3g.github.io/MolSimToolkitShared.jl/dev)
[![Tests](https://img.shields.io/badge/build-passing-green)](https://github.com/m3g/MolSimToolkitShared.jl/actions)
[![codecov](https://codecov.io/gh/m3g/MolSimToolkitShared.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/m3g/MolSimToolkitShared.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

# MolSimToolkitShared

This package is a placeholder for functions that are shared
among other packages. The functions are now defined in
[PDBToolsBase.jl](https://github.com/m3g/PDBToolsBase.jl), and imported
from there, such that packages that depend on MolSimToolkitShared and
packages that depend on PDBToolsBase extend the same functions. New code
should depend on PDBToolsBase directly.

See [the documentation](https://m3g.github.io/MolSimToolkitShared.jl)

