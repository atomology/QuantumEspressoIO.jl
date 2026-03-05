# QuantumEspressoIO.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://io.wannierjl.org/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://io.wannierjl.org/dev)
[![CI](https://github.com/atomology/QuantumEspressoIO.jl/workflows/CI/badge.svg)](https://github.com/atomology/QuantumEspressoIO.jl/actions/workflows/CI)
[![codecov](https://codecov.io/gh/atomology/QuantumEspressoIO.jl/graph/badge.svg?token=L3I63QL7IE)](https://codecov.io/gh/atomology/QuantumEspressoIO.jl)

A Julia package for reading/writing [Quantum ESPRESSO](https://gitlab.com/QEF/q-e) files.

## Quick examples

```julia
using QuantumEspressoIO

julia> pw = read_pw_in("pw.in")

julia> pw[:atomic_positions][:atoms]
2-element Vector{String}:
 "Si"
 "Si"

julia> write_pw_in("pw2.in", pw)
```
