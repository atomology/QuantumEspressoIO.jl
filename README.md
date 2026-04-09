# QuantumEspressoIO.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://atomology.github.io/QuantumEspressoIO.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://atomology.github.io/QuantumEspressoIO.jl/dev)
[![CI](https://github.com/atomology/QuantumEspressoIO.jl/workflows/CI/badge.svg)](https://github.com/atomology/QuantumEspressoIO.jl/actions/workflows/CI)
[![codecov](https://codecov.io/gh/atomology/QuantumEspressoIO.jl/graph/badge.svg?token=L3I63QL7IE)](https://codecov.io/gh/atomology/QuantumEspressoIO.jl)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)

A Julia package for reading/writing [Quantum ESPRESSO](https://gitlab.com/QEF/q-e) files.

## Quick examples

```julia
using QuantumEspressoIO

julia> pw = read_pw_in("pw.in")

julia> pw["atomic_positions"]["atoms"]
2-element Vector{String}:
 "Si"
 "Si"

julia> write_pw_in("pw2.in", pw)
```
