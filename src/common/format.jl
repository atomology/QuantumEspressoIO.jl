abstract type FileFormat end

"""Fortran formatted IO."""
struct FortranText <: FileFormat end

"""Fortran unformatted IO."""
struct FortranBinary <: FileFormat end

"""
Fortran unformatted IO with stream access.

For example, file written using these Fortran code:
```fortran
OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM", FORM="UNFORMATTED")
```
"""
struct FortranBinaryStream <: FileFormat end
