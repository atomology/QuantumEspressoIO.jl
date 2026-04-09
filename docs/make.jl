using QuantumEspressoIO
using Documenter

DocMeta.setdocmeta!(
    QuantumEspressoIO,
    :DocTestSetup, :(using QuantumEspressoIO);
    recursive=true,
)

makedocs(;
    sitename="QuantumEspressoIO.jl",
    authors="Aleksandr Poliukhin, Junfeng Qiao, Guoyuan Liu",
    modules=[QuantumEspressoIO],
    pages=[
        "Home" => "index.md",
        "API" => [
            "`pw.x`" => "api/pw.md",
            "`band.x`" => "api/band.md",
            "`projwfc.x`" => "api/projwfc.md",
            "`open_grid.x`" => "api/opengrid.md",
            "Fortran Namelist" => "api/namelist.md",
            "XML" => "api/xml.md",
            "Binaries" => "api/bin.md",
            "Utilities" => "api/utils.md",
        ]
    ],
    # doctest=:fix,  # update all the jldoctest
)

# Documenter will auto detect build environment; on local machine it will be
# skipped, so it's safe to run this script
deploydocs(; repo="github.com/atomology/QuantumEspressoIO.jl.git", devbranch="main")
