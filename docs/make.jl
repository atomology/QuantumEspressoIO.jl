using QuantumEspressoIO
using Documenter

DocMeta.setdocmeta!(QuantumEspressoIO, :DocTestSetup, :(using QuantumEspressoIO); recursive=true)

makedocs(;
    modules=[QuantumEspressoIO],
    authors="Junfeng Qiao, Guoyuan Liu, Aleksandr Poliukhin",
    sitename="QuantumEspressoIO.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
