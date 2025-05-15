using Documenter, QuantumEspressoIO
DocMeta.setdocmeta!(
    QuantumEspressoIO,
    :DocTestSetup, :(using QuantumEspressoIO);
    recursive=true,
)
doctest(
    QuantumEspressoIO,
    # fix=true,  # update all the output in `jldoctest`
)


using TestItemRunner
@run_package_tests verbose = true
