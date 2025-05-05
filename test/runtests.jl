using TestItemRunner

@testitem "doctest" begin
    # Test `jldoctest` in docstring
    using Documenter, QuantumEspressoIO

    doctest(
        QuantumEspressoIO,
        fix=true,  # update all the output in `jldoctest`
    )
end


@run_package_tests verbose = true
