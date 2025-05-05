using TestItemRunner

# Temporarily loading artifacts for all tests here
using LazyArtifacts

@run_package_tests verbose = true
