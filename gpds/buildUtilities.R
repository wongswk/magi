cat("test finished, now running test coverage report...\n")
covr::report(covr::package_coverage("."))
cat("test coverage report generated.")
