The scripts in this folder replicate the results of the benchmark comparison (Section 5).

The main script is "run-all.R" which will run 100 repetitions of the three methods (MAGI, deBInfer, and CollocInfer), and then compute the numerical summaries in Table 1. For simplicity the script runs everything sequentially on a single CPU core, but the repetitions can be run in parallel to reduce the time needed to complete the benchmark. The outputs of each method will be saved in the "results" folder.
