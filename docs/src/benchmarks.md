# Benchmarks

I have (very much not rigorously) benchmarked `Mice.jl` using the [test dataset](https://archive.ics.uci.edu/dataset/878) [dickson_prognosis_1989](@cite), and also performed an equivalent benchmark of the R package `mice`.

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric). Both used predictive mean matching for all variables that were to be imputed. For the Julia implementation, a new terminal was used for each test to count the time taken to compile the functions.

## Benchmark results
System info: Single-threaded execution, Apple M4 10-core CPU, 24GB LPDDR5 memory, running macOS Sequoia version 15.3.

R: version 4.4.2 running `mice` version 3.17.0.
Julia: version 1.11.3 running `Mice.jl` version 0.3.7.

| Number of imputations | R (`mice`) (s) | `Mice.jl` (s) |
| --- | --- | --- |
| 1 | 0.82 | 2.92 |
| 5 | 4.07 | 3.33 |
| 10 | 8.02 | 3.81 |
| 20 | 16.59 | 4.80 |
| 50 | 41.66 | 7.67 |
| 100 | 83.80 | 12.38 |
| 500 | not tested | 54.47 |

## Why is `Mice.jl` so slow for small jobs?

Julia is a compiled language. This means that the first time a function is run, it is compiled into machine code, which takes time. Therefore, the first iteration of `mice()` will be slower in Julia than in R, for example. However, subsequent iterations will be much faster, as all of the required functions are already compiled.

## Why is the first iteration so much slower than the rest?

See above.

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```