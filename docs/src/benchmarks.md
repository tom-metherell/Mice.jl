# Benchmarks

I have (very much not rigorously) benchmarked `Mice.jl` using the [test dataset](https://archive.ics.uci.edu/dataset/878) [dickson_prognosis_1989](@cite), and also performed an equivalent benchmark of the R package `mice`.

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric). Both used predictive mean matching for all variables that were to be imputed. In `Mice.jl`, `gcSchedule` was set to `0.3`.

## Benchmark results
### Linux
System info: Single-threaded execution, Intel® Core™ i7-12700H 2.30GHz CPU, 32GB 4800MHz DDR5 RAM, running Ubuntu (WSL) version 22.04.2.

R: version 4.3.2 running `mice` version 3.14.0.
Julia: version 1.9.4 running `Mice.jl` version 0.3.0.

| Number of imputations | R (`mice`) (s) | `Mice.jl` (s) |
| --- | --- | --- |
| 1 | 1.22 | 15.21 |
| 5 | 5.96 | 17.24 |
| 10 | 11.74 | 18.46 |
| 20 | 26.32 | 21.11 |
| 50 | 67.26 | 27.24 |
| 100 | 136.58 | 38.99 |

## Why is `Mice.jl` so slow for small jobs?

Julia is a compiled language. This means that the first time a function is run, it is compiled into machine code, which takes time. Therefore, the first iteration of `mice()` will be (much) slower in Julia than in R, for example. However, subsequent iterations will be much faster, as all of the required functions are already compiled.

## Why is the first iteration so much slower than the rest?

See above.