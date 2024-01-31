# Benchmarks

I have (very much not rigorously) benchmarked `Mice.jl` using the [test dataset](https://archive.ics.uci.edu/dataset/878) [dickson_prognosis_1989](@cite), and also performed an equivalent benchmark of the R package `mice`.

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric). Both used predictive mean matching for all variables that were to be imputed. In `Mice.jl`, `gcSchedule` was set to `0.3`.

## Benchmark results
### Windows
System info: Single-threaded execution, Intel® Core™ i7-12700H 2.30GHz CPU, 32GB 4800MHz DDR5 RAM, running Windows 11 version 10.0.22621.

R: version 4.3.2 running `mice` version 3.16.0.
Julia: version 1.10.0 running `Mice.jl` version 0.3.2.

| Number of imputations | R (`mice`) (s) | `Mice.jl` (s) |
| --- | --- | --- |
| 1 | 1.79 | 4.86 |
| 5 | 8.45 | 5.54 |
| 10 | 16.59 | 6.55 |
| 20 | 33.19 | 8.09 |
| 50 | 85.79 | 12.17 |
| 100 | 171.93 | 19.62 |

### Linux
System info: Single-threaded execution, Intel® Core™ i7-12700H 2.30GHz CPU, 32GB 4800MHz DDR5 RAM, running Ubuntu (WSL) version 22.04.3.

R: version 4.3.2 running `mice` version 3.16.0.
Julia: version 1.10.0 running `Mice.jl` version 0.3.2.

| Number of imputations | R (`mice`) (s) | `Mice.jl` (s) |
| --- | --- | --- |
| 1 | 1.24 | 4.93 |
| 5 | 5.92 | 5.63 |
| 10 | 11.74 | 6.56 |
| 20 | 23.87 | 8.44 |
| 50 | 63.83 | 12.54 |
| 100 | 125.65 | 21.30 |

## Why is `Mice.jl` so slow for small jobs?

Julia is a compiled language. This means that the first time a function is run, it is compiled into machine code, which takes time. Therefore, the first iteration of `mice()` will be (much) slower in Julia than in R, for example. However, subsequent iterations will be much faster, as all of the required functions are already compiled.

## Why is the first iteration so much slower than the rest?

See above.

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```