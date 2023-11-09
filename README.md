# Mice.jl

## What is Mice.jl?

`Mice.jl` is a Julia package for multiple imputation using chained equations. It is heavily based on the R package [mice](https://cran.r-project.org/web/packages/mice/index.html) by Stef van Buuren and Karin Groothuis-Oudshoorn [[1]](#1).

`Mice.jl`'s syntax is very similar to that of `mice` and it is intended to be used in conjunction with [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl).

## Installation

`Mice.jl` is not registered yet. To install the latest version:

```
] add https://github.com/tom-metherell/Mice.jl.git
```

or

```
using Pkg; Pkg.add("https://github.com/tom-metherell/Mice.jl.git")
```

## Quick-start guide

### Imputation (`mice()`)
Use the `mice()` function to perform multiple imputation on a DataFrame. The output is a multiply imputed dataset (`Mids`) object.

#### Usage
```
mice(df, m = 5, visitSequence = nothing, methods = nothing, predictorMatrix = nothing, iter = 10, progressReports = true, gcSchedule = 1.0, threads = true, ...)
```
where:

`df` is the DataFrame.

`m` is the number of imputations.

`visitSequence` is a vector of column names specifying the order in which the columns should be imputed. If `nothing`, the order is determined automatically. You can make a default visit sequence (which you can then edit) using the function `makeMonotoneSequence(data)`

`methods` is a vector of imputation methods. If `nothing`, the default methods are used. Currently, `Mice.jl` supports predictive mean matching (`pmm`) only. To skip a column, define the method for that column as an empty string (`""`). You can make a default methods vector (which you can then edit) using the function `makeMethods(data)`

`predictorMatrix` is a matrix of predictors for each column (with the predictors in the columns of the matrix). If `nothing`, all other columns are used for each column. To prevent one column from predicting another, ensure that value of the corresponding cell is set to `0`. You can make a default predictor matrix (which you can then edit) using the function `makePredictorMatrix(data)`

`iter` is the number of iterations to perform.

`progressReports` is a boolean indicating whether to print progress reports.

`gcSchedule` determines how often the garbage collector is invoked (additionally to when it would be anyway). This can improve performance for large datasets. The value is the fraction of your machine's RAM that remains free at which the garbage collector will be invoked. The default is `1.0` (i.e. after every iteration of every variable), but this may be excessive for smaller jobs.

`threads` determines whether the imputations are executed in parallel using multithreading. The default is `true`, but this may make performance worse for small jobs (see [Benchmarks](#benchmarks)). Note that unlike in R, [random number generation is thread-safe by default in Julia](https://julialang.org/blog/2021/11/julia-1.7-highlights/#new_rng_reproducible_rng_in_tasks).

#### Example
```
using CSV, DataFrames, Mice, Random

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df)
```

After imputation, you can use `plot(mids, variableNumber)` or `plot(mids, variableName)` to inspect the mean and variance trace plots.

### Repeated analysis (`with()`)
Use the `with()` function to perform data analysis functions on a `Mids` object. The output is a multiply imputed repeated analyses (`Mira`) object.

#### Usage
```
with(mids, function)
```
where `mids` is a `Mids` object.

`function` is a data analysis function. It should take the form `data -> function(arguments, data, moreArguments...)`, where the placeholder `data` (which can take other names if necessary) signifies the position of the data argument in the function.

#### Example
```
using CSV, DataFrames, GLM, Mice, Random

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df);

analyses = with(imputedData, data -> glm(@formula(y ~ x1 + x2), data, Poisson(), LogLink()))
```

### Pooling coefficients (`pool()`)
Use the `pool()` function to pool the parameter estimates from a `Mira` object. The output is a multiply imputed pooled outcomes (`Mipo`) object.

#### Usage
```
pool(mira)
```
where `mira` is a `Mira` object.

#### Example
```
using CSV, DataFrames, GLM, Mice, Random

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df);

analyses = with(imputedData, data -> glm(@formula(y ~ x1 + x2), data, Poisson(), LogLink()));

results = pool(analyses)
```

## Further information

For more information about multiple imputation by chained equations, and how to solve common problems, see [Flexible Imputation of Missing Data](https://stefvanbuuren.name/fimd/) by Stef van Buuren [[2]](#2).

## Benchmarks

I have (very much not rigorously) benchmarked `Mice.jl` using the [test dataset](https://archive.ics.uci.edu/dataset/878) [[3]](#3), and also performed an equivalent benchmark of the R package `mice`.

System info: Single-threaded execution, Intel® Core™ i7-12700H 2.30GHz CPU, 32GB 4800MHz DDR5 RAM, running Windows 11 version 10.0.22621.

R: version 4.3.1 running `mice` version 3.16.0.
Julia: version 1.9.2 running `Mice.jl` version 0.0.0.

### Imputation (`mice`)

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric).
In `Mice.jl`, `gcSchedule` was set to `0.3`.

| Number of imputations | R (`mice`) (s) | `Mice.jl` (s) |
| --- | --- | --- |
| 1 | 1.88 | 27.01 |
| 5 | 8.84 | 28.72 |
| 10 | 18.14 | 30.81 |
| 20 | 38.42 | 36.34 |
| 50 | 96.78 | 51.95 |
| 100 | 199.33 | 79.52 |

### Why is `Mice.jl` so slow for small jobs?

Julia is a compiled language. This means that the first time a function is run, it is compiled into machine code, which takes time. Therefore, the first iteration of `mice()` will be (much) slower in Julia than in R, for example. However, subsequent iterations will be much faster, as all of the required functions are already compiled.

### Why is the first iteration so much slower than the rest?

See [above](#why-is-micejl-so-slow-for-small-jobs).

## Issues

This is a very early work in progress: there will be bugs. If you find any (or the performance is not what you expect), please report them in the [Issues](https://github.com/tom-metherell/Mice.jl/issues) tab above.

## References
<a id="1">[1]</a>
van Buuren S, Groothuis-Oudshoorn K. 2011. mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software* **45**(3):1-67.

<a id="2">[2]</a>
van Buuren S (ed.). 2018. Flexible imputation of missing data. 2nd ed. New York: Chapman & Hall/CRC. 444 p. ISBN: 9780429492259.

<a id="3">[3]</a>
Dickson E, Grambsch P, Fleming T, Fisher L, Langworthy A. 1989. Prognosis in primary biliary cirrhosis: Model for decision making. *Hepatology* **10**(1):1-7.
