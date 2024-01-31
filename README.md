# Mice.jl
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://tom-metherell.github.io/Mice.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://tom-metherell.github.io/Mice.jl/dev)

## What is Mice.jl?

`Mice.jl` is a Julia package for multiple imputation using chained equations. It is heavily based on the R package [mice](https://cran.r-project.org/web/packages/mice/index.html) by Stef van Buuren and Karin Groothuis-Oudshoorn [[1]](#1).

`Mice.jl`'s syntax is very similar to that of `mice` and it is intended to be used in conjunction with [Tables.jl](https://github.com/JuliaData/Tables.jl)-compatible tables, including [DataFrames](https://github.com/JuliaData/DataFrames.jl).

## Installation

To install the latest stable version:

```julia
] add Mice
```

or

```julia
using Pkg; Pkg.add("Mice")
```

## Quick-start guide

### Imputation (`mice()`)
Use the `mice()` function to perform multiple imputation on a data table. The output is a multiply imputed dataset (`Mids`) object.

#### Usage
```julia
mice(data, ...)
```
where:

`data` is the data table.

`m` is the number of imputations.

`imputeWhere` is an `AxisVector` of vectors specifying which values of each variable are to be imputed. If not specified, all missing values are imputed. You can create a default `imputeWhere` vector (which you can then edit) using the function `findMissings(data)`.

`visitSequence` is a vector of column names specifying the order in which the columns should be imputed. If not specified, the order is determined automatically. You can skip the imputation of a column by removing it from the `visitSequence`.

`methods` is an `AxisVector` of imputation methods. If not specified, the default methods are used. Currently, `Mice.jl` supports only a few methods (`"pmm"`, `"norm"`, `"mean"` and `"sample"`). You can make a default methods vector (which you can then edit) using the function `makeMethods(data)`.

`predictorMatrix` is an `AxisMatrix` of predictors for each column (with the predictors in the columns of the matrix). If not specified, all other columns are used for each column. To prevent one column from predicting another, ensure that value of the corresponding cell is set to `0`. You can make a default predictor matrix (which you can then edit) using the function `makePredictorMatrix(data)`.

`iter` is the number of iterations to perform.

`progressReports` is a boolean indicating whether to print progress reports.

`gcSchedule` determines how often the garbage collector is invoked (additionally to when it would be anyway). This can improve performance for large datasets. The value is the fraction of your machine's RAM that remains free at which the garbage collector will be invoked. The default is `0.3`, but you may wish to increase this for very large jobs.

#### Example
```julia
using CSV, DataFrames, Mice, Random

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df)
```

After imputation, you can use `plot(mids, variableNumber)` or `plot(mids, variableName)` to inspect the mean and variance trace plots.

### Repeated analysis (`with()`)
Use the `with()` function to perform data analysis functions on a `Mids` object. The output is a multiply imputed repeated analyses (`Mira`) object.

#### Usage
```julia
with(mids, function)
```
where `mids` is a `Mids` object.

`function` is a data analysis function. It should take the form `data -> function(arguments, data, moreArguments...)`, where the placeholder `data` (which can take other names if necessary) signifies the position of the data argument in the function.

#### Example
```julia
using CSV, DataFrames, GLM, Mice, Random

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df);

analyses = with(imputedData, data -> glm(@formula(y ~ x1 + x2), data, Poisson(), LogLink()))
```

### Pooling coefficients (`pool()`)
Use the `pool()` function to pool the parameter estimates from a `Mira` object. The output is a multiply imputed pooled outcomes (`Mipo`) object.

#### Usage
```julia
pool(mira)
```
where `mira` is a `Mira` object.

#### Example
```julia
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

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric).
In `Mice.jl`, `gcSchedule` was set to `0.3`.

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

### Why is `Mice.jl` so slow for small jobs?

Julia is a just-in-time compiled language. This means that the first time a function is run, it is compiled into machine code, which takes time. Therefore, the first iteration of `mice()` will be (much) slower in Julia than in R, for example. However, subsequent iterations will be much faster, as all of the required functions are already compiled.

### Why is the first iteration so much slower than the rest?

See [above](#why-is-micejl-so-slow-for-small-jobs).

## Issues

This is still a work in progress: there will be bugs. If you find any (or the performance is not what you expect), please report them in the [Issues](https://github.com/tom-metherell/Mice.jl/issues) tab above.

## References
<a id="1">[1]</a>
van Buuren S, Groothuis-Oudshoorn K. 2011. mice: Multivariate Imputation by Chained Equations in R. *Journal of Statistical Software* **45**(3):1-67.

<a id="2">[2]</a>
van Buuren S (ed.). 2018. Flexible imputation of missing data. 2nd ed. New York: Chapman & Hall/CRC. 444 p. ISBN: 9780429492259.

<a id="3">[3]</a>
Dickson E, Grambsch P, Fleming T, Fisher L, Langworthy A. 1989. Prognosis in primary biliary cirrhosis: Model for decision making. *Hepatology* **10**(1):1-7.

<br>
<div align="right"> 
    <span style=""> Funded by Wellcome &nbsp;&nbsp;&nbsp; </span> 
    <img src="docs/wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> 
</div>