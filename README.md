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

`methods` is an `AxisVector` of imputation methods. If not specified, the default methods are used. Currently, `Mice.jl` supports only a few methods (`"pmm"`, `"rf"`, `"norm"`, `"logreg"`, `"mean"` and `"sample"`). You can make a default methods vector (which you can then edit) using the function `makeMethods(data)`.

`predictorMatrix` is an `AxisMatrix` of predictors for each column (with the predictors in the columns of the matrix). If not specified, all other columns are used for each column. To prevent one column from predicting another, ensure that value of the corresponding cell is set to `0`. You can make a default predictor matrix (which you can then edit) using the function `makePredictorMatrix(data)`.

`iter` is the number of iterations to perform.

`progressReports` is a boolean indicating whether to print progress reports.

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

15 iterations were completed to impute 12 variables (of which 4 binary categorical, 1 other categorical and 7 numeric) using a set of 18 predictors (those 12 variables plus 6 complete variables: 1 binary categorical, 2 other categorical and 3 numeric). Both used predictive mean matching for all variables that were to be imputed. For the Julia implementation, a new terminal was used for each test to count the time taken to compile the functions.

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
<table align="right" style="width: 20%; border: none;" cellspacing="0" cellpadding="0" border="0"> 
    <tr>
        <td valign="center"> Funded by Wellcome </td>
        <td valign="center"> <img src="docs/src/wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </td>
    </tr>
</table>
