# Mice.jl

## What is Mice.jl?

`Mice.jl` is a Julia package for multiple imputation using chained equations. It is heavily based on the R package [mice](https://cran.r-project.org/web/packages/mice/index.html) by Stef van Buuren and Karin Groothuis-Oudshoorn.

`Mice.jl`'s syntax is very similar to that of `mice` and it is intended to be used in conjunction with [DataFrames.jl](https://github.com/JuliaData/DataFrames.jl).

## Installation

```
] add Mice
```

or

```
using Pkg; Pkg.add("Mice")
```

## Quick-start guide

### Imputation (`mice()`)
Use the `mice()` function to perform multiple imputation on a DataFrame. The output is a multiply imputed dataset (`Mids`) object.

#### Usage
```
mice(df, m = 5, visitSequence = nothing, methods = nothing, predictorMatrix = nothing, iter = 10, progressReports = true, ...)
```
where:

`df` is the DataFrame.

`m` is the number of imputations.

`visitSequence` is a vector of column names specifying the order in which the columns should be imputed. If `nothing`, the order is determined automatically.

`methods` is a vector of imputation methods. If `nothing`, the default methods are used. Currently, `Mice.jl` supports predictive mean matching (`pmm`) only. To skip a column, define the method for that column as an empty string (`""`).

`predictorMatrix` is a matrix of predictors for each column (with the predictors in the columns of the matrix). If `nothing`, all other columns are used for each column. To prevent one column from predicting another, ensure that value of the corresponding cell is set to `0`.

`iter` is the number of iterations to perform.

`progressReports` is a boolean indicating whether to print progress reports.

#### Example
```
using CSV, Mice

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df)
```

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
using CSV, GLM, Mice

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
using CSV, GLM, Mice

Random.seed!(1234)

df = CSV.read("my_data.csv", DataFrame);

imputedData = mice(df);

analyses = with(imputedData, data -> glm(@formula(y ~ x1 + x2), data, Poisson(), LogLink()));

results = pool(analyses)
```