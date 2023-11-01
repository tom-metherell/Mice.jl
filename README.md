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

### Repeated analysis (`with()`)

### Pooling coefficients (`pool()`)