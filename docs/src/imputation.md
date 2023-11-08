# Imputation (`mice`)
The main function of the package is `mice`, which takes a `DataFrame` as its input. It returns a multiply imputed dataset (`Mids`) object with the imputed values.

```@docs
Mids
mice
```

```@docs
makeMonotoneSequence
```

```@docs
makePredictorMatrix
```

```@docs
makeMethods
```

```@docs
initialiseTraces
```

```@docs
initialiseImputations
```

## Diagnostics
After performing multiple imputation, you should inspect the trace plots of the imputed variables to verify convergence. `Mice.jl` includes the a plotting function to do this.

```@docs
plot
```