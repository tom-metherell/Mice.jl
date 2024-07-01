# Imputation (`mice`)

The main function of the package is `mice`, which takes a `Tables.jl`-compatible table as its input. It returns a multiply imputed dataset (`Mids`) object with the imputed values.

```@docs
Mids
mice
```