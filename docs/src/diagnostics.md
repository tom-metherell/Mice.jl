# Diagnostics

After performing multiple imputation, you should inspect the trace plots of the imputed variables to verify convergence. `Mice.jl` includes a plotting function to do this.

```@docs
plot
```

You do need to load the package `Plots.jl` to see the plots:

```julia
using Plots

# Not run
plot(myMids, 7)
```

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```