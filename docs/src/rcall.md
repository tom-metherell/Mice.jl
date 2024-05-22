# Integration with R

If you don't want to use Julia for your data analysis, but you'd like to try out `Mice.jl`, good news - you can call R functions from Julia! The Julia package [`RCall.jl`](https://juliainterop.github.io/RCall.jl/stable/) allows you to do this.

You can start in Julia, perform your data wrangling in R, use `Mice.jl` and then send your Mids object back to R and continue analysing it there. For example:

```r
julia> using Mice, Random, RCall

# You can switch from Julia to R by entering $

R> data <- read.csv("test/data/cirrhosis.csv")

R> data$Stage <- as.factor(data$Stage)

# Return to Julia by pressing backspace

julia> @rget data

julia> predictorMatrix = makePredictorMatrix(data);
    
julia> predictorMatrix[:, ["ID", "N_Days"]] .= false;

julia> Random.seed!(1234); # Set random seed for reproducibility

julia> imputedData = mice(data, predictorMatrix = predictorMatrix);

julia> @rput imputedData

R> library(mice)

R> analyses <- with(imputedData, lm(N_Days ~ Drug + Age + Stage + Bilirubin))

R> results <- summary(pool(analyses))
```

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```