# Integration with R

If you don't want to use Julia for your data analysis, but you'd like to try out `Mice.jl`, you can use it from R! The R package [`JuliaCall`](https://non-contradiction.github.io/JuliaCall/index.html) [li_juliacall_2019](@cite) and the Julia package [`RCall.jl`](https://juliainterop.github.io/RCall.jl/stable/) allow you to call Julia functions from R and vice-versa.

!!! note
    Using `Mice.jl` from R requires Julia v1.9 or higher to be installed on your computer. If you don't have it, you can install it from [here](https://julialang.org/downloads/). You also need to install the Julia package `RCall.jl` in addition to `Mice.jl` by entering in Julia:
    ```
    ] add https://github.com/tom-metherell/Mice.jl.git#v0.2.0, RCall
    ```

Using the combined power of these two packages, you can use `Mice.jl` without leaving R. Here's how:

```r
> library(JuliaCall)

> julia_setup()
# Julia version 1.9.4 at location /home/~USERNAME~/julia-1.9.4/bin will be used.
# Loading setup script for JuliaCall...
# Finish loading setup script for JuliaCall.

> julia_library("Mice")

> julia_library("RCall")

> julia_library("GLM")

> julia_library("Random")

> data <- read.csv("test/data/cirrhosis.csv")

> data$Stage <- as.factor(data$Stage)

> julia_assign("data", data)

> julia_command("myPredictorMatrix = makePredictorMatrix(data);")

> julia_command('myPredictorMatrix[:, ["ID", "N_Days"]] .= 0;')

> julia_command("Random.seed!(1234);") # Set random seed for reproducibility

> julia_command("imputedData = mice(data, predictorMatrix = myPredictorMatrix);")

> julia_command("analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));")

> julia_command("results = pool(analyses);")
```

You can also start in Julia, perform your data wrangling in R, use `Mice.jl` and then send your Mids object back to R and continue analysing it there. For example:

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