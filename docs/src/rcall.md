# Integration with R

If you don't want to use Julia for your data analysis, but you'd like to try out `Mice.jl`, you can use it from R! The R package [`JuliaCall`](https://non-contradiction.github.io/JuliaCall/index.html) [li_juliacall_2019](@cite) and the Julia package [`RCall.jl`](https://juliainterop.github.io/RCall.jl/stable/) allow you to call Julia functions from R and vice-versa.

`Mice.jl` is has an `RCall.jl` extension, allowing you to send your Mids objects to R and continue analysing them there. For example:

```r
julia> using CategoricalArrays, CSV, DataFrames, Mice, RCall

julia> myData = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA");

julia> myData.Stage = categorical(myData.Stage); # Making the Stage variable categorical

julia> predictorMatrix = makePredictorMatrix(data)
    
julia> predictorMatrix[:, ["ID", "N_Days"]] .= false

julia> imputedData = mice(myData, predictorMatrix = predictorMatrix, methods = myMethods);

julia> @rput imputedData

R> library(mice)

R> analyses <- with(imputedData, lm(N_Days ~ Drug + Age + Stage + Bilirubin))

R> results <- summary(pool(analyses))
```



