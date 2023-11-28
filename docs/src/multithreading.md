# Multithreading

`Mice.jl` supports multithreading. To start with, you need to make sure that your Julia session is started with multiple threads. See [here](https://docs.julialang.org/en/v1/manual/multi-threading/) for information on how to do this.

Once this is done, you can run `mice` in multithreaded mode by setting the `threads` argument to `true`. This causes the imputations to be computed in parallel.

If you instead want to run the entire `mice` function in parallel, you can do something like this:

```julia
using CSV, DataFrames, Mice, Random

myData = CSV.read("test/data/cirrhosis.csv", DataFrame, missingstring = "NA");

myData.Stage = categorical(myData.Stage); # Making the Stage variable categorical

myPredictorMatrix = makePredictorMatrix(myData);
myPredictorMatrix[:, ["ID", "N_Days"]] .= false;

Random.seed!(1234); # Set random seed for reproducibility

imputedData = Vector{Mids}(undef, 10); # Initialise vector of Mids outputs

Threads.@threads for i in 1:10 # Number of parallel runs
    # Produces 5 x 10 = 50 imputed datasets in 10 separate Mids objects
    imputedData[i] = mice(myData, m = 5, predictorMatrix = myPredictorMatrix, methods = myMethods, threads = false, progressReports = false)
end

imputedData = bindImputations(imputedData); # Binds the separate Mids objects into a single output
```