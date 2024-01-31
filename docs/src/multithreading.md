# Multithreading

To start with, you need to make sure that your Julia session is started with multiple threads. See [here](https://docs.julialang.org/en/v1/manual/multi-threading/) for information on how to do this.

As of v0.3.0, you need to run the entire `mice()` function in parallel to get the full benefit of multithreading. It's advisable to set `progressReports = false`. For example, you could do something like this:

```julia
using CategoricalArrays, CSV, DataFrames, Mice, Random

myData = CSV.read("test/data/cirrhosis.csv", DataFrame, missingstring = "NA");

myData.Stage = categorical(myData.Stage); # Making the Stage variable categorical

myPredictorMatrix = makePredictorMatrix(myData);
myPredictorMatrix[:, ["ID", "N_Days"]] .= 0;

Random.seed!(1234); # Set random seed for reproducibility

imputedData = Vector{Mids}(undef, 10); # Initialise vector of Mids outputs

Threads.@threads for i in 1:10 # Number of parallel runs
    # Produces 5 x 10 = 50 imputed datasets in 10 separate Mids objects
    imputedData[i] = mice(myData, m = 5, predictorMatrix = myPredictorMatrix, progressReports = false)
end

imputedData = bindImputations(imputedData); # Binds the separate Mids objects into a single output
```

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```