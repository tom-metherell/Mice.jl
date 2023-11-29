# Pooling coefficients (`pool`)

Once you have a `Mira` object containing the results of repeated analyses, you can use the `pool` function to pool the results. The `pool` function returns the pooled results wrapped in a `Mipo` object.

```@docs
Mipo
pool
```

The `pool` function should work on any `Mira` of model outputs that accept the StatsAPI functions `coef`, `stderror` and `nobs`. Otherwise, you will get an error and you will need to pool the results manually in accordance with Rubin's rules [rubin_multiple_1987](@cite).

For example:

```julia
using CSV, DataFrames, GLM, Mice, Random

myData = CSV.read("test/data/cirrhosis.csv", DataFrame, missingstring = "NA");

myData.Stage = categorical(myData.Stage); # Making the Stage variable categorical

myPredictorMatrix = makePredictorMatrix(myData);
myPredictorMatrix[:, ["ID", "N_Days"]] .= false;

Random.seed!(1234); # Set random seed for reproducibility

imputedData = mice(myData, predictorMatrix = myPredictorMatrix, methods = myMethods);

analysesLMs = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));
# returns Mira of linear model outputs from each imputed dataset

resultsLMs = pool(analysesLMs);
# returns Mipo of pooled linear model results
```