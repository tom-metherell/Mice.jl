# Analysis (`with`)
Once you have a `Mids` object containing imputed data, you can use it to perform repeated analyses.

## Inspecting imputed data

If you just want to inspect the outcome of the imputation process, you can use the `complete`/`listComplete` function to fill in the missing values in the original data frame.

```@docs
complete
listComplete
```

## Data analysis

To perform a data analysis procedure on each imputed dataset in turn, use the `with` function. The `with` function returns the results of the analyses wrapped in a `Mira` object.

```@docs
Mira
with
```

The `with` function requires the use of a closure, which then permits the function to run the specified analysis procedure on each imputed dataset in turn. For example:

```julia
using CSV, DataFrames, GLM, Mice, Random, Statistics

myData = CSV.read("test/data/cirrhosis.csv", DataFrame, missingstring = "NA");

myData.Stage = categorical(myData.Stage); # Making the Stage variable categorical

myPredictorMatrix = makePredictorMatrix(myData);

myPredictorMatrix[:, ["ID", "N_Days"]] .= false;

Random.seed!(1234); # Set random seed for reproducibility

imputedData = mice(myData, predictorMatrix = myPredictorMatrix, methods = myMethods);

analysesMeans = with(imputedData, data -> mean(data.Cholesterol));
# returns Mira of the mean of Bilirubin in each imputed dataset

analysesLMs = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));
# returns Mira of linear model outputs from each imputed dataset
```