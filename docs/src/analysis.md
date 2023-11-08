# Analysis (`with`)
Once you have a `Mids` object containing imputed data, you can use it to perform repeated analyses.

## Inspecting imputed data

If you just want to inspect the outcome of the imputation process, you can use the `complete` function to fill in the missing values in the original data frame.

```@docs
complete
```

## Data analysis

To perform a data analysis procedure on each imputed dataset in turn, use the `with` function. The `with` function returns the results of the analyses wrapped in a `Mira` object.

```@docs
Mira
with
```

The `with` function requires the use of a closure, which then permits the function to run the specified analysis procedure on each imputed dataset in turn. For example:

```julia
julia> using CSV, DataFrames, GLM, Mice, Random, Statistics

julia> myData = CSV.read("test/data/cirrhosis.csv", DataFrame);

# Defining missing values
julia> colsWithMissings = ["Drug", "Ascites", "Hepatomegaly", "Spiders", "Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin", "Stage"];
julia> myData[!, colsWithMissings] = allowmissing(myData[!, colsWithMissings]);
julia> for i in colsWithMissings
    replace!(myData[!, i], "NA" => missing)
end
julia> for i in ["Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]
    myData[!, i] = passmissing(x -> parse(Float64, x)).(myData[!, i])
end

julia> myMethods = makeMethods(myData);
julia> myMethods[["ID", "N_Days"]] .= "";

julia> myPredictorMatrix = makePredictorMatrix(myData);
julia> myPredictorMatrix[:, ["ID", "N_Days"]] .= false;

julia> Random.seed!(1234); # Set random seed for reproducibility

julia> imputedData = mice(myData, predictorMatrix = myPredictorMatrix, methods = myMethods);

julia> analysesMeans = with(imputedData, data -> mean(data.Cholesterol));
# returns Mira of the mean of Bilirubin in each imputed dataset

julia> analysesLMs = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));
# returns Mira of linear model outputs from each imputed dataset
```