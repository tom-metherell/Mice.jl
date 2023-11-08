# Pooling coefficients (`pool`)

Once you have a `Mira` object containing the results of repeated analyses, you can use the `pool` function to pool the results. The `pool` function returns the pooled results wrapped in a `Mipo` object.

```@docs
Mipo
pool
```

The `pool` function should work on any `Mira` of model outputs that accept the StatsAPI functions `coef`, `stderror` and `nobs`. Otherwise, you will get an error and you will need to pool the results manually in accordance with Rubin's rules [rubin_multiple_1987](@cite).

For example:

```julia
julia> using CSV, DataFrames, GLM, Mice, Random

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

julia> analysesLMs = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));
# returns Mira of linear model outputs from each imputed dataset

julia> resultsLMs = pool(analysesLMs);
# returns Mipo of pooled linear model results
```