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

myData = CSV.read("test/data/cirrhosis.csv", DataFrame);

# Defining missing values
colsWithMissings = ["Drug", "Ascites", "Hepatomegaly", "Spiders", "Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin", "Stage"];
myData[!, colsWithMissings] = allowmissing(myData[!, colsWithMissings]);
for i in colsWithMissings
    replace!(myData[!, i], "NA" => missing)
end
for i in ["Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]
    myData[!, i] = passmissing(x -> parse(Float64, x)).(myData[!, i])
end

myMethods = makeMethods(myData);
myMethods[["ID", "N_Days"]] .= "";

myPredictorMatrix = makePredictorMatrix(myData);
myPredictorMatrix[:, ["ID", "N_Days"]] .= false;

Random.seed!(1234); # Set random seed for reproducibility

imputedData = mice(myData, predictorMatrix = myPredictorMatrix, methods = myMethods);

analysesLMs = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data));
# returns Mira of linear model outputs from each imputed dataset

resultsLMs = pool(analysesLMs);
# returns Mipo of pooled linear model results
```