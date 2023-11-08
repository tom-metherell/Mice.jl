# Imputation (`mice`)
The main function of the package is `mice`, which takes a `DataFrame` as its input. It returns a multiply imputed dataset (`Mids`) object with the imputed values.

```@docs
Mids
mice
```

## Customising the imputation setup
You can customise various aspects of the imputation setup by passing keyword arguments to `mice`. These are described above. You can also use some of the functions below to define objects that you can customise to alter how `mice` handles the imputation.

### Visit sequence
The visit sequence is the order in which the variables are imputed. By default, `mice` sorts the variables in order of missingness (lowest to highest) via the function `makeMonotoneSequence`.

```@docs
makeMonotoneSequence
```

You can instead define your own visit sequence by creating a vector of variable names in your desired order and passing that to `mice`. For example:

```julia
using DataFrames, Mice, Random

myData = DataFrame(
    :col1 => Vector{Union{Missing, Float64}}([1.0, missing, 3.0, missing, 5.0]),
    :col2 => Vector{Union{Missing, Int64}}([1, 2, missing, 4, 5]),
    :col3 => Vector{Union{Missing, String}}([missing, "2", missing, "4", missing])
);

makeMonotoneSequence(myData)
# 3-element Vector{String}:
# "col2"
# "col1"
# "col3"

myVisitSequence1 = names(myData)
# 3-element Vector{String}:
# "col1"
# "col2"
# "col3"

Random.seed!(1234); # Set random seed for reproducibility

# Not run
mice(myData, visitSequence = myVisitSequence1)

myVisitSequence2 = ["col3", "col1", "col2"]
# 3-element Vector{String}:
# "col3"
# "col1"
# "col2"

# Not run
mice(myData, visitSequence = myVisitSequence2)
```

Assuming that the imputations converge normally, changing the visit sequence should not dramatically affect the output. However, it can be useful to change the visit sequence if you want to impute variables in a particular order for a specific reason. The sequence used by default in `Mice.jl` can make convergence faster in cases where the data follow a (near-)"monotone" missing data pattern [gelman_multivariate_2018](@cite).

### Predictor matrix
The predictor matrix defines which variables in the imputation model are used to predict which others. By default, every variable predicts every other variable, but there are a wide range of cases in which this is not desirable. For example, if your dataset includes an ID column, this is clearly useless for imputation and should be ignored.

To create a default predictor matrix that you can edit, you can use the function `makePredictorMatrix`.

```@docs
makePredictorMatrix
```

You can then edit the predictor matrix to remove any predictive relationships that you do not want to include in the imputation model. For example:

```julia
using DataFrames, Mice, Random

myData = DataFrame(
    :id => Vector{Int64}(1:5),
    :col1 => Vector{Union{Missing, Float64}}([1.0, missing, 3.0, missing, 5.0]),
    :col2 => Vector{Union{Missing, Int64}}([1, 2, missing, 4, 5]),
    :col3 => Vector{Union{Missing, String}}([missing, "2", missing, "4", missing])
);

myPredictorMatrix = makePredictorMatrix(myData)
# 4x4 Named Matrix{Bool}
# A \ B |    id   col1   col2   col3
# ------|---------------------------
# id    | false   true   true   true
# col1  |  true  false   true   true
# col2  |  true   true  false   true
# col3  |  true   true   true  false

# To stop the ID column from predicting any other variable
myPredictorMatrix[:, "id"] .= false;
myPredictorMatrix
# 4x4 Named Matrix{Bool}
# A \ B |    id   col1   col2   col3
# ------|---------------------------
# id    | false   true   true   true
# col1  | false  false   true   true
# col2  | false   true  false   true
# col3  | false   true   true  false

# To stop col1 from predicting col3
myPredictorMatrix["col3", "col1"] = false;
myPredictorMatrix
# 4x4 Named Matrix{Bool}
# A \ B |    id   col1   col2   col3
# ------|---------------------------
# id    | false   true   true   true
# col1  | false  false   true   true
# col2  | false   true  false   true
# col3  | false  false   true  false

Random.seed!(1234); # Set random seed for reproducibility

# Not run
mice(myData, predictorMatrix = myPredictorMatrix)
```

### Methods
The imputation methods are the functions that are used to impute each variable. By default, `mice` uses predictive mean matching (`"pmm"`) for all variables (and currently PMM is the only method that `Mice.jl` supports). However, you can use the methods vector to specify any variables that should not be imputed.

To create a default methods vector, use the function `makeMethods`.

```@docs
makeMethods
```

You can then customise the vector as needed. For example:

```julia
using DataFrames, Mice, Random

myData = DataFrame(
    :id => Vector{Int64}(1:5),
    :col1 => Vector{Union{Missing, Float64}}([1.0, missing, 3.0, missing, 5.0]),
    :col2 => Vector{Union{Missing, Int64}}([1, 2, missing, 4, 5]),
    :col3 => Vector{Union{Missing, String}}([missing, "2", missing, "4", missing])
);

myMethods = makeMethods(myData)
# 4-element Named Vector{String}
# A    |
# -----|------
# id   | "pmm"
# col1 | "pmm"
# col2 | "pmm"
# col3 | "pmm"

# To stop the ID column from being imputed
myMethods["id"] = "";
myMethods
# 4-element Named Vector{String}
# A    |
# -----|------
# id   |    ""
# col1 | "pmm"
# col2 | "pmm"
# col3 | "pmm"

Random.seed!(1234); # Set random seed for reproducibility

# Not run
mice(myData, methods = myMethods)
```

## Diagnostics
After performing multiple imputation, you should inspect the trace plots of the imputed variables to verify convergence. `Mice.jl` includes the a plotting function to do this.

```@docs
plot
```