# Customising the imputation setup

You can customise various aspects of the imputation setup by passing keyword arguments to `mice`. These are described above. You can also use some of the functions below to define objects that you can customise to alter how `mice` handles the imputation.

## Locations to impute
You can customise which data points are imputed by manipulating the `imputeWhere` argument. By default, this will specify that all missing data are to be imputed (using the function `findMissings()`).

```@docs
findMissings
```

You can over-impute existing data by setting the locations of non-missing data to `true` in the relevant vector in `imputeWhere`. For example, to over-impute the value of `col1` for the first row, you could do the following:

```julia
using DataFrames, Mice, Random

myData = DataFrame(
    :col1 => Vector{Union{Missing, Float64}}([1.0, missing, 3.0, missing, 5.0]),
    :col2 => Vector{Union{Missing, Int64}}([1, 2, missing, 4, 5]),
    :col3 => Vector{Union{Missing, String}}([missing, "2", missing, "4", missing])
);

myImputeWhere = findMissings(myData)
# 1-dimensional AxisArray{Vector{Bool},1,...} with axes:
#     :row, ["col1", "col2", "col3"]
# And data, a 3-element Vector{Vector{Bool}}:
#  [0, 1, 0, 1, 0]
#  [0, 0, 1, 0, 0]
#  [1, 0, 1, 0, 1]

myImputeWhere["col1"][1] = true;
myImputeWhere
# 1-dimensional AxisArray{Vector{Bool},1,...} with axes:
#     :row, ["col1", "col2", "col3"]
# And data, a 3-element Vector{Vector{Bool}}:
#  [1, 1, 0, 1, 0]
#  [0, 0, 1, 0, 0]
#  [1, 0, 1, 0, 1]

# Not run
mice(myData, imputeWhere = myImputeWhere)
```

## Visit sequence
The visit sequence is the order in which the variables are imputed. By default, `mice` sorts the variables in order of missingness (lowest to highest) via the internal function `makeMonotoneSequence`. You can instead define your own visit sequence by creating a vector of variable names in your desired order and passing that to `mice`. For example:

```julia
using DataFrames, Mice, Random

myData = DataFrame(
    :col1 => Vector{Union{Missing, Float64}}([1.0, missing, 3.0, missing, 5.0]),
    :col2 => Vector{Union{Missing, Int64}}([1, 2, missing, 4, 5]),
    :col3 => Vector{Union{Missing, String}}([missing, "2", missing, "4", missing])
);

Mice.makeMonotoneSequence(findMissings(myData))
# 3-element Vector{String}:
#  "col2"
#  "col1"
#  "col3"

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

Assuming that the imputations converge normally, changing the visit sequence should not dramatically affect the output. However, it can be useful to change the visit sequence if you want to impute variables in a particular order for a specific reason. The sequence used by default in `Mice.jl` can make convergence faster in cases where the data follow a (near-)"monotone" missing data pattern [van_buuren_flexible_2018](@cite).

You can leave variables out of the `visitSequence` to cause `mice()` to not impute them.

## Predictor matrix
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
# 2-dimensional AxisArray{Int64,2,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
#     :col, ["id", "col1", "col2", "col3"]
# And data, a 4x4 Matrix{Int64}:
#  0  1  1  1
#  1  0  1  1
#  1  1  0  1
#  1  1  1  0

# To stop the ID column from predicting any other variable
myPredictorMatrix[:, "id"] .= 0;
myPredictorMatrix
# 2-dimensional AxisArray{Int64,2,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
#     :col, ["id", "col1", "col2", "col3"]
# And data, a 4x4 Matrix{Int64}:
#  0  1  1  1
#  0  0  1  1
#  0  1  0  1
#  0  1  1  0

# To stop col1 from predicting col3
myPredictorMatrix["col3", "col1"] = false;
myPredictorMatrix
# 2-dimensional AxisArray{Int64,2,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
#     :col, ["id", "col1", "col2", "col3"]
# And data, a 4x4 Matrix{Int64}:
#  0  1  1  1
#  0  0  1  1
#  0  1  0  1
#  0  0  1  0

Random.seed!(1234); # Set random seed for reproducibility

# Not run
mice(myData, predictorMatrix = myPredictorMatrix)
```

## Methods
The imputation methods are the functions that are used to impute each variable. By default, `mice` uses predictive mean matching (`"pmm"`) for all variables. Currently `Mice.jl` supports the following methods:

| Method | Description | Variable type |
| ------ | ----------- | ------------- |
| `pmm` | Predictive mean matching | Any |
| `rf` | Random forest | Any (but see [below](#rf-warning)) |
| `sample` | Random sample from observed values | Any |
| `mean` | Mean of observed values | Numeric (float) |
| `norm` | Bayesian linear regression | Numeric (float) |
| `logreg` | Logistic regression | Binary |

```@raw html
<a name="rf-warning">
</a> 
```

!!! warning

    If you use `rf` on a variable with integer values, the imputed values will be rounded to the nearest integer in the output. If you want to prevent this behaviour, you have two options:

    * Convert the variable to a float before imputation, so it is treated as continuous or

    * Convert the variable to a categorical/string array so it is treated as discrete.

The `mean` and `sample` methods should not generally be used.

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
# 1-dimensional AxisArray{String,1,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
# And data, a 4-element Vector{String}:
#  "pmm"
#  "pmm"
#  "pmm"
#  "pmm"

# To stop the ID column from being imputed (but you can also achieve this by leaving "id"
# out of the visit sequence)
myMethods["id"] = "";
myMethods
# 1-dimensional AxisArray{String,1,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
# And data, a 4-element Vector{String}:
#  ""
#  "pmm"
#  "pmm"
#  "pmm"

# To use Bayesian linear regression to impute col1
myMethods["col1"] = "norm";
myMethods
# 1-dimensional AxisArray{String,1,...} with axes:
#     :row, ["id", "col1", "col2", "col3"]
# And data, a 4-element Vector{String}:
#  ""
#  "norm"
#  "pmm"
#  "pmm"

Random.seed!(1234); # Set random seed for reproducibility

# Not run
mice(myData, methods = myMethods)
```

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```