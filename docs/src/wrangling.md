# Data wrangling in Julia
The procedures used for data wrangling in Julia are often very similar to those used in R. However, there are some fundamental differences between the languages that it is important to be aware of.

## Better sources of information than this
`Mice.jl` is intended to be used with `DataFrames.jl`. The [documentation](https://dataframes.juliadata.org/stable/man/missing/) for that package contains a comprehensive overview of using DataFrames.

## Data types
Julia has a number of different data types. These include integers, floats, strings (known in R as "characters"), booleans and missing values.

In a `DataFrame`, each column will have at least one type assigned to it. For example:

```julia
using DataFrames

myData = DataFrame(
    :col1 => [1.0, 2.0, 3.0, 4.0],
    :col2 => [1, 2, 3, 4],
    :col3 => ["1", "2", "3", "4"],
    :col4 => [1.0, 2.0, missing, 4.0],
    :col5 => [1, 2.0, "3", missing]
);

myData.col1
# 4-element Vector{Float64}:
#  1.0
#  2.0
#  3.0
#  4.0

myData.col2
# 4-element Vector{Int64}:
#  1
#  2
#  3
#  4

myData.col3
# 4-element Vector{String}:
#  "1"
#  "2"
#  "3"
#  "4"

myData.col4
# 4-element Vector{Union{Missing, Float64}}:
#  1.0
#  2.0
#   missing
#  4.0

myData.col5
# 4-element Vector{Any}:
#  1
#  2.0
#   "3"
#   missing
```

The first three columns have a single type assigned to them. The fourth column has a `Union` type, which means that it can be either a `Float64` or `missing`. The fifth column has an `Any` type, which means that it can be any type.

If you try to replace a value in a column with a value of the wrong type, one of two things will happen:

```julia
myData[3, :col3] = 5
# ERROR: MethodError: Cannot `convert` an object of type Int64 to an object of type String
# ...

myData[3, :col1] = 5;
myData.col1
# 4-element Vector{Float64}:
#  1.0
#  2.0
#  5.0
#  4.0
```

As you can see, the first example throws an error, because you cannot convert an integer to a string. The second example works, because you can convert an integer to a float.

You can amend a column to allow missing values using the `allowmissing` function:

```julia
myData.col1 = allowmissing(myData.col1)
# 4-element Vector{Union{Missing, Float64}}:
#  1.0
#  2.0
#  5.0
#  4.0
```

or you can achieve this in place using the `allowmissing!` function:

```julia
allowmissing!(myData, 3);
myData.col3
# 4-element Vector{Union{Missing, String}}:
#  "1"
#  "2"
#  "3"
#  "4"
```

## Preparing your data for `mice`

The `mice` function in `Mice.jl` requires a `DataFrame` as its first argument. Columns which contain values of the `missing` type will be imputed (unless you specify otherwise: more on this later).

(Continuous) numeric variables don't require any special treatment. However, categorical variables are handled differently and may require some additional preparation.

`Mice.jl` is compatible with [`CategoricalArrays.jl`](https://categoricalarrays.juliadata.org/stable/), which allows you to specify that a column is categorical. For example:

```julia
using CategoricalArrays, DataFrames

myData = DataFrame(
    :col1 => [1, 2, 2, 3],
    :col2 => categorical([1, 2, 2, 3])
);

myData.col1
# 4-element Vector{Int64}:
#  1
#  2
#  2
#  3

myData.col2
# 4-element CategoricalArray{Int64,1,UInt32}:
#  1
#  2
#  2
#  3
```

As you can see, the `col2` column is now a `CategoricalArray`. This is a special type of array that allows you to specify that a column is categorical. This is important, because `mice` will treat categorical variables differently to continuous variables.

If you have a column that contains strings (with or without missings, but with no other types), `Mice.jl` will handle it as a categorical variable automatically.