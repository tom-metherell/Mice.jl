module Mice

    # Dependencies
    using Distributions, LinearAlgebra, Random, StatsBase, Statistics

    # All functions and structs except the main 'mice' function are defined in this file
    include("micehelperfunctions.jl")

"""
    mice(data[, m, visitSequence, methods, predictorMatrix, iter])

Imputes missing values in a dataset using the MICE algorithm. Heavily based on the R package `mice` (Buuren & Groothuis-Oudshoorn, 2011).

Currently, only predictive mean matching (pmm) is supported. All variables in the dataset will be imputed unless the `method` for that variable is specified as an empty string ("").

# Arguments
- `data`: The dataset to impute missing values in. Must be supplied as a `DataFrame`.
- `m`: The number of imputations to generate. Default is 5.
- `visitSequence`: Vector of variables names from the dataset in the order in which they should be imputed. Default is the order in which they appear in the dataset.
- `methods`: Vector of imputation methods to use for each variable (in the order specified by `visitSequence`). Default is `pmm` for all variables. "" will skip imputation for that variable.
- `predictorMatrix`: Boolean matrix indicating which variables should be used as predictors for each variable (in the order specified by `visitSequence`). Default is to use all other variables as predictors for each variable.
- `iter`: Number of iterations. Default is 10.
"""
    function mice(
        data::DataFrame, 
        m::Int = 5,
        visitSequence::AbstractVector = nothing,
        methods::AbstractVector = nothing,
        predictorMatrix::AbstractMatrix = nothing,
        iter::Int = 10
        )

        if visitSequence === nothing
            visitSequence = names(data)
        end

        if methods === nothing
            methods = makeMethods()
        end

        if predictorMatrix === nothing
            predictorMatrix = makePredictorMatrix()
        end

        imputations, meanTraces, varTraces = sampler()

        midsObj = Mids(
            data,
            imputations,
            meanTraces,
            varTraces
        )

        return midsObj
    end

    export Mids, mice
end