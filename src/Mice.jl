module Mice

    # Dependencies
    using Distributions, LinearAlgebra, Random, StatsBase, Statistics

    # All functions and structs except the main 'mice' function are defined in this file
    include("micehelperfunctions.jl")

"""
    mice(
        data::DataFrame, 
        m::Int = 5, 
        visitSequence::Vector{String} = names(data), 
        methods::Vector{String} = nothing, 
        predictorMatrix::Matrix{Bool} = nothing, 
        iter::Int = 10
    )

Imputes missing values in a dataset using the MICE algorithm. 
Heavily based on the R package `mice` (Buuren & Groothuis-Oudshoorn, 2011).

The number of imputations created is specified by `m`.

The variables will be imputed in the order specified by `visitSequence`. 
The default is the order in which they appear in the dataset; 
the order can be customised using a vector of variable names in the desired order.

The imputation method for each variable is specified by `methods`. 
The default is to use predictive mean matching (`pmm`) for all variables. 
Currently only `pmm` is supported. 
Any variable not to be imputed can be marked as such using an empty string ("").

The predictor matrix is specified by `predictorMatrix`. 
The default is to use all other variables as predictors for each variable. 
Any variable not predicting another variable can be marked as such in the matrix using a 0.

The number of iterations is specified by `iter`.
"""
    function mice(
        data::DataFrame, 
        m::Int = 5,
        visitSequence::Vector{String} = names(data),
        methods::Vector{String} = nothing,
        predictorMatrix::Matrix{Bool} = nothing,
        iter::Int = 10
        )

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