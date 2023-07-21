module Mice

    using Distributions, LinearAlgebra, Random, StatsBase, Statistics

    include("micehelperfunctions.jl")

    function mice(
        data::DataFrame, 
        m::Int = 5,
        visitSequence::AbstractVector = nothing,
        methods::AbstractVector = nothing,
        predictorMatrix::AbstractMatrix = nothing, # methods and predictorMatrix must be arranged in visitSequence order!
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
end