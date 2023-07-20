module Mice

using Distributions, LinearAlgebra, Random, StatsBase, Statistics

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5::Int,
    visitSequence = nothing,
    methods = nothing,
    predictorMatrix = nothing, # methods and predictorMatrix must be arranged in visitSequence order!
    iter = 10::Int
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
end