module Mice

using StatsBase, Statistics

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5::Int,
    visitSequence = nothing,
    methods = nothing,
    predictorMatrix = nothing,
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