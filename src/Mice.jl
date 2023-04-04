module Mice

using StatsBase

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5,
    visitSequence = nothing,
    methods = nothing,
    predictorMatrix = nothing,
    iter = 10
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