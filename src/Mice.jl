module Mice

using StatsBase

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5,
    method = nothing,
    method = nothing,
    predictorMatrix = nothing,
    visitSequence = nothing,
    iter = 10)
    visitSequence = nothing,
    iter = 10)

    if method === nothing
        method = makeMethod()
        method = makeMethod()
    end

    if predictorMatrix === nothing
        predictorMatrix = makePredictorMatrix()
    end

    if visitSequence === nothing
        visitSequence = names(data)
        predictorMatrix = makePredictorMatrix()
    end

    if visitSequence === nothing
        visitSequence = names(data)
    end

    imputations = initialiseImputations()

    imputations = initialiseImputations()
end