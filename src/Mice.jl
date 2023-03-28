module Mice

using StatsBase

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5,
    method = nothing,
    predictorMatrix = nothing,
    visitSequence = nothing,
    iter = 10)

    if method === nothing
        method = makeMethod()
    end

    if predictorMatrix === nothing
        predictorMatrix = makePredictorMatrix()
    end

    if visitSequence === nothing
        visitSequence = names(data)
    end

    imputations = initialiseImputations()
end