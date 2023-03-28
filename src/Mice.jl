module Mice

using DataFrames

include("helperfunctions.jl")

function mice(
    data::DataFrame, 
    m = 5,
    method = nothing,
    predictorMatrix = nothing,
    iter = 10)

    if method === nothing
        method = makeMethod(data)
    end

    if predictorMatrix === nothing
        predictorMatrix = makePredictorMatrix(data)
    end
end