function makeMethod(data::DataFrame)
    method = fill("pmm", ncol(data))
end

function makePredictorMatrix(data::DataFrame)
    predictorMatrix = DataFrame(Matrix{Int8}(undef, ncol(data), ncol(data)), names(data))
    predictorMatrix[:, :] .= 1
    setindex!.(Ref(predictorMatrix), 0, 1:ncol(data), 1:ncol(data))
end