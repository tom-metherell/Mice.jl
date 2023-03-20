function makeMethod(data::DataFrame)
    method = similar(data, 1)

    for i in 1:ncol(data)
        if isordered(data[:, i])
            method[1, i] = "polr"
        elseif data[:, i] isa CategoricalArray && length(levels(data[:, i])) > 2
            method[1, i] = "polyreg"
        elseif data[:, i] isa CategoricalArray && length(levels(data[:, i])) == 2
            method[1, i] = "logreg"
        else
            method[1, i] = "pmm"
        end
    end

    return method
end

function makePredictorMatrix(data::DataFrame)
    predictorMatrix = Matrix{Int64}(1, ncol(data), ncol(data))
    setindex!.(Ref(predictorMatrix), 0, 1:ncol(data), 1:ncol(data))
end
