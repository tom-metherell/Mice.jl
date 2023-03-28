function makeMethod(data = data)
    method = fill("pmm", ncol(data))
end

function makePredictorMatrix(data = data)
    predictorMatrix = Matrix{Int8}(fill(1, ncol(data), ncol(data)))
    for i in 1:ncol(data)
        predictorMatrix[i, i] = 0
    end
end

function initialiseImputations(
    data = data,
    m = m,
    visitSequence = visitSequence,
    method = method
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    presentData = ismissing.(data) .== 0

    for i in 1:length(visitSequence)
        var = visitSequence[i]
        if method[i] != ""
            relevantData = data[:, [var]][:, 1]
            presentLocations = BitVector(presentData[:, [var]][:, 1])
            imputations[i] = Matrix{nonmissingtype(eltype(relevantData))}(undef, sum(.!presentLocations), m)
            if sum(presentLocations) > 0
                for j in 1:m
                    imputations[i][:, j] = sample(relevantData[presentLocations], sum(.!presentLocations))
                end
            else
                if isa(relevantData, CategoricalArray)
                    for j in 1:m
                        imputations[i][:, j] = CategoricalArray{nonmissingtype(eltype(relevantData))}(sample(levels(relevantData), sum(.!presentLocations)))
                    end
                else
                    for j in 1:m
                        imputations[i][:, j] .= randn(sum(.!presentLocations))
                    end
                end
            end
        end
    end

    return imputations
end