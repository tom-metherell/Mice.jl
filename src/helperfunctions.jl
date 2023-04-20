function makeMethods(data = data)
    methods = fill("pmm", ncol(data))

    return methods
end

function makePredictorMatrix(data = data)
    predictorMatrix = Matrix{Bool}(fill(1, ncol(data), ncol(data)))
    for i in 1:ncol(data)
        predictorMatrix[i, i] = 0
    end

    return predictorMatrix
end

function initialiseImputations(
    data = data,
    m = m,
    visitSequence = visitSequence,
    methods = methods
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    presentData = ismissing.(data) .== 0

    for i in eachindex(visitSequence)
        var = visitSequence[i]
        if methods[i] != ""
            relevantData = data[:, [var]][:, 1]
            presentLocations = BitVector(presentData[:, [var]][:, 1])
            presentDataCount = sum(.!presentLocations)
            imputations[i] = Matrix{nonmissingtype(eltype(relevantData))}(undef, missingDataCount, m)
            if sum(presentLocations) > 0
                for j in 1:m
                    imputations[i][:, j] = sample(relevantData[presentLocations], missingDataCount)
                end
            else
                if isa(relevantData, CategoricalArray)
                    for j in 1:m
                        imputations[i][:, j] = CategoricalArray{nonmissingtype(eltype(relevantData))}(sample(levels(relevantData), presentDataCount))
                    end
                else
                    for j in 1:m
                        imputations[i][:, j] .= randn(presentDataCount)
                    end
                end
            end
        end
    end

    return imputations
end

function initialiseTraces(
    visitSequence = visitSequence,
    iter = iter,
    m = m
    )

    traces = [[Vector{Union{Missing, Float64}}(undef, m) for _ = 1:iter] for _ = eachindex(visitSequence)]

    return traces
end

function sampler(
    data = data,
    m = m,
    imputations = initialiseImputations(),
    methods = methods,
    visitSequence = visitSequence,
    predictorMatrix = predictorMatrix,
    iter = iter
    )

    presentData = ismissing.(data) .== 0

    meanTraces = initialiseTraces()
    varTraces = deepcopy(meanTraces)

    iterCounter = 1

    while iterCounter <= iter
        for i in eachindex(visitSequence)
            var = visitSequence[i]
            y = data[:, [var]][:, 1]
            X = data[:, predictorMatrix[i, :]][:, 1]

            if methods[i] == "pmm" && any(ismissing(y))

                for j in 1:m
                    imputations[i][:, j] = pmmImpute()

                    data[presentData .== 0][:, [var]] = imputations[i][:, j]
                end

                if y isa CategoricalArray
                    mapping = Dict(levels(data[:, [var]][:, 1])[i] => i for i in eachindex(levels(data[:, [var]][:, 1])))
    
                    data[:, [var]][:, 1] = [mapping[v] for v in data[:, [var]][:, 1]]
                end
    
                meanTraces[i][iterCounter][j] = mean(data[:, [var]][:, 1])
                varTraces[i][iterCounter][j] = var(data[:, [var]][:, 1])
            end
        end

        iterCounter += 1
    end

    return imputations, meanTraces, varTraces
end

function pmmImpute(
    y = y,
    X = X,
    donors = 5::Int
    )

    X = Matrix(X)

    X = vcat(1, X)

    if y isa CategoricalArray
        mapping = Dict(levels(y)[i] => i for i in eachindex(levels(y)))

        y = [mapping[v] for v in y]
    end

    coefs, β = blrDraw()
    
end

function blrDraw(
    y = y,
    X = X
    )

end