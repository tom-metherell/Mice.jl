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

function initialiseChain(
    visitSequence = visitSequence,
    iter = iter,
    m = m
    )

    chain = [[Vector{Union{Missing, Float64}}(undef, length(visitSequence)) for _ = 1:iter] for _ = 1:m]

    return chain
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

    meanChain = initialiseChain()
    varChain = deepcopy(meanChain)

    iterCounter = 1

    while iterCounter <= iter
        for i in 1:m
            for j in eachindex(visitSequence)
                var = visitSequence[j]

                if methods[j] != ""
                    imputations[j][:, i] = univariateSampler(method = methods[j])

                end
            end
        end
    end
end

function univariateSampler(
    data = data,
    presentData = presentData,
    method::String,
    var = var
    )    

    presentLocations = BitVector(presentData[:, [var]][:, 1])