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

    meanTraces = initialiseTraces()
    varTraces = deepcopy(meanTraces)

    iterCounter = 1

    while iterCounter <= iter
        for i in eachindex(visitSequence)
            var = visitSequence[i]
            y = data[:, [var]][:, 1]
            X = data[:, predictorMatrix[i, :]]

            if methods[i] == "pmm" && any(ismissing(y))

                for j in 1:m
                    for k in findall(predictorMatrix[i, :])
                        xVar = visitSequence[k]
                        X[ismissing.(X) .== 1, [xVar]] = imputations[k][:, j]
                    end

                    imputations[i][:, j] = pmmImpute()
                    
                    plottingData = deepcopy(data[:, [var]][:, 1])
                    plottingData[ismissing.(plottingData) .== 1] = imputations[i][:, j]
                end

                if plottingData isa CategoricalArray
                    mapping = Dict(levels(plottingData)[i] => i for i in eachindex(levels(plottingData)))
    
                    plottingData = [mapping[v] for v in plottingData]
                end
    
                meanTraces[i][iterCounter][j] = mean(plottingData)
                varTraces[i][iterCounter][j] = var(plottingData)
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

    X = hcat(repeat([1], size(X, 1)), X)

    if y isa CategoricalArray
        mapping = Dict(levels(y)[i] => i for i in eachindex(levels(y)))

        y = [mapping[v] for v in y]
    end

    Xₒ = X[ismissing.(y) .== 0, :]
    Xₘ = X[ismissing.(y) .== 1, :]
    yₒ = y[ismissing.(y) .== 0]

    β̂, β̇, = blrDraw()

    ŷₒ = Xₒ * β̂
    ŷₘ = Xₘ * β̇
    
end

function blrDraw(
    yₒ = yₒ,
    Xₒ = Xₒ, 
    κ = 0.0001::Float
    )

    S = transpose(Xₒ) * Xₒ

    V = inv(S + diag(S) * κ)

    β̂ = V * transpose(Xₒ) * yₒ

    ġ = rand(Chisq(size(Xₒ, 1) - size(Xₒ, 2)))

    σ̇ = sqrt((transpose(yₒ - Xₒ * β̂) * (yₒ - Xₒ * β̂)) / ġ)

    ż = randn(size(Xₒ, 2))

    sqrtV = cholesky(V)

    β̇ = β̂ + σ̇ * ż * sqrtV

    return β̂, β̇
end

using Random, StatsBase

function matchindex(yₒ::Vector{Real}, yₘ::Vector{Real}, donors = 5::Int)
    # Shuffle records to remove effects of ties
    nₒ = length(yₒ)
    ishuf = randperm(nₒ)
    yshuf = yₒ[ishuf]

    # Obtain sorting order on shuffled data
    isort = sortperm(yshuf)

    # Calculate index on input data and sort
    id = ishuf[isort]
    ysort = yₒ[id]

    # Pre-sample n0 values between 1 and k
    nₘ = length(yₘ)
    donors = min(donors, nₒ)
    donors = max(donors, 1)
    presample = sample(1:donors, nₘ, replace = true)

    indices = similar(yₘ, Int)

    # Loop over the target units
    for i in eachindex(yₘ)
        value = yₘ[i]
        donor = presample[i]
        count = 0

        # Find the two adjacent neighbours
        r = searchsortedfirst(ysort, value)
        l = r - 1

        ### CHECKED UP TO HERE

        # Find the h_i'th nearest neighbour
        # Store the index of that neighbour
        while count < hi && l >= 1 && r <= n1
            if val - ysort[l] < ysort[r] - val
                idx[i] = id[l]
                l -= 1
            else
                idx[i] = id[r]
                r += 1
            end
            count += 1
        end

        # If right side is exhausted, take left elements
        while count < hi && l >= 1
            idx[i] = id[l]
            l -= 1
            count += 1
        end

        # If left side is exhausted, take right elements
        while count < hi && r <= n1
            idx[i] = id[r]
            r += 1
            count += 1
        end
    end

    return idx
end