function makeMonotoneSequence(data::DataFrame = data)
    missingness = Vector{Int64}(undef, size(data, 2))

    # Counting missing data in each column
    for i in axes(data, 2)
        missingness[i] = sum(ismissing.(data[:, i]))
    end

    # Sort the missingness vector in descending order
    missingness = sortperm(missingness, rev = true)

    # Sort the data frame names vector by missingness
    visitSequence = names(data)[missingness]

    return visitSequence
end

function makeMethods(data::DataFrame = data)
    methods = fill("pmm", ncol(data))

    return methods
end

function makePredictorMatrix(data::DataFrame = data)
    predictorMatrix = Matrix{Bool}(fill(1, ncol(data), ncol(data)))
    for i in 1:ncol(data)
        predictorMatrix[i, i] = 0
    end

    return predictorMatrix
end

function initialiseImputations(
    data::DataFrame = data,
    m::Int = m,
    visitSequence::AbstractVector = visitSequence,
    methods::AbstractVector = methods
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    presentData = ismissing.(data) .== 0

    for i in eachindex(visitSequence)
        var = visitSequence[i]
        if methods[i] != ""
            relevantData = data[:, [var]][:, 1]
            presentLocations = BitVector(presentData[:, [var]][:, 1])
            presentDataCount = sum(presentLocations)
            missingDataCount = sum(.!presentLocations)
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
    visitSequence::AbstractVector = visitSequence,
    iter::Int = iter,
    m::Int = m
    )

    traces = [[Vector{Union{Missing, Float64}}(undef, m) for _ = 1:iter] for _ = eachindex(visitSequence)]

    return traces
end

function sampler(
    data::DataFrame = data,
    m::Int = m,
    methods::AbstractVector = methods,
    visitSequence::AbstractVector = visitSequence,
    predictorMatrix::AbstractMatrix = predictorMatrix,
    iter::Int = iter
    )

    imputations = initialiseImputations(data, m, visitSequence, methods)

    meanTraces = initialiseTraces(visitSequence, iter, m)
    varTraces = deepcopy(meanTraces)

    iterCounter = 1

    while iterCounter <= iter
        for i in eachindex(visitSequence)
            var = visitSequence[i]
            y = data[:, [var]][:, 1]
            X = data[:, predictorMatrix[i, :]] # This is wrong

            if methods[i] == "pmm" && any(ismissing(y))

                for j in 1:m

                    # Doesn't work
                    for k in findall(predictorMatrix[i, :])
                        xVar = visitSequence[k]
                        replacements = Vector{Union{Missing,nonmissingtype(eltype(X[:, [xVar]][:, 1]))}}(undef, size(X, 1))
                        for i in 1:size(X, 1)
                            counter = 1
                            if ismissing(X[:, [xVar]][i, 1])
                                replacements[i] = imputations[k][counter, j]
                                counter += 1
                            end
                        end
                        X[:, [xVar]][:, 1] = coalesce(X[:, [xVar]][:, 1], replacements)
                    end

                    imputations[i][:, j] = pmmImpute(y, X, 5)
                    
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
    y::AbstractVector = y,
    X::DataFrame = X,
    donors::Int = 5
    )

    for i in names(X)
        if X[:, i] isa CategoricalArray
            xMapping = Dict(levels(X[:, i])[j] => j for j in eachindex(levels(X[:, i])))
            X[:, i] = [xMapping[v] for v in X[:, i]]
        end
    end

    X = hcat(repeat([1], size(X, 1)), Matrix(X))

    Xₒ = X[ismissing.(y) .== 0, :]
    Xₘ = X[ismissing.(y) .== 1, :]
    yₒ = y[ismissing.(y) .== 0]

    β̂, β̇, = blrDraw(yₒ, Xₒ, 0.0001)

    ŷₒ = Xₒ * β̂
    ŷₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, 5)

    return yₒ[indices]    
end

function pmmImpute(
    y::CategoricalArray = y,
    X::DataFrame = X,
    donors::Int = 5
    )

    for i in names(X)
        if X[:, i] isa CategoricalArray
            xMapping = Dict(levels(X[:, i])[j] => j for j in eachindex(levels(X[:, i])))
            X[:, i] = [xMapping[v] for v in X[:, i]]
        end
    end

    X = hcat(repeat([1], size(X, 1)), Matrix(X))

    mapping = Dict(levels(y)[i] => i for i in eachindex(levels(y)))
    yNum = [mapping[v] for v in y]

    Xₒ = X[ismissing.(y) .== 0, :]
    Xₘ = X[ismissing.(y) .== 1, :]
    yₒ = yNum[ismissing.(y) .== 0]

    β̂, β̇, = blrDraw(yₒ, Xₒ, 0.0001)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, 5)

    return y[ismissing.(y) .== 0][indices]
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

function matchIndex(
    ŷₒ::Vector{Real} = ŷₒ, 
    ẏₘ::Vector{Real} = ẏₘ,
    donors::Int = 5
    )

    # Shuffle records to remove effects of ties
    nₒ = length(ŷₒ)
    ishuf = randperm(nₒ)
    yshuf = ŷₒ[ishuf]

    # Obtain sorting order on shuffled data
    isort = sortperm(yshuf)

    # Calculate index on input data and sort
    id = ishuf[isort]
    ysort = ŷₒ[id]

    # Pre-sample n0 values between 1 and k
    nₘ = length(ẏₘ)
    donors = min(donors, nₒ)
    donors = max(donors, 1)
    presample = sample(1:donors, nₘ, replace = true)

    indices = similar(ẏₘ, Int)

    # Loop over the target units
    for i in eachindex(ẏₘ)
        value = ẏₘ[i]
        donorID = presample[i]
        count = 0

        # Find the two adjacent neighbours
        r = searchsortedfirst(ysort, value)
        l = r - 1

        # Find the h_i'th nearest neighbour
        # Store the index of that neighbour
        while count < donorID && l >= 1 && r <= nₒ
            if value - ysort[l] < ysort[r] - value
                indices[i] = id[l]
                l -= 1
            else
                indices[i] = id[r]
                r += 1
            end
            count += 1
        end

        # If right side is exhausted, take left elements
        while count < donorID && l >= 1
            indices[i] = id[l]
            l -= 1
            count += 1
        end

        # If left side is exhausted, take right elements
        while count < donorID && r <= nₒ
            indices[i] = id[r]
            r += 1
            count += 1
        end
    end

    return indices
end

export makeMonotoneSequence, makeMethods, makePredictorMatrix, initialiseImputations, initialiseTraces, sampler, pmmImpute, blrDraw, matchIndex