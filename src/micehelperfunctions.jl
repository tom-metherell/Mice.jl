function makeMonotoneSequence(data::DataFrame)
    missingness = Vector{Int}(undef, size(data, 2))

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

function makeMethods(data::DataFrame)
    methods = NamedArray(Vector{String}(fill("pmm", ncol(data))))

    setnames!(methods, names(data), 1)

    return methods
end

function makePredictorMatrix(data::DataFrame)
    predictorMatrix = NamedArray(Matrix{Bool}(fill(1, ncol(data), ncol(data))))
    for i in 1:ncol(data)
        predictorMatrix[i, i] = 0
    end

    setnames!(predictorMatrix, names(data), 1)
    setnames!(predictorMatrix, names(data), 2)

    return predictorMatrix
end

function initialiseImputations(
    data::DataFrame,
    m::Int,
    visitSequence::AbstractVector,
    methods::AbstractVector
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    presentData = ismissing.(data) .== 0

    for i in eachindex(visitSequence)
        yVar = visitSequence[i]
        if methods[yVar] != ""
            relevantData = data[:, yVar]
            presentLocations = BitVector(presentData[:, yVar])
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
    visitSequence::AbstractVector,
    iter::Int,
    m::Int
    )

    traces = [Matrix{AbstractFloat}(undef, iter, m) for _ = eachindex(visitSequence)]

    return traces
end

function pmmImpute(
    y::AbstractVector,
    X::DataFrame,
    donors::Int
    )

    for z in axes(X, 2)
        if X[:, z] isa CategoricalArray
            name = names(X)[z]
            xArray = deepcopy(X[:, z])
            select!(X, Not(z))
            xMapping = Dict(levels(xArray)[j] => j for j in eachindex(levels(xArray)))
            insertcols!(X, z, name => Vector{Int64}([xMapping[v] for v in xArray]))
        end
    end

    X = hcat(repeat([1], size(X, 1)), Matrix(X))

    Xₒ = X[ismissing.(y) .== 0, :]
    Xₘ = X[ismissing.(y) .== 1, :]
    yₒ = y[ismissing.(y) .== 0]

    β̂, β̇ = Mice.blrDraw(yₒ, Xₒ, 0.0001)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return yₒ[indices]    
end

function pmmImpute(
    y::CategoricalArray,
    X::DataFrame,
    donors::Int
    )

    for z in axes(X, 2)
        if X[:, z] isa CategoricalArray
            name = names(X)[z]
            xArray = deepcopy(X[:, z])
            select!(X, Not(z))
            xMapping = Dict(levels(xArray)[j] => j for j in eachindex(levels(xArray)))
            insertcols!(X, z, name => Vector{Int}([xMapping[v] for v in xArray]))
        end
    end

    X = hcat(repeat([1], size(X, 1)), Matrix(X))

    Xₒ = X[ismissing.(y) .== 0, :]
    Xₘ = X[ismissing.(y) .== 1, :]

    yₒ = y[ismissing.(y) .== 0]
    mapping = Dict(levels(yₒ)[i] => i for i in eachindex(levels(yₒ)))
    yₒ = [mapping[v] for v in yₒ]

    β̂, β̇ = blrDraw(yₒ, Xₒ, 0.0001)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[ismissing.(y) .== 0][indices]
end

function blrDraw(
    yₒ,
    Xₒ, 
    κ::AbstractFloat
    )

    S = transpose(Xₒ) * Xₒ

    V = inv(S + diagm(diag(S)) * κ)

    β̂ = V * transpose(Xₒ) * yₒ

    ġ = rand(Chisq(size(Xₒ, 1) - size(Xₒ, 2)))

    σ̇ = sqrt((transpose(yₒ - Xₒ * β̂) * (yₒ - Xₒ * β̂)) / ġ)

    ż = randn(size(Xₒ, 2))

    sqrtV = cholesky((V + transpose(V)) / 2).factors

    β̇ = β̂ + σ̇ * sqrtV * ż 

    return β̂, β̇
end

function matchIndex(
    ŷₒ::Vector, 
    ẏₘ::Vector,
    donors::Int
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

export makeMonotoneSequence, makeMethods, makePredictorMatrix, initialiseImputations, pmmImpute, blrDraw, matchIndex