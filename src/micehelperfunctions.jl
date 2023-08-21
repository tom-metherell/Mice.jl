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
    visitSequence::Vector{String},
    methods::Vector{String}
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    for i in eachindex(visitSequence)
        yVar = visitSequence[i]
        if methods[yVar] != ""
            y = data[:, yVar]
            yMissings = ismissing.(data[:, yVar])
            missingDataCount = sum(yMissings)
            imputations[i] = Matrix{nonmissingtype(eltype(relevantData))}(undef, missingDataCount, m)
            if !all(yMissings)
                for j in 1:m
                    imputations[i][:, j] = sample(y[.!yMissings], missingDataCount)
                end
            else
                if y isa CategoricalArray
                    for j in 1:m
                        imputations[i][:, j] = CategoricalArray{nonmissingtype(eltype(relevantData))}(sample(levels(y), length(y)))
                    end
                else
                    for j in 1:m
                        imputations[i][:, j] .= randn(length(y))
                    end
                end
            end
        end
    end

    return imputations
end

function initialiseTraces(
    visitSequence::Vector{String},
    iter::Int,
    m::Int
    )

    traces = [Matrix{AbstractFloat}(undef, iter, m) for _ = eachindex(visitSequence)]

    return traces
end

function sampler!(
    imputations::Vector{Matrix},
    meanTraces::Vector{Matrix{AbstractFloat}},
    varTraces::Vector{Matrix{AbstractFloat}},
    data::DataFrame,
    m::Int,
    methods::Vector{String},
    predictorMatrix::Matrix{Bool},
    iterCounter::Int,
    i::Int;
    kwargs...
    )

    yVar = visitSequence[i]
    y = data[:, yVar]
    predictorVector = predictorMatrix[:, yVar]
    predictors = names(predictorMatrix)[2][predictorVector]
    if length(predictors) > 0
        X = data[:, predictors]
        pacify!(X, predictors)
        if methods[yVar] == "pmm" && any(ismissing.(y))
            for j in 1:m
                fillXMissings!(X, predictors, visitSequence, imputations, j)
                
                imputedData = pmmImpute(y, X, 5, 1e-5)

                updateTraces!(meanTraces, varTraces, data, yVar, imputedData, i, iterCounter, j)

                imputations[i][:, j] = imputedData

                if(progressReports)
                    progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                    miceEmojis = string(repeat("🐁", floor(Int8, progress/10)), repeat("🐭", ceil(Int8, (100 - progress)/10)))
                    @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress
                end
            end
        end
    end
end

function pacify!(
    X::DataFrame,
    predictors::Vector{String}
    )

    for p in predictors
        if X[:, p] isa CategoricalArray || nonmissingtype(eltype(X[:, p])) <: AbstractString
            x = X[:, p]
            position = findfirst(names(X) .== p)
            select!(X, Not(p))
            xLevels = levels(x)
            [insertcols!(X, position+q-2, p * string(xLevels[q]) => Vector{Float64}(x .== xLevels[q])) for q in eachindex(xLevels)[2:end]]
        end
    end
end

function pacify(y::Vector)
    yLevels = levels(y)
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))
    for q in eachindex(yLevels)
        yDummies[:, q] = y .== yLevels[q]
    end

    return yDummies
end

function pacify(y::CategoricalArray)
    yLevels = levels(y)
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))
    for q in eachindex(yLevels)
        yDummies[:, q] = y .== yLevels[q]
    end

    return yDummies
end

function fillXMissings!(
        X::DataFrame,
        predictors::Vector{String},
        visitSequence::Vector{String},
        imputations::Vector{Matrix},
        j::Int
    )

    for k in predictors
        kVS = findfirst(visitSequence .== k)
        xMissings = ismissing.(X[:, k])
        if any(xMissings)
            X[xMissings, k] = imputations[kVS][:, j]
        end
    end
end

function updateTraces!(
        meanTraces::Vector{Matrix{AbstractFloat}},
        varTraces::Vector{Matrix{AbstractFloat}},
        data::DataFrame,
        yVar::String,
        imputedData::Vector,
        i::Int,
        iterCounter::Int,
        j::Int
    )

    plottingData = deepcopy(data[:, yVar])
    plottingData[ismissing.(plottingData)] = imputedData

    if plottingData isa CategoricalArray || nonmissingtype(eltype(plottingData)) <: AbstractString
        mapping = Dict(levels(plottingData)[i] => i-1 for i in eachindex(levels(plottingData)))
        plottingData = [mapping[v] for v in plottingData]
    end

    meanTraces[i][iterCounter, j] = mean(plottingData)
    varTraces[i][iterCounter, j] = var(plottingData)
end

function pmmImpute(
    y::Vector,
    X::DataFrame,
    donors::Int,
    ridge::AbstractFloat
    )

    yMissings = ismissing.(y)

    Xₒ = Matrix{Float64}(hcat(repeat[1], sum(.!yMissings), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat[1], sum(yMissings), X[yMissings, :]))

    if nonmissingtype(eltype(y)) <: AbstractString
        yₒ = y[.!yMissings]
        quantify!(yₒ, Xₒ)
    else
        yₒ = y[.!yMissings]
    end

    β̂, β̇ = blrDraw(yₒ, Xₒ, ridge)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]    
end

function pmmImpute(
    y::CategoricalArray,
    X::DataFrame,
    donors::Int,
    ridge::AbstractFloat
    )

    yMissings = ismissing.(y)

    Xₒ = Matrix{Float64}(hcat(repeat[1], sum(.!yMissings), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat[1], sum(yMissings), X[yMissings, :]))

    yₒ = y[.!yMissings]
    quantify!(yₒ, Xₒ)

    β̂, β̇ = blrDraw(yₒ, Xₒ, ridge)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]
end

####### QUANTIFY! FUNCTIONS NOT WORKING!!!!!!!

function quantify!(
    yₒ::Vector,
    Xₒ::Matrix
    )

    yDummies = pacify(yₒ)
    cca = fit(CCA, transpose(Xₒ), transpose(yDummies))
    yₒ = 

function quantify!(
    yₒ::CategoricalArray,
    Xₒ::Matrix
    )

    yDummies = pacify(yₒ)
    cca = fit(CCA, Xₒ, yDummies)

function blrDraw(
    yₒ::Vector,
    Xₒ::Matrix, 
    κ::AbstractFloat
    )

    β̂ = Xₒ \ yₒ 
    R = qr(Xₒ).R

    V = try
        inv(transpose(R) * R)
    catch
        S = transpose(R) * R;
        inv(S + diagm(diag(S)) * κ)
    end

    σ̇ = sqrt(sum((yₒ - Xₒ * β̂).^2)) / rand(Chisq(max(length(yₒ) - size(Xₒ, 2), 1)))
    β̇ = β̂ + σ̇ * cholesky((V + transpose(V)) / 2).factors * randn(size(Xₒ, 2))

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
    selections = sample(1:donors, nₘ, replace = true)

    indices = similar(ẏₘ, Int)

    # Loop over the target units
    for i in eachindex(ẏₘ)
        value = ẏₘ[i]
        donorID = selections[i]
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

export makeMonotoneSequence, makeMethods, makePredictorMatrix, initialiseImputations, sampler!