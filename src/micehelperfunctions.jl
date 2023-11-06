function makeMonotoneSequence(data::DataFrame)
    missingness = Vector{Int}(undef, size(data, 2))

    # Counting missing data in each column
    for i in axes(data, 2)
        missingness[i] = sum(ismissing.(data[:, i]))
    end

    # Sort the missingness vector in ascending order
    missingness = sortperm(missingness)

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
    methods::NamedVector{String}
    )

    imputations = Vector{Matrix}(undef, ncol(data))

    for i in eachindex(visitSequence)
        yVar = visitSequence[i]
        if methods[yVar] != ""
            y = data[:, yVar]
            yMissings = ismissing.(data[:, yVar])
            missingDataCount = sum(yMissings)
            imputations[i] = Matrix{nonmissingtype(eltype(y))}(undef, missingDataCount, m)
            if !all(yMissings)
                for j in 1:m
                    imputations[i][:, j] = sample(y[.!yMissings], missingDataCount)
                end
            else
                if y isa CategoricalArray
                    for j in 1:m
                        imputations[i][:, j] = CategoricalArray{nonmissingtype(eltype(y))}(sample(levels(y), length(y)))
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

const initializeImputations = initialiseImputations

function initialiseTraces(
    visitSequence::Vector{String},
    iter::Int,
    m::Int
    )

    traces = [Matrix{Float64}(undef, iter, m) for _ = eachindex(visitSequence)]

    return traces
end

const initializeTraces = initialiseTraces

function sampler!(
    imputations::Vector{Matrix},
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    data::DataFrame,
    m::Int,
    visitSequence::Vector{String},
    methods::NamedVector{String},
    predictorMatrix::NamedArray{Bool},
    iter::Int,
    iterCounter::Int,
    i::Int,
    progressReports::Bool,
    loggedEvents::Vector{String},
    threads::Bool;
    kwargs...
    )

    yVar = visitSequence[i]
    y = data[:, yVar]
    predictorVector = predictorMatrix[:, yVar]
    predictors = names(predictorVector)[1][predictorVector]
    if length(predictors) > 0
        if methods[yVar] == "pmm" && any(ismissing.(y))
            if threads
                lk = ReentrantLock()

                Threads.@threads for j in 1:m
                    tempLog = Vector{String}([])
                    X = data[:, predictors]
                    fillXMissings!(X, predictors, visitSequence, imputations, j)
                    X = pacify(X, predictors)
                    origNCol = size(X, 2)
                    removeLinDeps!(X, y)

                    if size(X, 2) > 0
                        if size(X, 2) < origNCol
                            diff = OrigNCol - size(X, 2)
                            push!(tempLog, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end
                        imputedData = pmmImpute!(y, X, 5, 1e-5, yVar, iterCounter, j, tempLog)
                    else
                        push!(tempLog, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                        imputedData = imputations[i][:, j]
                    end

                    lock(lk)
                    try
                        updateTraces!(meanTraces, varTraces, imputedData, i, iterCounter, j)
                        imputations[i][:, j] = imputedData
                        append!(loggedEvents, tempLog)
                    finally
                        unlock(lk)
                    end                    
                    
                    if progressReports
                        progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter) * 100
                        progressRound = floor(Int8, progress / 10)
                        miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                        @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar miceEmojis progress length(loggedEvents)
                    end
                end
            else
                for j in 1:m
                    X = data[:, predictors]
                    fillXMissings!(X, predictors, visitSequence, imputations, j)
                    X = pacify(X, predictors)
                    origNCol = size(X, 2)
                    removeLinDeps!(X, y)

                    if size(X, 2) > 0
                        if size(X, 2) < origNCol
                            diff = origNCol - size(X, 2)
                            push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end
                        imputedData = pmmImpute!(y, X, 5, 1e-5, yVar, iterCounter, j, loggedEvents)
                    else
                        push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                        imputedData = imputations[i][:, j]
                    end

                    updateTraces!(meanTraces, varTraces, imputedData, i, iterCounter, j)

                    imputations[i][:, j] = imputedData

                    if progressReports
                        progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                        progressRound = floor(Int8, progress / 10)
                        miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                        @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress length(loggedEvents)
                    end
                end
            end
        else
            if methods[yVar] == ""
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no method specified.")
            elseif methods[yVar] != "pmm"
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - method not supported.")
            else
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no missing data.")
            end
        end
    end
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

function pacify(
    X::DataFrame,
    predictors::Vector{String}
    )

    categoricalPredictors = Vector{String}([])

    for xVar in predictors
        x = X[:, xVar]
        if x isa CategoricalArray || nonmissingtype(eltype(x)) <: AbstractString
            if length(levels(x)) > 2
                push!(categoricalPredictors, xVar)
            end
        elseif nonmissingtype(eltype(x)) <: Real
            X[!, xVar] = convert.(Float64, X[:, xVar])
            X[:, xVar] = standardize(UnitRangeTransform, X[:, xVar])
        end
    end

    mf = ModelFrame(term(0) ~ sum(term.(predictors)), X)
    setcontrasts!(mf, Dict([Symbol(xVar) => PolynomialCoding() for xVar in categoricalPredictors]))

    X = ModelMatrix(mf).m[:, 2:end]

    return X
end

mutable struct PolynomialCoding <: AbstractContrasts
end

import StatsModels.contrasts_matrix

function contrasts_matrix(C::PolynomialCoding, _, n)
    X = reduce(hcat, [((1:n) .- mean(1:n)) .^ i for i in 0:n-1])
    qrX = qr(X)
    Z = qrX.Q * Diagonal(qrX.R)
    for i in axes(Z, 2)
        Z[:, i] = Z[:, i] ./ sqrt(sum(Z[:, i].^2))
    end
    return Z[:, 2:end]
end

import StatsModels.termnames

function termnames(C::PolynomialCoding, levels::AbstractVector, _::Integer)
    return Vector{String}([".^$i" for i in 1:length(levels)])
end

function pacify(y::AbstractArray)
    yLevels = levels(y)
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))
    for q in eachindex(yLevels)
        yDummies[:, q] = y .== yLevels[q]
    end

    return yDummies
end

function removeLinDeps!(
    X::Matrix{Float64},
    y::AbstractArray
    )

    if all(ismissing.(y))
        return
    end

    Xₒ = Matrix{Float64}(X[.!ismissing.(y), :])
    
    if y isa CategoricalArray || nonmissingtype(eltype(y)) <: AbstractString
        yₒ = y[.!ismissing.(y)]
        mapping = Dict(levels(yₒ)[i] => i-1 for i in eachindex(levels(yₒ)))
        yₒ = Vector{Float64}([mapping[v] for v in yₒ])
    else
        yₒ = Vector{Float64}(y[.!ismissing.(y)])
    end

    if var(yₒ) < 1e-4
        select!(X, [])
        return
    end

    keep = var.(eachcol(Xₒ)) .> 1e-4 .&& cor.(eachcol(Xₒ), [yₒ]) .< 0.99

    keepSum = sum(keep)
    if keepSum < 2
        X = X[:, keep]
        return
    end

    xCors = cor(Xₒ)
    nxCors = xCors[findall(keep), findall(keep)]
    eigenCors = eigen(nxCors)

    eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
    sortedeigvals = eigenCors.values[eigvalsorder]
    sortedeigvecs = eigenCors.vectors[:, eigvalsorder]

    while sortedeigvals[keepSum] / sortedeigvals[1] < 1e-4
        w = sortperm(abs.(sortedeigvecs[:, keepSum]), rev = true)[1]
        keep[findall(keep)[w]] = false
        nxCors = xCors[findall(keep), findall(keep)]
        keepSum -= 1
        eigenCors = eigen(nxCors)
        eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
        sortedeigvals = eigenCors.values[eigvalsorder]
        sortedeigvecs = eigenCors.vectors[:, eigvalsorder]
    end

    X = X[:, keep]
end

function updateTraces!(
        meanTraces::Vector{Matrix{Float64}},
        varTraces::Vector{Matrix{Float64}},
        imputedData::AbstractArray,
        i::Int,
        iterCounter::Int,
        j::Int
    )

    if imputedData isa CategoricalArray || nonmissingtype(eltype(imputedData)) <: AbstractString
        mapping = Dict(levels(imputedData)[i] => i-1 for i in eachindex(levels(imputedData)))
        imputedData = [mapping[v] for v in imputedData]
    end

    meanTraces[i][iterCounter, j] = mean(imputedData)
    varTraces[i][iterCounter, j] = var(imputedData)
end

function pmmImpute!(
    y::AbstractArray,
    X::Matrix{Float64},
    donors::Int,
    ridge::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    yMissings = ismissing.(y)

    Xₒ = Matrix{Float64}(hcat(repeat([1], sum(.!yMissings)), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], sum(yMissings)), X[yMissings, :]))

    if nonmissingtype(eltype(y)) <: AbstractString
        mapping = Dict(levels(y[.!yMissings])[i] => i-1 for i in eachindex(levels(y[.!yMissings])))
        yₒ = Vector{Float64}([mapping[v] for v in y[.!yMissings]])
        # yₒ = quantify(y[.!yMissings], Xₒ)
    else
        yₒ = Vector{Float64}(y[.!yMissings])
    end

    β̂, β̇ = blrDraw!(yₒ, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]    
end

function pmmImpute!(
    y::CategoricalArray,
    X::Matrix{Float64},
    donors::Int,
    ridge::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    yMissings = ismissing.(y)

    Xₒ = Matrix{Float64}(hcat(repeat([1], sum(.!yMissings)), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], sum(yMissings)), X[yMissings, :]))

    yₒ = quantify(y[.!yMissings], Xₒ)

    β̂, β̇ = blrDraw!(yₒ, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]
end

function quantify(
    yₒ::AbstractArray,
    Xₒ::Matrix{Float64}
    )

    yDummies = pacify(yₒ)
    Ycoef = miceCCA(Xₒ, yDummies)
    yₒ = Vector{Float64}(zscore(yDummies * Ycoef[:, 2]))

    return yₒ
end

function miceCCA(
    X::Matrix{Float64},
    Y::Matrix{Float64}
    )

    qrX = qr(X)
    qrY = qr(Y)
    nr = size(X, 1)
    dx = rank(X)
    dy = rank(Y)

    Z = svd((transpose(qrY.Q) * (qrX.Q * diagm(nr, dx, repeat([1], min(nr, dx)))))[1:dy, :])

    Ycoef = qrY.R \ Z.U

    return Ycoef
end

function blrDraw!(
    yₒ::Vector{Float64},
    Xₒ::Matrix{Float64}, 
    κ::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    β̂ = Xₒ \ yₒ 
    R = qr(Xₒ).R

    V = try
        inv(transpose(R) * R)
    catch
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: ridge penalty applied (unstable results) - predictors are highly multicollinear.")
        S = transpose(R) * R;
        inv(S + Diagonal(S) * κ)
    end

    σ̇ = sqrt(sum((yₒ - Xₒ * β̂).^2)) / rand(Chisq(max(length(yₒ) - size(Xₒ, 2), 1)))
    β̇ = β̂ + σ̇ * cholesky((V + transpose(V)) / 2).factors * randn(size(Xₒ, 2))

    return β̂, β̇
end

function matchIndex(
    ŷₒ::Vector{Float64}, 
    ẏₘ::Vector{Float64},
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