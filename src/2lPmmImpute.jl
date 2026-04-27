# The twoLevelPmmImpute! function includes a ! as it updates loggedEvents in place
function twoLevelPmmImpute!(
    yData::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    types::Vector{Int},
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    intercept::Bool = true,
    nIterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    matchSampledPars::Bool = true,
    noise::Float64 = 1e5,
    unusedKwargs...
    )

    if whereCount == 0
        return Vector{eltype(yData)}(undef, 0)
    end

    if intercept
        X = hcat(ones(Float64, size(X, 1)), X)
        types = vcat(2, types)
    end

    classCols = findall(types .== -2)
    randomCols = findall(types .== 2)

    if isempty(classCols)
        throw(ArgumentError("Two-level PMM method specified, but no class variable (coded -2) found."))
    end
    if isempty(randomCols)
        throw(ArgumentError("Two-level PMM requires at least one random predictor (coded 2)."))
    end

    classMatrix = Matrix{Float64}(X[:, classCols])
    classKeys = [Tuple(classMatrix[r, c] for c in axes(classMatrix, 2)) for r in axes(classMatrix, 1)]
    classLevels = unique(classKeys)
    nClasses = length(classLevels)
    classMap = Dict(level => idx for (idx, level) in enumerate(classLevels))
    gfFull = [classMap[key] for key in classKeys]
    gf = gfFull[.!whereY]

    yₒRaw = yData[.!whereY]
    yₒ = Vector{Float64}(yₒRaw)

    if length(unique(gf)) < nClasses
        throw(ArgumentError("Two-level PMM requires at least one observed outcome per class."))
    end

    Xₒ = Matrix{Float64}(X[.!whereY, randomCols])
    Xₘ = Matrix{Float64}(X[whereY, randomCols])

    XG = [Matrix{Float64}(Xₒ[gf .== class, :]) for class in 1:nClasses]
    Xₛₛ = [transpose(A) * A for A in XG]
    yG = [yₒ[gf .== class] for class in 1:nClasses]
    nRc = size(XG[1], 2)

    β = zeros(nClasses, nRc)
    ss = zeros(nClasses)
    μ = zeros(nRc)
    invΨ = Matrix{Float64}(Diagonal(ones(nRc)))
    invσ² = ones(nClasses)
    σ²₀ = 1.0
    θ = 1.0

    for _ in 1:nIterations
        for class in 1:nClasses
            vv = symridge(invσ²[class] * Xₛₛ[class] + invΨ, ridge)
            βVar = chol2inv(cholesky(vv))

            meanPart = βVar * (
                transpose(invσ²[class] * XG[class]) * yG[class] + invΨ * μ
            )

            noisePart = cholesky(symridge(βVar, ridge)).U' * randn(nRc)
            β[class, :] = meanPart + noisePart

            resid = yG[class] - XG[class] * β[class, :]
            ss[class] = sum(resid .^ 2)
        end

        μVar = chol2inv(cholesky(symridge(invΨ, ridge))) / nClasses
        μ = vec(mean(β, dims = 1)) + cholesky(μVar).U' * randn(nRc)

        centred = β .- transpose(μ)
        psiScale = chol2inv(cholesky(symridge(transpose(centred) * centred, ridge)))
        invΨ = randWishart(nClasses - nRc - 1, cholesky(psiScale).U)

        for class in 1:nClasses
            shape = sum(gf .== class) / 2 + 1 / (2 * θ)
            scale = 2 * θ / (ss[class] * θ + σ²₀)
            invσ²[class] = rand(Gamma(shape, scale))
        end

        H = 1 / mean(invσ²)
        σ²₀ = rand(Gamma(nClasses / (2 * θ) + 1, 2 * θ * H / nClasses))

        G = exp(mean(log.(1 ./ invσ²)))
        thetaShape = nClasses / 2 - 1
        thetaScale = 2 / (nClasses * (σ²₀ / H - log(σ²₀) + log(G) - 1))
        θ = 1 / rand(Gamma(thetaShape, thetaScale))
    end

    gfₘ = gfFull[whereY]

    # Predicted values for missing units based on sampled parameters.
    ŷₘ = [dot(Xₘ[r, :], β[gfₘ[r], :]) for r in eachindex(gfₘ)]

    # miceadds uses non-sampled parameters for observed units unless explicitly overridden.
    βₒ = matchSampledPars ? [vec(β[c, :]) for c in 1:nClasses] : [ridgeCoef(XG[c], yG[c], ridge) for c in 1:nClasses]
    ŷₒ = [dot(Xₒ[r, :], βₒ[gf[r]]) for r in eachindex(gf)]

    donorIdx = pmm5DonorIndices(ŷₒ, ŷₘ, donors; noise = noise)
    return yₒRaw[donorIdx]
end

function ridgeCoef(X::Matrix{Float64}, y::Vector{Float64}, ridge::Float64)
    nPred = size(X, 2)
    S = transpose(X) * X
    return (S + ridge * I(nPred)) \ (transpose(X) * y)
end

function pmm5DonorIndices(
    ŷₒ::Vector{Float64},
    ŷₘ::Vector{Float64},
    donors::Int;
    noise::Float64 = 1e5
    )

    nObs = length(ŷₒ)
    nMis = length(ŷₘ)

    if nObs == 0
        throw(ArgumentError("Cannot do PMM without observed donor values."))
    end

    donors = max(min(donors, nObs), 1)

    # rows: [obs_flag, original_index_within_obs_or_mis, yhat]
    obsRows = [Float64[1.0, i, ŷₒ[i]] for i in 1:nObs]
    misRows = [Float64[0.0, i, ŷₘ[i]] for i in 1:nMis]
    rows = vcat(obsRows, misRows)

    # Add tiny noise to break ties, following miceadds pmm5.
    yHats = [row[3] for row in rows]
    diffs = abs.(diff(yHats))
    posDiffs = diffs[diffs .> 0]
    fg1 = isempty(posDiffs) ? 1.0 : minimum(posDiffs)
    for row in rows
        row[3] += rand() * fg1 / noise
    end

    sortedRows = sort(rows, by = row -> row[3])
    n = length(sortedRows)

    obsCum = cumsum([Int(row[1]) for row in sortedRows])
    nY = nObs

    obsLow = clamp.(obsCum, 1, nY)
    obsRev = reverse([Int(row[1]) for row in sortedRows])
    c1 = nY .- cumsum(obsRev) .+ 1
    obsUpp = clamp.(reverse(c1), 1, nY)

    sortedMis = [sortedRows[idx] for idx in 1:n if sortedRows[idx][1] == 0.0]
    misPos = findall(row -> row[1] == 0.0, sortedRows)

    # Keep missing rows in their original missing-order.
    sortpermMis = sortperm([Int(row[2]) for row in sortedMis])
    misPos = misPos[sortpermMis]

    obsIndices = [Int(sortedRows[idx][2]) for idx in 1:n if sortedRows[idx][1] == 1.0]

    donorMat = Matrix{Int}(undef, nMis, 2 * donors)

    for dd in 1:donors
        lowIdx = clamp.([obsLow[pos] - dd + 1 for pos in misPos], 1, nY)
        uppIdx = clamp.([obsUpp[pos] + dd - 1 for pos in misPos], 1, nY)

        donorMat[:, dd] = obsIndices[lowIdx]
        donorMat[:, dd + donors] = obsIndices[uppIdx]
    end

    sampledCol = rand(1:(2 * donors), nMis)
    return [donorMat[i, sampledCol[i]] for i in 1:nMis]
end

const TWO_LEVEL_PMM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents;
    intercept::Bool = true,
    nIterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    matchSampledPars::Bool = true,
    noise::Float64 = 1e5,
    kwargs...) -> begin
    twoLevelPmmImpute!(
        yData,
        X,
        whereY,
        whereCount,
        types,
        yVar,
        iterCounter,
        j,
        loggedEvents;
        intercept = intercept,
        nIterations = nIterations,
        ridge = ridge,
        donors = donors,
        matchSampledPars = matchSampledPars,
        noise = noise,
        kwargs...
    )
end; twoLevel = true)

registerImputer!("2l.pmm", TWO_LEVEL_PMM_IMPUTER)
