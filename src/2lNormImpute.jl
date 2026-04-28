# The twoLevelNormImpute! function includes a ! as it updates loggedEvents in place
function twoLevelNormImpute!(
    yData::Vector{Float64},
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
    unusedKwargs...
    )

    prep = prepareTwoLevelImputationInputs(X, whereY, types; intercept = intercept)
    gfFull = prep.gfFull
    gf = prep.gf
    nClasses = prep.nClasses
    yₒ = Vector{Float64}(yData[.!whereY])

    if length(unique(gf)) < nClasses
        throw(ArgumentError("Two-level imputation requires at least one observed outcome per class."))
    end

    Xₒ = prep.Xₒ
    Xₘ = prep.Xₘ
    XG = [Matrix{Float64}(Xₒ[gf .== class, :]) for class in 1:nClasses]
    Xₛₛ = [transpose(A) * A for A in XG]
    yG = [yₒ[gf .== class] for class in 1:nClasses]
    nG = [sum(gf .== class) for class in 1:nClasses]
    nRc = size(Xₒ, 2)

    β, invσ² = drawTwoLevelParams(XG, Xₛₛ, yG, nG, nIterations, ridge)

    if whereCount == 0
        return Vector{Float64}(undef, 0)
    end

    gfₘ = gfFull[whereY]
    imps = Vector{Float64}(undef, whereCount)

    for r in eachindex(gfₘ)
        class = gfₘ[r]
        imps[r] = dot(Xₘ[r, :], β[class, :]) + randn() * sqrt(1 / invσ²[class])
    end

    return imps
end

function prepareTwoLevelImputationInputs(
    X::Matrix{Float64},
    whereY::Vector{Bool},
    types::Vector{Int};
    intercept::Bool = true
    )

    XWork = X
    typesWork = types
    if intercept
        XWork = hcat(ones(Float64, size(XWork, 1)), XWork)
        typesWork = vcat(2, typesWork)
    end

    classCols = findall(typesWork .== -2)
    randomCols = findall(typesWork .== 2)
    fixedCols = findall(typesWork .> 0)

    if isempty(classCols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end
    if isempty(randomCols)
        throw(ArgumentError("Two-level imputation requires at least one random predictor (coded 2)."))
    end

    classMatrix = Matrix{Float64}(XWork[:, classCols])
    classKeys = [Tuple(classMatrix[r, c] for c in axes(classMatrix, 2)) for r in axes(classMatrix, 1)]
    classLevels = unique(classKeys)
    nClasses = length(classLevels)
    classMap = Dict(level => idx for (idx, level) in enumerate(classLevels))
    gfFull = [classMap[key] for key in classKeys]
    gf = gfFull[.!whereY]

    Xₒ = Matrix{Float64}(XWork[.!whereY, randomCols])
    Xₘ = Matrix{Float64}(XWork[whereY, randomCols])

    return (
        X = XWork,
        types = typesWork,
        randomCols = randomCols,
        fixedCols = fixedCols,
        gfFull = gfFull,
        gf = gf,
        nClasses = nClasses,
        Xₒ = Xₒ,
        Xₘ = Xₘ,
    )
end

function drawTwoLevelParams(
    XG::Vector{Matrix{Float64}},
    Xₛₛ::Vector{Matrix{Float64}},
    yG::Vector{Vector{Float64}},
    nG::Vector{Int},
    nIterations::Int,
    ridge::Float64
    )

    nClasses = length(XG)
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
            ss[class] = dot(resid, resid)
        end

        μVar = chol2inv(cholesky(symridge(invΨ, ridge))) / nClasses
        μ = vec(mean(β, dims = 1)) + cholesky(μVar).U' * randn(nRc)

        centred = β .- transpose(μ)
        ψScale = chol2inv(cholesky(symridge(centred' * centred, ridge)))
        invΨ = randWishart(nClasses - nRc - 1, cholesky(ψScale).U)

        for class in 1:nClasses
            shape = nG[class] / 2 + 1 / (2 * θ)
            scale = 2 * θ / (ss[class] * θ + σ²₀)
            invσ²[class] = rand(Gamma(shape, scale))
        end

        H = 1 / mean(invσ²)
        σ²₀ = rand(Gamma(nClasses / (2 * θ) + 1, 2 * θ * H / nClasses))

        G = exp(mean(log.(1 ./ invσ²)))
        θShape = nClasses / 2 - 1
        θScale = 2 / (nClasses * (σ²₀ / H - log(σ²₀) + log(G) - 1))
        θ = 1 / rand(Gamma(θShape, θScale))
    end

    return β, invσ²
end

function symridge(A::AbstractMatrix{Float64}, ridge::Float64 = 1e-4)
    Asym = (A + transpose(A)) ./ 2
    if size(Asym, 1) == 1
        return Matrix{Float64}(Asym)
    end
    return Matrix{Float64}(Asym + Diagonal(diag(Asym) .* ridge))
end

function chol2inv(C)
    return inv(C.U) * inv(C.U')
end

function randWishart(df::Int, sqrtΣ::AbstractMatrix{Float64})
    p = size(sqrtΣ, 1)
    Z = zeros(p, p)

    for i in 1:p
        Z[i, i] = sqrt(rand(Chisq(df - i + 1)))
    end
    for i in 2:p
        for j in 1:i-1
            Z[i, j] = randn()
        end
    end

    B = Z * sqrtΣ
    return transpose(B) * B
end

const TWO_LEVEL_NORM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept::Bool = true, nIterations::Int = 100, ridge::Float64 = 1e-4, kwargs...) -> begin
    twoLevelNormImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept = intercept, nIterations = nIterations, ridge = ridge, kwargs...)
end; twoLevel = true)

registerImputer!("2l.norm", TWO_LEVEL_NORM_IMPUTER)