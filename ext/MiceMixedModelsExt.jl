module MiceMixedModelsExt
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Distributions
    using LinearAlgebra: Diagonal, Symmetric, cholesky, dot, eigen, inv, transpose
    using Mice: Imputer, makeMethods, makePredictorMatrix, mice, randWishart, registerImputer!, symridge
    using MixedModels: GeneralizedLinearMixedModel, fit!, fixef, ranef
    using PrecompileTools: @compile_workload
    using Random: rand, randn
    using Statistics: var
    import StatsAPI
    import StatsModels

    # The twoLevelBinImpute! function includes a ! as it updates loggedEvents in place
    function twoLevelBinImpute!(
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
        unusedKwargs...
        )

        if whereCount == 0
            return Vector{eltype(yData)}(undef, 0)
        end

        # Keep keyword for API compatibility.
        _ = nIterations

        if intercept
            X = hcat(ones(Float64, size(X, 1)), X)
            types = vcat(2, types)
        end

        classCols = findall(types .== -2)
        randomCols = findall(types .== 2)
        fixedCols = findall(types .> 0)

        if isempty(classCols)
            throw(ArgumentError("Two-level binary method specified, but no class variable (coded -2) found."))
        end
        if isempty(randomCols)
            throw(ArgumentError("Two-level binary requires at least one random predictor (coded 2)."))
        end

        classMatrix = Matrix{Float64}(X[:, classCols])
        classKeys = [Tuple(classMatrix[r, c] for c in axes(classMatrix, 2)) for r in axes(classMatrix, 1)]
        classLevels = unique(classKeys)
        nClasses = length(classLevels)
        classMap = Dict(level => idx for (idx, level) in enumerate(classLevels))
        gfFull = [classMap[key] for key in classKeys]
        gf = gfFull[.!whereY]

        yₒRaw = yData[.!whereY]
        yₒ = Vector{Float64}([y isa CategoricalValue ? y.ref - 1 : convert(Float64, y) for y in yₒRaw])

        if length(unique(gf)) < nClasses
            throw(ArgumentError("Two-level binary requires at least one observed outcome per class."))
        end

        # Check that y is binary
        uniqueVals = unique(yₒ)
        if length(uniqueVals) != 2 || !all(v ∈ [0.0, 1.0] for v in uniqueVals)
            throw(ArgumentError("Two-level binary imputation requires a binary outcome variable. Got unique values: $uniqueVals"))
        end

        yObs = Int.(round.(yₒ))
        obsMask = .!whereY

        Xₒ = Matrix{Float64}(X[obsMask, fixedCols])
        Zₒ = Matrix{Float64}(X[obsMask, randomCols])
        Xₘ = Matrix{Float64}(X[whereY, fixedCols])
        Zₘ = Matrix{Float64}(X[whereY, randomCols])

        if size(Zₒ, 2) == 0
            throw(ArgumentError("Two-level binary requires at least one random-effect column in the design matrix."))
        end

        randColsNoIntercept = [k for k in 1:size(Zₒ, 2) if !all(abs.(Zₒ[:, k] .- 1.0) .< 1e-12)]

        colSyms = Symbol[:y, :cluster]
        colData = Vector{Any}[yObs, gf]

        fixedSyms = Symbol[]
        for k in 1:size(Xₒ, 2)
            s = Symbol("x", k)
            push!(fixedSyms, s)
            push!(colSyms, s)
            push!(colData, Xₒ[:, k])
        end

        randomSyms = Symbol[]
        for k in randColsNoIntercept
            s = Symbol("z", k)
            push!(randomSyms, s)
            push!(colSyms, s)
            push!(colData, Zₒ[:, k])
        end

        tableData = NamedTuple{Tuple(colSyms)}(Tuple(colData))

        fixedPart = isempty(fixedSyms) ? "1" : "1 + " * join(string.(fixedSyms), " + ")
        randomPart = isempty(randomSyms) ? "(1 | cluster)" : "(1 + " * join(string.(randomSyms), " + ") * " | cluster)"
        formulaString = "y ~ " * fixedPart * " + " * randomPart
        formulaExpr = Meta.parse(formulaString)
        f = Core.eval(@__MODULE__, :(StatsModels.@formula($formulaExpr)))

        model = nothing
        try
            model = fit!(GeneralizedLinearMixedModel(f, tableData, Distributions.Bernoulli()))
        catch err
            push!(loggedEvents, "2l.bin: MixedModels fit failed for $(yVar) (iter $(iterCounter), chain $(j)): $(err)")
            return yData[whereY]
        end

        βhat = fixef(model)
        βcov = Matrix{Float64}(StatsAPI.vcov(model))
        βcov = symridge(βcov, ridge)
        βdraw = βhat + cholesky(βcov).U' * randn(length(βhat))

        rancoef = Matrix{Float64}(ranef(model)[1])
        q = size(rancoef, 2)
        ψhat = q == 1 ? reshape(var(vec(rancoef)), 1, 1) : Statistics.cov(rancoef)
        ψhat = symridge(Matrix{Float64}(ψhat), ridge)

        λ = transpose(rancoef) * rancoef
        s = q * ψhat
        temp = symridge(λ + s, ridge)
        ev = eigen(Symmetric(temp))
        evValues = max.(ev.values, 0.0)
        deco = ev.vectors * Diagonal(sqrt.(evValues))

        ν = size(rancoef, 1) + q
        ψstarInvScale = randWishart(ν, Diagonal(ones(q)))
        ψstar = inv(symridge(deco * ψstarInvScale * transpose(deco), ridge))
        ψstar = symridge(ψstar, ridge)

        # Impute missing values using final coefficients
        gfMissing = gfFull[whereY]
        imputedValues = Vector{Float64}(undef, whereCount)
        randEffectDraws = Dict{Int, Vector{Float64}}()

        for i in 1:whereCount
            classIdx = gfMissing[i]
            if !haskey(randEffectDraws, classIdx)
                randEffectDraws[classIdx] = vec(rand(Distributions.MvNormal(zeros(q), ψstar)))
            end
            η = dot(Xₘ[i, :], βdraw) + dot(Zₘ[i, :], randEffectDraws[classIdx])
            p = 1.0 / (1.0 + exp(-η))
            imputedValues[i] = rand() < p ? 1.0 : 0.0
        end

        # Convert back to original data type
        if eltype(yₒRaw) <: CategoricalValue
            levels_list = levels(yₒRaw[1])
            return CategoricalArray([levels_list[Int(val) + 1] for val in imputedValues], levels_list)
        else
            return imputedValues
        end
    end

    const TWO_LEVEL_BIN_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept::Bool = true, nIterations::Int = 100, ridge::Float64 = 1e-4, kwargs...) -> begin
        twoLevelBinImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept = intercept, nIterations = nIterations, ridge = ridge, kwargs...)
    end; twoLevel = true)

    function __init__()
        registerImputer!("2l.bin", TWO_LEVEL_BIN_IMPUTER)
    end

    @compile_workload begin
        n = 48
        cluster = repeat(1:8, inner = 6)
        x = randn(n)
        η = @. -0.2 + 0.8 * x + 0.1 * cluster
        y = Vector{Union{Missing, Int}}(Int.(rand(n) .< (1.0 ./ (1.0 .+ exp.(-η)))))

        y[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

        ct = (
            y = y,
            x = x,
            cluster = cluster
        )

        methods = makeMethods(ct)
        methods .= ""
        methods["y"] = "2l.bin"

        predictorMatrix = makePredictorMatrix(ct)
        predictorMatrix[:, :] .= 0
        predictorMatrix["y", "x"] = 2
        predictorMatrix["y", "cluster"] = -2

        mice(ct, m = 1, iter = 1, methods = methods, predictorMatrix = predictorMatrix, progressReports = false, nIterations = 1, ridge = 1e-3)
    end
end
