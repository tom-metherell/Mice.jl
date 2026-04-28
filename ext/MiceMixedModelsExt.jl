module MiceMixedModelsExt
    using CategoricalArrays: categorical, CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Distributions
    using LinearAlgebra: Diagonal, Symmetric, cholesky, dot, eigen, inv, transpose, Hermitian, I
    using Mice: Imputer, makeMethods, makePredictorMatrix, mice, prepareTwoLevelImputationInputs, randWishart, registerImputer!, symridge
    using MixedModels: GeneralizedLinearMixedModel, fit!, fixef, ranef
    using PrecompileTools: @compile_workload
    using Random: rand, randn
    using Statistics: cov, var
    using StatsBase: countmap
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
        ridge::Float64 = 1e-4,
        unusedKwargs...
        )

        if whereCount == 0
            return Vector{eltype(yData)}(undef, 0)
        end

        if isempty(findall(types .== 2))
            throw(ArgumentError("Two-level binary requires at least one random predictor (coded 2)."))
        end

        prep = prepareTwoLevelImputationInputs(X, whereY, types; intercept = false)
        gfFull = prep.gfFull
        gf = prep.gf
        nClasses = prep.nClasses
        randomCols = prep.randomCols
        fixedCols = prep.fixedCols

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

        colSyms = Symbol[:y, :cluster]
        colData = Any[yObs, categorical(gf)]  # Changed from Vector{Any}[...]

        fixedSyms = Symbol[]
        for k in 1:size(Xₒ, 2)
            s = Symbol("x", k)
            push!(fixedSyms, s)
            push!(colSyms, s)
            push!(colData, Xₒ[:, k])  # Already Float64 from the matrix
        end

        randomSyms = Symbol[]
        for k in 1:size(Zₒ, 2)
            if !all(abs.(Zₒ[:, k] .- 1.0) .< 1e-12)
                s = Symbol("z", k)
                push!(randomSyms, s)
                push!(colSyms, s)
                push!(colData, Zₒ[:, k])
            end  # Already Float64 from the matrix
        end

        tableData = NamedTuple{Tuple(colSyms)}(Tuple(colData))

        if intercept
            fixedPart = isempty(fixedSyms) ? "1" : "1 + " * join(string.(fixedSyms), " + ")
        else
            fixedPart = isempty(fixedSyms) ? "0" : "0 + " * join(string.(fixedSyms), " + ")
        end

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

        β̂ = fixef(model)
        βcov = Matrix{Float64}(StatsAPI.vcov(model))

        # Sample fixed effects with regularization if needed
        β̇ = try
            β̂ + cholesky(Hermitian(βcov)).U' * randn(length(β̂))
        catch e1
            try
                βcov_reg = βcov + ridge * I(length(β̂))
                β̂ + cholesky(Hermitian(βcov_reg)).U' * randn(length(β̂))
            catch e2
                push!(loggedEvents, "2l.bin: Fixed effects sampling failed for $(yVar), using point estimates")
                β̂
            end
        end

        rancoef = Matrix{Float64}(ranef(model)[1])'
        q = size(rancoef, 2)
        ψ̂ = q == 1 ? reshape(var(vec(rancoef)), 1, 1) : cov(rancoef, dims = 1)

        # Ensure ψ̂ is positive definite
        ψ̂ = Hermitian(ψ̂ + ridge * I(q))

        λ = rancoef' * rancoef
        s = q * ψ̂
        ev = eigen(Hermitian(λ + s))
        eigenvalues = max.(ev.values, ridge)  # Ensure positive eigenvalues
        deco = ev.vectors * Diagonal(sqrt.(eigenvalues))

        ν = size(rancoef, 1) + q
        ψ̇InvScale = randWishart(ν, Diagonal(ones(q)))

        # Use pseudo-inverse with ridge for numerical stability
        ψ̇ = try
            inv(Hermitian(deco * ψ̇InvScale * transpose(deco) + ridge * I(q)))
        catch
            # Fallback to regularized estimate
            Hermitian(ψ̂ + ridge * I(q))
        end

        # Impute missing values using final coefficients
        gfₘ = gfFull[whereY]
        imps = Vector{Float64}(undef, whereCount)
        randEffectDraws = Dict{Int, Vector{Float64}}()

        # Prepare design matrices for prediction
        # If intercept was included in model, add column of 1s
        Xₘ_pred = intercept ? hcat(ones(size(Xₘ, 1)), Xₘ) : Xₘ
        Zₘ_pred = hcat(ones(size(Zₘ, 1)), Zₘ)  # Always add random intercept column

        for i in 1:whereCount
            classIdx = gfₘ[i]
            if !haskey(randEffectDraws, classIdx)
                randEffectDraws[classIdx] = vec(rand(Distributions.MvNormal(zeros(q), ψ̇)))
            end
            η = dot(Xₘ_pred[i, :], β̇) + dot(Zₘ_pred[i, :], randEffectDraws[classIdx])
            p = 1.0 / (1.0 + exp(-η))
            imps[i] = rand(Binomial(1, p))
        end

        # Convert back to original data type
        if eltype(yₒRaw) <: CategoricalValue
            levels_list = levels(yₒRaw[1])
            return CategoricalArray([levels_list[Int(val) + 1] for val in imps], levels_list)
        else
            return imps
        end
    end

    const TWO_LEVEL_BIN_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept::Bool = true, ridge::Float64 = 1e-4, kwargs...) -> begin
        twoLevelBinImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; intercept = intercept, ridge = ridge, kwargs...)
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

        mice(ct, m = 1, iter = 1, methods = methods, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)
    end
end
