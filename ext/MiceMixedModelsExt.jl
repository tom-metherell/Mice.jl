module MiceMixedModelsExt
    using CategoricalArrays: categorical, CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Distributions
    using LinearAlgebra: Diagonal, Symmetric, cholesky, dot, eigen, inv, transpose, Hermitian, I
    using Mice: Imputer, makemethods, makepredictormatrix, mice, preparetwolevelimputationinputs, randwishart, registerimputer!, symridge
    using MixedModels: GeneralizedLinearMixedModel, fit!, fixef, ranef
    using PrecompileTools: @compile_workload
    using Random: rand, randn
    using Statistics: cov, var
    using StatsBase: countmap
    import StatsAPI
    import StatsModels

    # The twolevelbin_impute! function includes a ! as it updates loggedevents in place
    function twolevelbin_impute!(
        ydata::AbstractArray,
        X::Matrix{Float64},
        where_y::Vector{Bool},
        wherecount::Int,
        types::Vector{Int},
        yvar::String,
        itercounter::Int,
        j::Int,
        loggedevents::Vector{String};
        intercept::Bool = true,
        ridge::Float64 = 1e-4,
        unusedkwargs...
        )

        if wherecount == 0
            return Vector{eltype(ydata)}(undef, 0)
        end

        if isempty(findall(types .== 2))
            throw(ArgumentError("Two-level binary requires at least one random predictor (coded 2)."))
        end

        prep = preparetwolevelimputationinputs(X, where_y, types; intercept = false)
        gf_full = prep.gf_full
        gf = prep.gf
        nclasses = prep.nclasses
        randomcols = prep.randomcols
        fixedcols = prep.fixedcols

        yₒ_raw = ydata[.!where_y]
        yₒ = Vector{Float64}([y isa CategoricalValue ? y.ref - 1 : convert(Float64, y) for y in yₒ_raw])

        if length(unique(gf)) < nclasses
            throw(ArgumentError("Two-level binary requires at least one observed outcome per class."))
        end

        # Check that y is binary
        uniquevals = unique(yₒ)
        if length(uniquevals) != 2 || !all(v ∈ [0.0, 1.0] for v in uniquevals)
            throw(ArgumentError("Two-level binary imputation requires a binary outcome variable. Got unique values: $uniquevals"))
        end

        yObs = Int.(round.(yₒ))
        obsMask = .!where_y

        Xₒ = Matrix{Float64}(X[obsMask, fixedcols])
        Zₒ = Matrix{Float64}(X[obsMask, randomcols])
        Xₘ = Matrix{Float64}(X[where_y, fixedcols])
        Zₘ = Matrix{Float64}(X[where_y, randomcols])

        if size(Zₒ, 2) == 0
            throw(ArgumentError("Two-level binary requires at least one random-effect column in the design matrix."))
        end

        fixedsyms = Symbol[]
        randomsyms = Symbol[]
        colsyms = Symbol[:y, :cluster]
        coldata = Any[yObs, categorical(gf)]

        labelmap = Dict{Int, Symbol}()
        addedcols = Set{Int}()

        allcols = unique(vcat(fixedcols, randomcols))

        for col in allcols
            if !haskey(labelmap, col)
                labelmap[col] = Symbol("x", col)
            end
            s = labelmap[col]
            
            # Add to colData only once
            if !(col in addedcols)
                push!(colsyms, s)
                push!(coldata, X[obsMask, col])
                push!(addedcols, col)
            end
        end
        
        # Build fixed and random symbol lists using the shared labelmap
        for k in eachindex(fixedcols)
            push!(fixedsyms, labelmap[fixedcols[k]])
        end
        
        for k in eachindex(randomcols)
            push!(randomsyms, labelmap[randomcols[k]])
        end

        tabledata = NamedTuple{Tuple(colsyms)}(Tuple(coldata))

        if intercept
            fixedpart = isempty(fixedsyms) ? "1" : "1 + " * join(string.(fixedsyms), " + ")
        else
            fixedpart = isempty(fixedsyms) ? "0" : "0 + " * join(string.(fixedsyms), " + ")
        end

        randompart = isempty(randomsyms) ? "(1 | cluster)" : "(1 + " * join(string.(randomsyms), " + ") * " | cluster)"

        formulaString = "y ~ " * fixedpart * " + " * randompart
        formulaExpr = Meta.parse(formulaString)
        f = Core.eval(@__MODULE__, :(StatsModels.@formula($formulaExpr)))

        model = nothing
        try
            model = fit!(GeneralizedLinearMixedModel(f, tabledata, Distributions.Bernoulli()))
        catch err
            push!(loggedevents, "2l.bin: MixedModels fit failed for $(yvar) (iter $(itercounter), chain $(j)): $(err)")
            return ydata[where_y]
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
                push!(loggedevents, "2l.bin: Fixed effects sampling failed for $(yvar), using point estimates")
                β̂
            end
        end

        rancoef = Matrix{Float64}(ranef(model)[1])'
        q = size(rancoef, 2)
        ψ̂ = q == 1 ? reshape(var(vec(rancoef)), 1, 1) : cov(rancoef, dims = 1)

        λ = rancoef' * rancoef
        s = q * Hermitian(ψ̂)
        ev = eigen(Hermitian(λ + s))
        eigenvalues = max.(ev.values, ridge)  # Ensure positive eigenvalues
        deco = ev.vectors * Diagonal(sqrt.(eigenvalues))

        ν = size(rancoef, 1) + q
        ψ̇InvScale = randwishart(ν, Diagonal(ones(q)))

        # Use pseudo-inverse with ridge for numerical stability
        ψ̇ = try
            inv(Hermitian(deco * ψ̇InvScale * transpose(deco)))
        catch
            # Fallback to regularized estimate
            Hermitian(ψ̂ + ridge * I(q))
        end

        # Impute missing values using final coefficients
        gfₘ = gf_full[where_y]
        imps = Vector{Float64}(undef, wherecount)
        randeffectdraws = Dict{Int, Vector{Float64}}()

        # Prepare design matrices for prediction
        # If intercept was included in model, add column of 1s
        Xₘ_pred = intercept ? hcat(ones(size(Xₘ, 1)), Xₘ) : Xₘ
        Zₘ_pred = hcat(ones(size(Zₘ, 1)), Zₘ)  # Always add random intercept column

        for i in 1:wherecount
            classindex = gfₘ[i]
            if !haskey(randeffectdraws, classindex)
                randeffectdraws[classindex] = vec(rand(Distributions.MvNormal(zeros(q), ψ̇)))
            end
            η = dot(Xₘ_pred[i, :], β̇) + dot(Zₘ_pred[i, :], randeffectdraws[classindex])
            p = 1.0 / (1.0 + exp(-η))
            imps[i] = rand(Binomial(1, p))
        end

        # Convert back to original data type
        if eltype(yₒ_raw) <: CategoricalValue
            levelslist = levels(yₒ_raw[1])
            return CategoricalArray([levelslist[Int(val) + 1] for val in imps], levelslist)
        else
            return imps
        end
    end

    const TWO_LEVEL_BIN_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; intercept::Bool = true, ridge::Float64 = 1e-4, kwargs...) -> begin
        twolevelbin_impute!(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; intercept = intercept, ridge = ridge, kwargs...)
    end; twolevel = true)

    function __init__()
        registerimputer!("2l.bin", TWO_LEVEL_BIN_IMPUTER)
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

        methods = makemethods(ct)
        methods .= ""
        methods["y"] = "2l.bin"

        predictormatrix = makepredictormatrix(ct)
        predictormatrix[:, :] .= 0
        predictormatrix["y", "x"] = 2
        predictormatrix["y", "cluster"] = -2

        mice(ct, m = 1, iter = 1, methods = methods, predictormatrix = predictormatrix, progressreports = false, ridge = 1e-3)
    end
end
