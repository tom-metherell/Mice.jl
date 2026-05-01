module MiceBetaMLExt
    using BetaML: fit!, RandomForestImputer, NONE
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Mice: Imputer, makemethods, mice, registerimputer!
    using PrecompileTools: @compile_workload
    using Random: rand, randperm

    function _rf_imputer_with_supported_kwargs(; n_trees::Int = 10, verbosity = NONE, kwargs...)
        imputerKw = (n_trees = n_trees, verbosity = verbosity)

        for (k, v) in pairs(kwargs)
            candidateKw = merge(imputerKw, NamedTuple{(k,)}((v,)))
            try
                RandomForestImputer(; candidateKw...)
                imputerKw = candidateKw
            catch e
                if !(e isa MethodError)
                    rethrow(e)
                end
            end
        end

        return RandomForestImputer(; imputerKw...)
    end

    function rf_impute(
        y::AbstractArray,
        X::Matrix{Float64},
        where_y::Vector{Bool};
        n_trees::Int = 10,
        verbosity = NONE,
        kwargs...
        )

        yDecat = y isa CategoricalArray || eltype(y) <: CategoricalValue ? Vector{String}(string.(y)) : y

        yX = Matrix{Union{Missing, eltype(yDecat), Float64}}(hcat(yDecat, X))

        yX[where_y, 1] .= missing

        rfImputer = _rf_imputer_with_supported_kwargs(; n_trees = n_trees, verbosity = verbosity, kwargs...)
        ŷX = fit!(rfImputer, yX)

        if y == yDecat
            return eltype(y) <: Integer ? round.(ŷX[where_y, 1], digits = 0) : ŷX[where_y, 1]
        end

        levelLookup = Dict(string.(levels(y)) .=> levels(y))
        return convert.(nonmissingtype(eltype(y)), getindex.(Ref(levelLookup), ŷX[where_y, 1]))
    end

    const RF_IMPUTER = Imputer((ydata, X, where_y, wherecount, yvar, itercounter, j, loggedevents; n_trees::Int = 10, verbosity = NONE, kwargs...) -> begin
        rf_impute(ydata, X, where_y; n_trees = n_trees, verbosity = verbosity, kwargs...)
    end)

    function __init__()
        registerimputer!("rf", RF_IMPUTER)
    end

    @compile_workload begin
        catPool = CategoricalPool(["a", "b", "c"])
        ct = (
            a = Vector{Union{Missing, Int}}(randperm(20)),
            b = Vector{Union{Missing, Float64}}(randperm(20)),
            c = Vector{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
            d = Vector{Union{Missing, Bool}}(rand(Bool, 20)),
            e = CategoricalArray{Union{Missing, Int}}(rand([1, 2, 3], 20)),
            f = CategoricalArray{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
            g = Vector{Union{Missing, CategoricalValue}}(rand([CategoricalValue(catPool, 1), CategoricalValue(catPool, 2), CategoricalValue(catPool, 3)], 20))
        )

        for col in ct
            col[rand(1:20, 1)] .= missing
        end

        rf_methods = makemethods(ct)
        rf_methods["b"] = "rf"
        imputedDataRf = mice(ct, m = 1, iter = 1, methods = rf_methods, progressreports = false)
    end

    export rf_impute
end