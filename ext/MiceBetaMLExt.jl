module MiceBetaMLExt
    using BetaML: fit!, RandomForestImputer, NONE
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Mice: makeMethods, mice
    import Mice: rfImpute!
    using PrecompileTools: @compile_workload
    using Random: rand, randperm

    function rfImpute!(
        y::AbstractArray,
        X::Matrix{Float64},
        whereY::Vector{Bool};
        n_trees::Int = 10,
        verbosity = NONE,
        kwargs...
        )

        yDecat = y isa CategoricalArray || eltype(y) <: CategoricalValue ? Vector{String}(string.(y)) : y

        yX = Matrix{Union{Missing, eltype(yDecat), Float64}}(hcat(yDecat, X))

        yX[whereY, 1] .= missing

        ŷX = fit!(RandomForestImputer(; n_trees = n_trees, verbosity = verbosity, kwargs...), yX)

        return y == yDecat ? (eltype(y) <: Integer ? round.(ŷX[whereY, 1], digits = 0) : ŷX[whereY, 1]) : parse.(eltype(levels(y)), ŷX[whereY, 1])
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

        rfMethods = makeMethods(ct)
        rfMethods["b"] = "rf"
        imputedDataRf = mice(ct, m = 1, iter = 1, methods = rfMethods, progressReports = false)
    end

    export rfImpute!
end