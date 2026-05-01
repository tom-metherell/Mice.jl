module MiceDataFramesExt
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue
    using DataFrames: DataFrame
    using Mice: bindimputations, complete, listcomplete, makemethods, mice
    using PrecompileTools: @compile_workload
    using Random: rand, randperm

    @compile_workload begin
        catPool = CategoricalPool(["a", "b", "c"])
        df = DataFrame(
            a = Vector{Union{Missing, Int}}(randperm(20)),
            b = Vector{Union{Missing, Float64}}(randperm(20)),
            c = Vector{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
            d = Vector{Union{Missing, Bool}}(rand(Bool, 20)),
            e = CategoricalArray{Union{Missing, Int}}(rand([1, 2, 3], 20)),
            f = CategoricalArray{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
            g = Vector{Union{Missing, CategoricalValue}}(rand([CategoricalValue(catPool, 1), CategoricalValue(catPool, 2), CategoricalValue(catPool, 3)], 20))
        )

        for col in axes(df, 2)
            df[rand(1:20, 1), col] .= missing
        end

        imputeddatapmm = mice(df, m = 1, iter = 1, progressreports = false)

        meanmethods = makemethods(df)
        meanmethods["b"] = "mean"
        imputeddatamean = mice(df, m = 1, iter = 1, methods = meanmethods, progressreports = false)

        normmethods = meanmethods
        normmethods["b"] = "norm"
        imputeddatanorm = mice(df, m = 1, iter = 1, methods = normmethods, progressreports = false)

        imputeddatanorm = normmethods
        imputeddatanorm[:] .= "sample"
        imputeddatasample = mice(df, m = 1, iter = 1, methods = imputeddatanorm, progressreports = false)

        bindimputations(imputeddatapmm, imputeddatapmm)

        complete(imputeddatapmm, 1)
        listcomplete(imputeddatapmm)
    end
end