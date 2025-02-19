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

    imputedDataPmm = mice(ct, m = 1, iter = 1, progressReports = false)

    meanLRMethods = makeMethods(ct)
    meanLRMethods["b"] = "mean"
    meanLRMethods["d"] = "logreg"
    imputedDataMeanLR = mice(ct, m = 1, iter = 1, methods = meanLRMethods, progressReports = false)

    normMethods = makeMethods(ct)
    normMethods["b"] = "norm"
    imputedDataNorm = mice(ct, m = 1, iter = 1, methods = normMethods, progressReports = false)

    sampleMethods = normMethods
    sampleMethods[:] .= "sample"
    imputedDataSample = mice(ct, m = 1, iter = 1, methods = sampleMethods, progressReports = false)

    bindImputations(imputedDataPmm, imputedDataPmm)

    complete(imputedDataPmm, 1)
    listComplete(imputedDataPmm)
end
