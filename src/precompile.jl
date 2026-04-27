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

    meanMethods = makeMethods(ct)
    meanMethods["b"] = "mean"
    imputedDataMean = mice(ct, m = 1, iter = 1, methods = meanMethods, progressReports = false)

    normMethods = meanMethods
    normMethods["b"] = "norm"
    imputedDataNorm = mice(ct, m = 1, iter = 1, methods = normMethods, progressReports = false)

    sampleMethods = normMethods
    sampleMethods[:] .= "sample"
    imputedDataSample = mice(ct, m = 1, iter = 1, methods = sampleMethods, progressReports = false)

    bindImputations(imputedDataPmm, imputedDataPmm)

    complete(imputedDataPmm, 1)
    listComplete(imputedDataPmm)
end

@compile_workload begin
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    yNorm = Vector{Union{Missing, Float64}}(@. 0.3 + 0.9 * x + 0.1 * cluster + 0.2 * randn())

    yNorm[[3, 9, 14, 21, 27, 34, 40, 46]] .= missing

    ct2l = (
        y = yNorm,
        x = x,
        cluster = cluster
    )

    predictorMatrix = makePredictorMatrix(ct2l)
    predictorMatrix[:, :] .= 0
    predictorMatrix["y", "x"] = 2
    predictorMatrix["y", "cluster"] = -2

    methods2lNorm = makeMethods(ct2l)
    methods2lNorm .= ""
    methods2lNorm["y"] = "2l.norm"
    mice(ct2l, m = 1, iter = 1, methods = methods2lNorm, predictorMatrix = predictorMatrix, progressReports = false, nIterations = 0, ridge = 1e-3)

    methods2lPmm = makeMethods(ct2l)
    methods2lPmm .= ""
    methods2lPmm["y"] = "2l.pmm"
    mice(ct2l, m = 1, iter = 1, methods = methods2lPmm, predictorMatrix = predictorMatrix, progressReports = false, nIterations = 0, ridge = 1e-3)

    methods2lOnlyMean = makeMethods(ct2l)
    methods2lOnlyMean .= ""
    methods2lOnlyMean["y"] = "2lonly.mean"
    mice(ct2l, m = 1, iter = 1, methods = methods2lOnlyMean, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)

    methods2lOnlyNorm = makeMethods(ct2l)
    methods2lOnlyNorm .= ""
    methods2lOnlyNorm["y"] = "2lonly.norm"
    mice(ct2l, m = 1, iter = 1, methods = methods2lOnlyNorm, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)

    methods2lOnlyPmm = makeMethods(ct2l)
    methods2lOnlyPmm .= ""
    methods2lOnlyPmm["y"] = "2lonly.pmm"
    mice(ct2l, m = 1, iter = 1, methods = methods2lOnlyPmm, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)
end
