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

    yNormMissingWithinClusters = deepcopy(yNorm)
    yNormMissingWithinClusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    yNormMissingEntireClusters = deepcopy(yNorm)
    yNormMissingEntireClusters[cluster .== 3] .= missing

    ct2lMissingWithinClusters = (
        y = yNormMissingWithinClusters,
        x = x,
        cluster = cluster
    )

    ct2lMissingEntireClusters = (
        y = yNormMissingEntireClusters,
        x = x,
        cluster = cluster
    )

    predictorMatrix = makePredictorMatrix(ct2lMissingWithinClusters)
    predictorMatrix[:, :] .= 0
    predictorMatrix["y", "x"] = 2
    predictorMatrix["y", "cluster"] = -2

    methods2lNorm = makeMethods(ct2lMissingWithinClusters)
    methods2lNorm .= ""
    methods2lNorm["y"] = "2l.norm"
    mice(ct2lMissingWithinClusters, m = 1, iter = 1, methods = methods2lNorm, predictorMatrix = predictorMatrix, progressReports = false, nIterations = 0, ridge = 1e-3)

    methods2lPmm = makeMethods(ct2lMissingWithinClusters)
    methods2lPmm .= ""
    methods2lPmm["y"] = "2l.pmm"
    mice(ct2lMissingWithinClusters, m = 1, iter = 1, methods = methods2lPmm, predictorMatrix = predictorMatrix, progressReports = false, nIterations = 0, ridge = 1e-3)

    methods2lOnlyMean = makeMethods(ct2lMissingWithinClusters)
    methods2lOnlyMean .= ""
    methods2lOnlyMean["y"] = "2lonly.mean"
    mice(ct2lMissingWithinClusters, m = 1, iter = 1, methods = methods2lOnlyMean, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)

    methods2lOnlyMode = makeMethods(ct2lMissingWithinClusters)
    methods2lOnlyMode .= ""
    methods2lOnlyMode["y"] = "2lonly.mode"
    mice(ct2lMissingWithinClusters, m = 1, iter = 1, methods = methods2lOnlyMode, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)

    methods2lOnlyNorm = makeMethods(ct2lMissingEntireClusters)
    methods2lOnlyNorm .= ""
    methods2lOnlyNorm["y"] = "2lonly.norm"
    mice(ct2lMissingEntireClusters, m = 1, iter = 1, methods = methods2lOnlyNorm, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)

    methods2lOnlyPmm = makeMethods(ct2lMissingEntireClusters)
    methods2lOnlyPmm .= ""
    methods2lOnlyPmm["y"] = "2lonly.pmm"
    mice(ct2lMissingEntireClusters, m = 1, iter = 1, methods = methods2lOnlyPmm, predictorMatrix = predictorMatrix, progressReports = false, ridge = 1e-3)
end