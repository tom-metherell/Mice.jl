@compile_workload begin
    catpool = CategoricalPool(["a", "b", "c"])
    ct = (
        a = Vector{Union{Missing, Int}}(randperm(20)),
        b = Vector{Union{Missing, Float64}}(randperm(20)),
        c = Vector{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
        d = Vector{Union{Missing, Bool}}(rand(Bool, 20)),
        e = CategoricalArray{Union{Missing, Int}}(rand([1, 2, 3], 20)),
        f = CategoricalArray{Union{Missing, String}}(rand(["a", "b", "c"], 20)),
        g = Vector{Union{Missing, CategoricalValue}}(rand([CategoricalValue(catpool, 1), CategoricalValue(catpool, 2), CategoricalValue(catpool, 3)], 20))
    )

    for col in ct
        col[rand(1:20, 1)] .= missing
    end

    imputeddatapmm = mice(ct, m = 1, iter = 1, progressreports = false)

    meanmethods = makemethods(ct)
    meanmethods["b"] = "mean"
    imputeddatamean = mice(ct, m = 1, iter = 1, methods = meanmethods, progressreports = false)

    normmethods = meanmethods
    normmethods["b"] = "norm"
    imputeddatanorm = mice(ct, m = 1, iter = 1, methods = normmethods, progressreports = false)

    imputeddatanorm = normmethods
    imputeddatanorm[:] .= "sample"
    imputeddatasample = mice(ct, m = 1, iter = 1, methods = imputeddatanorm, progressreports = false)

    bindimputations(imputeddatapmm, imputeddatapmm)

    complete(imputeddatapmm, 1)
    listcomplete(imputeddatapmm)
end

@compile_workload begin
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    yNorm = Vector{Union{Missing, Float64}}(@. 0.3 + 0.9 * x + 0.1 * cluster + 0.2 * randn())

    ynormmissingwithinclusters = deepcopy(yNorm)
    ynormmissingwithinclusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    ynormmissingentireclusters = deepcopy(yNorm)
    ynormmissingentireclusters[cluster .== 3] .= missing

    ct2lmissingwithinclusters = (
        y = ynormmissingwithinclusters,
        x = x,
        cluster = cluster
    )

    ct2lmissingentireclusters = (
        y = ynormmissingentireclusters,
        x = x,
        cluster = cluster
    )

    predictormatrix = makepredictormatrix(ct2lmissingwithinclusters)
    predictormatrix[:, :] .= 0
    predictormatrix["y", "x"] = 2
    predictormatrix["y", "cluster"] = -2

    methods2lnorm = makemethods(ct2lmissingwithinclusters)
    methods2lnorm .= ""
    methods2lnorm["y"] = "2l.norm"
    mice(ct2lmissingwithinclusters, m = 1, iter = 1, methods = methods2lnorm, predictormatrix = predictormatrix, progressreports = false, n_iterations = 0, ridge = 1e-3)

    methods2lpmm = makemethods(ct2lmissingwithinclusters)
    methods2lpmm .= ""
    methods2lpmm["y"] = "2l.pmm"
    mice(ct2lmissingwithinclusters, m = 1, iter = 1, methods = methods2lpmm, predictormatrix = predictormatrix, progressreports = false, n_iterations = 0, ridge = 1e-3)

    methods2lonlymean = makemethods(ct2lmissingwithinclusters)
    methods2lonlymean .= ""
    methods2lonlymean["y"] = "2lonly.mean"
    mice(ct2lmissingwithinclusters, m = 1, iter = 1, methods = methods2lonlymean, predictormatrix = predictormatrix, progressreports = false, ridge = 1e-3)

    methods2lonlymode = makemethods(ct2lmissingwithinclusters)
    methods2lonlymode .= ""
    methods2lonlymode["y"] = "2lonly.mode"
    mice(ct2lmissingwithinclusters, m = 1, iter = 1, methods = methods2lonlymode, predictormatrix = predictormatrix, progressreports = false, ridge = 1e-3)

    methods2lonlynorm = makemethods(ct2lmissingentireclusters)
    methods2lonlynorm .= ""
    methods2lonlynorm["y"] = "2lonly.norm"
    mice(ct2lmissingentireclusters, m = 1, iter = 1, methods = methods2lonlynorm, predictormatrix = predictormatrix, progressreports = false, ridge = 1e-3)

    methods2lonlypmm = makemethods(ct2lmissingentireclusters)
    methods2lonlypmm .= ""
    methods2lonlypmm["y"] = "2lonly.pmm"
    mice(ct2lmissingentireclusters, m = 1, iter = 1, methods = methods2lonlypmm, predictormatrix = predictormatrix, progressreports = false, ridge = 1e-3)
end