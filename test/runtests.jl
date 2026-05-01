using BetaML, CategoricalArrays, CSV, DataFrames, GLM, Mice, MixedModels, Tables, Test, TypedTables

const CIRRHOSIS = joinpath(@__DIR__, "data", "cirrhosis.csv")

function make_two_level_workload()
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    z = randn(n)
    y = Vector{Union{Missing, Float64}}(@. 0.3 + 0.9 * x + 0.1 * cluster + 0.2 * randn())

    ymissingwithinclusters = copy(y)
    ymissingwithinclusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    ymissingentireclusters = copy(y)
    ymissingentireclusters[cluster .<= 3] .= missing

    return (
        within = (; y = ymissingwithinclusters, x = x, z = z, cluster = cluster),
        whole = (; y = ymissingentireclusters, x = x, z = z, cluster = cluster),
    )
end

function make_two_level_binary_workload()
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    y = Vector{Union{Missing, Int}}(repeat([0, 1], outer = 24))

    ymissingwithinclusters = copy(y)
    ymissingwithinclusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    return (; y = ymissingwithinclusters, x = x, cluster = cluster)
end

function missing_count(table)
    ct = Tables.columntable(table)
    return sum(sum(ismissing.(Tables.getcolumn(ct, name))) for name in propertynames(ct))
end

function assert_complete_tablelike(imputeddata)
    @test length(imputeddata.loggedevents) == 0

    imputeddatalist = listcomplete(imputeddata)

    @test sum(missing_count(ct) for ct in imputeddatalist) == 0
end

function assert_standard_workflow(data; methods = nothing, predictormatrix = nothing, second_pass = false)
    imputeddata = methods === nothing ?
        (predictormatrix === nothing ?
            mice(data; progressreports = false) :
            mice(data; predictormatrix = predictormatrix, progressreports = false)) :
        (predictormatrix === nothing ?
            mice(data; methods = methods, progressreports = false) :
            mice(data; methods = methods, predictormatrix = predictormatrix, progressreports = false))

    assert_complete_tablelike(imputeddata)

    analyses = with(imputeddata, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))
    @test length(analyses.analyses) == 5

    results = pool(analyses)
    @test length(results.coefs) == 7

    if second_pass
        imputeddata = methods === nothing ?
            (predictormatrix === nothing ?
                mice(imputeddata; progressreports = false) :
                mice(imputeddata; predictormatrix = predictormatrix, progressreports = false)) :
            (predictormatrix === nothing ?
                mice(imputeddata; methods = methods, progressreports = false) :
                mice(imputeddata; methods = methods, predictormatrix = predictormatrix, progressreports = false))
        assert_complete_tablelike(imputeddata)

        analyses = with(imputeddata, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))
        @test length(analyses.analyses) == 5

        results = pool(analyses)
        @test length(results.coefs) == 7
    end
end

function run_two_level_method(data, method; n_iterations = 0, intercept = true)
    predictormatrix = makepredictormatrix(data)
    predictormatrix[:, :] .= 0
    predictormatrix["y", "x"] = 2
    predictormatrix["y", "cluster"] = -2

    if hasproperty(data, :z)
        predictormatrix["y", "z"] = 1
    end

    methods = makemethods(data)
    methods .= ""
    methods["y"] = method

    return mice(data, m = 1, iter = 1, methods = methods, predictormatrix = predictormatrix, progressreports = false, n_iterations = n_iterations, ridge = 1e-3, intercept = intercept)
end

function cirrhosis_dataframe()
    data = CSV.read(CIRRHOSIS, DataFrame, missingstring = "NA")
    data.Stage = categorical(data.Stage)
    return data
end

function cirrhosis_table()
    data = CSV.read(CIRRHOSIS, Table, missingstring = "NA")
    return Table((; zip([i for i in Tables.columnnames(data) if i ≠ :Stage], [Tables.getcolumn(data, i) for i in Tables.columnnames(data) if i ≠ :Stage])...), Stage = categorical(data.Stage))
end

function cirrhosis_norm_data()
    data = cirrhosis_dataframe()
    data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]] = convert.(Union{Float64, Missing}, data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]])
    return data
end

function cirrhosis_method_config(data, method)
    methods = makemethods(data)
    methods .= method
    return methods
end

function cirrhosis_predictor_matrix(data)
    predictormatrix = makepredictormatrix(data)
    predictormatrix[:, ["ID", "N_Days"]] .= false
    return predictormatrix
end

@testset "Mice (PMM, DF)" begin
    data = cirrhosis_dataframe()
    predictormatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; predictormatrix = predictormatrix, second_pass = true)
end

@testset "Mice (PMM, TT)" begin
    data = cirrhosis_table()
    predictormatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; predictormatrix = predictormatrix)
end

@testset "Mice (norm, DF)" begin
    data = cirrhosis_norm_data()
    methods = cirrhosis_method_config(data, "sample")
    methods[["Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]] .= "norm"
    predictormatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; methods = methods, predictormatrix = predictormatrix)
end

@testset "Mice (sample, DF)" begin
    data = cirrhosis_dataframe()
    methods = cirrhosis_method_config(data, "sample")
    assert_standard_workflow(data; methods = methods)
end

@testset "Mice (mean, DF)" begin
    data = cirrhosis_norm_data()
    methods = cirrhosis_method_config(data, "sample")
    methods[["Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]] .= "mean"
    assert_standard_workflow(data; methods = methods)
end

@testset "Mice (RF, DF)" begin
    data = cirrhosis_dataframe()
    methods = cirrhosis_method_config(data, "rf")
    predictormatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; methods = methods, predictormatrix = predictormatrix)
end

@testset "Mice (2l, DF)" begin
    workload = make_two_level_workload()

    for method in ("2l.norm", "2l.pmm", "2lonly.mean", "2lonly.mode")
        imputeddata = run_two_level_method(workload.within, method)
        assert_complete_tablelike(imputeddata)
    end

    for method in ("2lonly.norm", "2lonly.pmm")
        imputeddata = run_two_level_method(workload.whole, method)
        assert_complete_tablelike(imputeddata)
    end
end

@testset "Mice (2l, TT)" begin
    workload = make_two_level_workload()
    datawithin = Table(workload.within)
    datawhole = Table(workload.whole)

    for method in ("2l.norm", "2l.pmm", "2lonly.mean", "2lonly.mode")
        imputeddata = run_two_level_method(datawithin, method)
        assert_complete_tablelike(imputeddata)
    end

    for method in ("2lonly.norm", "2lonly.pmm")
        imputeddata = run_two_level_method(datawhole, method)
        assert_complete_tablelike(imputeddata)
    end
end

@testset "Mice (2l.bin, DF)" begin
    binaryworkload = make_two_level_binary_workload()
    imputeddata = run_two_level_method(binaryworkload, "2l.bin")
    print(imputeddata.loggedevents)
    @test length(imputeddata.loggedevents) == 0
end

@testset "Mice (2l.bin, TT, no intercept)" begin
    binaryworkload = Table(make_two_level_binary_workload())
    imputeddata = run_two_level_method(binaryworkload, "2l.bin", intercept = false)
    @test length(imputeddata.loggedevents) == 0
end