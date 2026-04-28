using BetaML, CategoricalArrays, CSV, DataFrames, GLM, Mice, MixedModels, Tables, Test, TypedTables

const CIRRHOSIS = joinpath(@__DIR__, "data", "cirrhosis.csv")

function make_two_level_workload()
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    z = randn(n)
    y = Vector{Union{Missing, Float64}}(@. 0.3 + 0.9 * x + 0.1 * cluster + 0.2 * randn())

    yMissingWithinClusters = copy(y)
    yMissingWithinClusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    yMissingEntireClusters = copy(y)
    yMissingEntireClusters[cluster .<= 3] .= missing

    return (
        within = (; y = yMissingWithinClusters, x = x, z = z, cluster = cluster),
        whole = (; y = yMissingEntireClusters, x = x, z = z, cluster = cluster),
    )
end

function make_two_level_binary_workload()
    n = 48
    cluster = repeat(1:8, inner = 6)
    x = randn(n)
    y = Vector{Union{Missing, Int}}(repeat([0, 1], inner = 1, outer = 24))

    yMissingWithinClusters = copy(y)
    yMissingWithinClusters[[2, 8, 13, 19, 25, 31, 37, 43]] .= missing

    return (; y = yMissingWithinClusters, x = x, cluster = cluster)
end

function missing_count(table)
    columnTable = Tables.columntable(table)
    return sum(sum(ismissing.(Tables.getcolumn(columnTable, name))) for name in propertynames(columnTable))
end

function assert_complete_tablelike(imputedData)
    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(missing_count(ct) for ct in imputedDataList) == 0
end

function assert_standard_workflow(data; methods = nothing, predictorMatrix = nothing, second_pass = false)
    imputedData = methods === nothing ?
        (predictorMatrix === nothing ?
            mice(data; progressReports = false) :
            mice(data; predictorMatrix = predictorMatrix, progressReports = false)) :
        (predictorMatrix === nothing ?
            mice(data; methods = methods, progressReports = false) :
            mice(data; methods = methods, predictorMatrix = predictorMatrix, progressReports = false))

    assert_complete_tablelike(imputedData)

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))
    @test length(analyses.analyses) == 5

    results = pool(analyses)
    @test length(results.coefs) == 7

    if second_pass
        imputedData = methods === nothing ?
            (predictorMatrix === nothing ?
                mice(imputedData; progressReports = false) :
                mice(imputedData; predictorMatrix = predictorMatrix, progressReports = false)) :
            (predictorMatrix === nothing ?
                mice(imputedData; methods = methods, progressReports = false) :
                mice(imputedData; methods = methods, predictorMatrix = predictorMatrix, progressReports = false))
        assert_complete_tablelike(imputedData)

        analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))
        @test length(analyses.analyses) == 5

        results = pool(analyses)
        @test length(results.coefs) == 7
    end
end

function run_two_level_method(data, method; nIterations = 0)
    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, :] .= 0
    predictorMatrix["y", "x"] = 2
    predictorMatrix["y", "cluster"] = -2

    if hasproperty(data, :z)
        predictorMatrix["y", "z"] = 1
    end

    methods = makeMethods(data)
    methods .= ""
    methods["y"] = method

    return mice(data, m = 1, iter = 1, methods = methods, predictorMatrix = predictorMatrix, progressReports = false, nIterations = nIterations, ridge = 1e-3)
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
    methods = makeMethods(data)
    methods .= method
    return methods
end

function cirrhosis_predictor_matrix(data)
    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false
    return predictorMatrix
end

@testset "Mice (PMM, DF)" begin
    data = cirrhosis_dataframe()
    predictorMatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; predictorMatrix = predictorMatrix, second_pass = true)
end

@testset "Mice (PMM, TT)" begin
    data = cirrhosis_table()
    predictorMatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; predictorMatrix = predictorMatrix)
end

@testset "Mice (norm, DF)" begin
    data = cirrhosis_norm_data()
    methods = cirrhosis_method_config(data, "sample")
    methods[["Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]] .= "norm"
    predictorMatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; methods = methods, predictorMatrix = predictorMatrix)
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
    predictorMatrix = cirrhosis_predictor_matrix(data)
    assert_standard_workflow(data; methods = methods, predictorMatrix = predictorMatrix)
end

@testset "Mice (2l, DF)" begin
    workload = make_two_level_workload()

    for method in ("2l.norm", "2l.pmm", "2lonly.mean", "2lonly.mode")
        imputedData = run_two_level_method(workload.within, method)
        assert_complete_tablelike(imputedData)
    end

    for method in ("2lonly.norm", "2lonly.pmm")
        imputedData = run_two_level_method(workload.whole, method)
        assert_complete_tablelike(imputedData)
    end
end

@testset "Mice (2l, TT)" begin
    workload = make_two_level_workload()
    dataWithin = Table(workload.within)
    dataWhole = Table(workload.whole)

    for method in ("2l.norm", "2l.pmm", "2lonly.mean", "2lonly.mode")
        imputedData = run_two_level_method(dataWithin, method)
        assert_complete_tablelike(imputedData)
    end

    for method in ("2lonly.norm", "2lonly.pmm")
        imputedData = run_two_level_method(dataWhole, method)
        assert_complete_tablelike(imputedData)
    end
end

@testset "Mice (2l.bin, DF)" begin
    binaryWorkload = make_two_level_binary_workload()
    imputedData = run_two_level_method(binaryWorkload, "2l.bin")
    print(imputedData.loggedEvents)
    @test length(imputedData.loggedEvents) == 0
end

@testset "Mice (2l.bin, TT)" begin
    binaryWorkload = Table(make_two_level_binary_workload())
    imputedData = run_two_level_method(binaryWorkload, "2l.bin")
    @test length(imputedData.loggedEvents) == 0
end