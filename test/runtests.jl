using CategoricalArrays, CSV, DataFrames, GLM, Mice, Tables, Test, TypedTables

@testset "Mice (DF)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = false, gcSchedule = 0.0, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (DF, threaded)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = true, gcSchedule = 0.0, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (TT)" begin
    data = CSV.read("data/cirrhosis.csv", Table, missingstring = "NA")

    data = Table((; zip([i for i in Tables.columnnames(data) if i != :Stage], [Tables.getcolumn(data, i) for i in Tables.columnnames(data) if i != :Stage])...), Stage = categorical(data.Stage))

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = false, gcSchedule = 0.0, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum([sum([sum(ismissing.(ct[i])) for i in eachindex(ct)]) for ct in columntable.(imputedDataList)]) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end