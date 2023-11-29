using CategoricalArrays, CSV, DataFrames, GLM, Mice, Tables, Test, TypedTables

@testset "Mice (PMM, DF)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = false, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (PMM, DF, threaded)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = true, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (PMM, TT)" begin
    data = CSV.read("data/cirrhosis.csv", Table, missingstring = "NA")

    data = Table((; zip([i for i in Tables.columnnames(data) if i != :Stage], [Tables.getcolumn(data, i) for i in Tables.columnnames(data) if i != :Stage])...), Stage = categorical(data.Stage))

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, predictorMatrix = predictorMatrix, threads = false, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum([sum([sum(ismissing.(ct[i])) for i in eachindex(ct)]) for ct in columntable.(imputedDataList)]) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (norm, DF)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]] = convert.(Union{Float64, Missing}, data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]])

    theMethods = makeMethods(data)
    theMethods[["Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]] .= "norm"

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, methods = theMethods, predictorMatrix = predictorMatrix, threads = false, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (sample, DF)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    theMethods = makeMethods(data)
    theMethods .= "sample"

    imputedData = mice(data, iter = 1, methods = theMethods, threads = false, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end

@testset "Mice (mean, DF)" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame, missingstring = "NA")

    data.Stage = categorical(data.Stage)

    data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]] = convert.(Union{Float64, Missing}, data[!, ["Age", "Cholesterol", "Copper", "Tryglicerides", "Platelets"]])

    theMethods = makeMethods(data)
    theMethods .= "sample"
    theMethods[["Age", "Bilirubin", "Cholesterol", "Albumin", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]] .= "mean"

    imputedData = mice(data, iter = 1, methods = theMethods, threads = false, progressReports = false)

    @test length(imputedData.loggedEvents) == 0

    imputedDataList = listComplete(imputedData)

    @test sum(sum.(ismissing.(Matrix.(imputedDataList)))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 5

    results = pool(analyses)

    @test length(results.coefs) == 7
end