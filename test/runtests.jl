using CSV, DataFrames, GLM, Mice, Test

@testset "Mice" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame)
    colsWithMissings = ["Drug", "Ascites", "Hepatomegaly", "Spiders", "Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin", "Stage"]
    data[!, colsWithMissings] = allowmissing(data[!, colsWithMissings])
    for i in colsWithMissings
        replace!(data[!, i], "NA" => missing)
    end
    for i in ["Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]
        data[!, i] = passmissing(x -> parse(Float64, x)).(data[!, i])
    end

    theMethods = makeMethods(data)
    theMethods[["ID", "N_Days"]] .= ""

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, m = 20, iter = 15, methods = theMethods, predictorMatrix = predictorMatrix, threads = false, gcSchedule = 0.3)

    @test length(imputedData.loggedEvents) == 120

    imputedDataLong = complete(imputedData, "long")

    @test sum(ismissing.(Matrix(imputedDataLong))) == 0

    analyses = with(imputedData, data -> lm(@formula(N_Days ~ Drug + Age + Stage + Bilirubin), data))

    @test length(analyses.analyses) == 20

    results = pool(analyses)

    @test length(results.coefs) == 7
end