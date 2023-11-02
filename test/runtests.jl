using CSV, DataFrames, Mice, Test

@testset "mice" begin
    data = CSV.read("data/cirrhosis.csv", DataFrame)
    colsWithMissings = ["Drug", "Ascites", "Hepatomegaly", "Spiders", "Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin", "Stage"]
    data[!, colsWithMissings] = allowmissing(data[!, colsWithMissings])
    for i in colsWithMissings
        replace!(data[!, i], "NA" => missing)
    end
    for i in ["Cholesterol", "Copper", "Alk_Phos", "SGOT", "Tryglicerides", "Platelets", "Prothrombin"]
        data[!, i] = passmissing(x -> parse(Float64, x)).(data[!, i])
    end

    methods = makeMethods(data)
    methods[["ID", "N_Days"]] .= ""

    predictorMatrix = makePredictorMatrix(data)
    predictorMatrix[:, ["ID", "N_Days"]] .= false

    imputedData = mice(data, m = 20, iter = 15)
end
