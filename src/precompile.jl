@compile_workload begin
    using CSV

    data = CSV.read("test/data/cirrhosis.csv", DataFrame)

    imputedData = mice(data, progressReports = false)

    with(imputedData, data -> mean(data.Bilirubin))
end

precompile(plot, (Mids, String))
precompile(plot, (Mids, Int))