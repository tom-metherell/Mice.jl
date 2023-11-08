@compile_workload begin
    data = DataFrame(
        :col1 => Vector{Union{Missing, String}}(["a", "b", "c", missing, "e"]),
        :col2 => Vector{Union{Missing, AbstractFloat}}(vcat(randn(2), missing, randn(2))),
        :col3 => Vector{Union{Missing, Int}}(vcat(rand(1:5, 2), missing, rand(1:5, 2)))
    )

    mice(data)

    with(data, data -> mean(data.col2))
end

precompile(plot, (Mids, String))