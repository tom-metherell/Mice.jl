module MicePlotsExt
    using Mice: makemethods, mice
    using Plots
    using PrecompileTools: @compile_workload
    using Random: rand, randn

    @compile_workload begin
        n = 40
        x = randn(n)
        z = randn(n)
        y = Vector{Union{Missing, Float64}}(@. 0.5 + 1.1 * x - 0.4 * z + 0.3 * randn())
        y[[3, 9, 17, 25, 33]] .= missing

        ct = (
            y = y,
            x = x,
            z = z
        )

        methods = makemethods(ct)
        methods .= ""
        methods["y"] = "norm"

        mids = mice(ct, m = 1, iter = 1, methods = methods, progressreports = false)
        Plots.plot(mids, "y")
    end
end
