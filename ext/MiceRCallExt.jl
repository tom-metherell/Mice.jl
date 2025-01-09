module MiceRCallExt
    using AxisArrays: axes, AxisArray
    using CategoricalArrays: CategoricalValue
    using Mice
    using RCall: protect, unprotect, RClass, setclass!
    import RCall: @R_str, sexp, sexpclass
    using Tables: columnnames

    function sexp(::Type{RClass{:list}}, mids::Mids)
        r = protect(sexp(Dict(
            "data" => mids.data,
            "imp" => AxisArray(
                [
                    try 
                        AxisArray(
                            eltype(mids.imputations[findfirst(mids.visitSequence .== i)]) <: CategoricalValue ? get.(mids.imputations[findfirst(mids.visitSequence .== i)]) : mids.imputations[findfirst(mids.visitSequence .== i)],
                            string.(findall(mids.imputeWhere[mids.visitSequence[findfirst(mids.visitSequence .== i)]])),
                            string.(1:mids.m)
                        )
                    catch
                        AxisArray(
                            Matrix{Any}(undef, 0, mids.m),
                            [],
                            string.(1:mids.m)
                        )
                    end for i in collect(string.(columnnames(mids.data)))
                ], 
                collect(string.(columnnames(mids.data)))
            ),
            "m" => mids.m,
            "where" => AxisArray(reduce(hcat, mids.imputeWhere), Base.axes(reduce(hcat, mids.imputeWhere), 1), axes(mids.imputeWhere)[1][:]),
            "blocks" => nothing,
            "call" => R"match.call()",
            "nmis" => nothing,
            "method" => AxisArray(mids.methods, collect(string.(columnnames(mids.data)))),
            "predictorMatrix" => AxisArray(Matrix{Int}(mids.predictorMatrix), collect(string.(columnnames(mids.data))), collect(string.(columnnames(mids.data)))),
            "visitSequence" => mids.visitSequence,
            "formulas" => nothing,
            "post" => nothing,
            "blots" => nothing,
            "ignore" => nothing,
            "seed" => nothing,
            "iteration" => mids.iter,
            "lastSeedValue" => nothing,
            "chainMean" => AxisArray(
                cat(dims = 3, [reduce(hcat, [[mids.meanTraces[findfirst(mids.visitSequence .== i)][j, k] for i in mids.visitSequence] for j in 1:mids.iter]) for k in 1:mids.m]...),
                mids.visitSequence,
                string.(1:mids.iter),
                ["Chain $k" for k in 1:mids.m]
            ),
            "chainVar" => AxisArray(
                cat(dims = 3, [reduce(hcat, [[mids.varTraces[findfirst(mids.visitSequence .== i)][j, k] for i in mids.visitSequence] for j in 1:mids.iter]) for k in 1:mids.m]...),
                mids.visitSequence,
                string.(1:mids.iter),
                ["Chain $k" for k in 1:mids.m]
            ),
            "loggedEvents" => nothing,
            "version" => nothing,
            "date" => nothing
        )))
        setclass!(r, sexp("mids"))
        unprotect(1)
        return r
    end

    sexpclass(mids::Mids) = RClass{:list}

    export sexp, sexpclass
end