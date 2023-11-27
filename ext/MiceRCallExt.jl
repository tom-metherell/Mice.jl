module MiceRCallExt
    using CategoricalArrays: CategoricalValue
    using Mice
    using NamedArrays: NamedArray
    using RCall: protect, unprotect, RClass, setclass!
    import RCall: sexp, sexpclass

    function sexp(::Type{RClass{:list}}, mids::Mids)
        r = protect(sexp(Dict(
            "data" => mids.data,
            "imp" => [NamedArray(
                eltype(mids.imputations[i]) <: CategoricalValue ? get.(mids.imputations[i]) : mids.imputations[i],
                (string.(findall(mids.imputeWhere[mids.visitSequence[i]])), string.(1:mids.m)),
                ("Rows", "Cols")
            ) for i in eachindex(mids.imputations)],
            "m" => mids.m,
            "where" => mids.imputeWhere,
            "blocks" => nothing,
            "call" => nothing,
            "nmis" => nothing,
            "method" => mids.methods,
            "predictorMatrix" => mids.predictorMatrix,
            "visitSequence" => mids.visitSequence,
            "formulas" => nothing,
            "post" => nothing,
            "blots" => nothing,
            "ignore" => nothing,
            "seed" => nothing,
            "iteration" => mids.iter,
            "lastSeedValue" => nothing,
            "chainMean" => mids.meanTraces,
            "chainVar" => mids.varTraces,
            "loggedEvents" => mids.loggedEvents,
            "version" => nothing,
            "date" => nothing
        )))
        setclass!(r, sexp("mids"))
        unprotect(1)
        return r
    end

    sexpclass(mids::Mids) = RClass{:list}
end