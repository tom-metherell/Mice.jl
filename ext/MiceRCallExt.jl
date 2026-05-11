module MiceRCallExt
    using AxisArrays: axes, AxisArray
    using CategoricalArrays: CategoricalValue
    using Mice
    using RCall: protect, RClass, @R_str, setclass!, unprotect, VecSxp
    import RCall: rcopy, rcopytype, sexp, sexpclass
    using Tables: columnnames

    function rcopy(::Type{Mids}, s::Ptr{VecSxp})
    protect(s)
        try
            nt = rcopy(NamedTuple, s)
            Mids(
                nt.data,
                [Matrix(nt.imp[Symbol(i)]) for i in nt.visitSequence],
                AxisArray(
                    [Vector{Bool}(nt.where[:, i]) for i in eachindex(1:size(nt.data, 2))], 
                    names(nt.data)
                ),
                Int(nt.m),
                AxisArray(nt.method, names(nt.data)),
                AxisArray(Matrix{Int}(nt.predictorMatrix), names(nt.data), names(nt.data)),
                nt.visitSequence,
                Int(nt.iteration),
                [nt.chainMean[findfirst(names(nt.data) .== i), :, :] for i in nt.visitSequence],
                [nt.chainVar[findfirst(names(nt.data) .== i), :, :] for i in nt.visitSequence],
                ["This Mids object originated in R. Logged events have not been transferred."]
            )
        finally
            unprotect(1)
        end
    end

    rcopytype(::Type{RClass{:mids}}, x::Ptr{VecSxp}) = Mids

    function sexp(::Type{RClass{:mids}}, mids::Mids)
        r = protect(sexp(Dict(
            "data" => mids.data,
            "imp" => AxisArray(
                [
                    try 
                        AxisArray(
                            eltype(mids.imputations[findfirst(mids.visitsequence .== i)]) <: CategoricalValue ? get.(mids.imputations[findfirst(mids.visitsequence .== i)]) : mids.imputations[findfirst(mids.visitsequence .== i)],
                            string.(findall(mids.imputewhere[mids.visitsequence[findfirst(mids.visitsequence .== i)]])),
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
            "where" => AxisArray(reduce(hcat, mids.imputewhere), Base.axes(reduce(hcat, mids.imputewhere), 1), axes(mids.imputewhere)[1][:]),
            "blocks" => nothing,
            "call" => R"match.call()",
            "nmis" => nothing,
            "method" => AxisArray(mids.methods, collect(string.(columnnames(mids.data)))),
            "predictorMatrix" => AxisArray(Matrix{Int}(mids.predictormatrix), collect(string.(columnnames(mids.data))), collect(string.(columnnames(mids.data)))),
            "visitSequence" => mids.visitsequence,
            "formulas" => nothing,
            "post" => nothing,
            "blots" => nothing,
            "ignore" => nothing,
            "seed" => nothing,
            "iteration" => mids.iter,
            "lastSeedValue" => nothing,
            "chainMean" => AxisArray(
                cat(dims = 3, [reduce(hcat, [[mids.meantraces[findfirst(mids.visitsequence .== i)][j, k] for i in mids.visitsequence] for j in 1:mids.iter]) for k in 1:mids.m]...),
                mids.visitsequence,
                string.(1:mids.iter),
                ["Chain $k" for k in 1:mids.m]
            ),
            "chainVar" => AxisArray(
                cat(dims = 3, [reduce(hcat, [[mids.vartraces[findfirst(mids.visitsequence .== i)][j, k] for i in mids.visitsequence] for j in 1:mids.iter]) for k in 1:mids.m]...),
                mids.visitsequence,
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

    sexpclass(mids::Mids) = RClass{:mids}

    export rcopy, rcopytype, sexp, sexpclass
end