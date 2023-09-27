"""
    complete(
        mids::Mids,
        imputation::Union{Int, Nothing} = nothing;
        action::Union{String, Nothing} = nothing
    )

Produces (a) DataFrame(s) with missings replaced with multiply imputed values.

The multiply imputed dataset (`Mids`) object must be supplied first.

The `imputation` argument is optional.
If specified, the function will return a single DataFrame, based on that
imputation alone.

The `action` argument is optional.
If specified as "long", the function will return a single DataFrame, containing
the results of each imputation in succession with an identifier (`imp`).
If specified as "list", the function will return a vector of individual
DataFrames.

Exactly one of `imputation` and `action` must have a non-`nothing` value.
"""
function complete(
    mids::Mids,
    imputation::Union{Int, Nothing} = nothing;
    action::Union{String, Nothing} = nothing
)

    if imputation === nothing && action === nothing
        throw(ArgumentError("Too few arguments."))
    end

    if imputation != nothing
        data = mids.data

        for i in eachindex(mids.visitSequence)
            var = mids.visitSequence[i]
            data[ismissing.(data[:, var]), var] = data.imputations[i][:, imputation]
        end

        return data
    else
        if !(action ∈ ["list", "long"])
            throw(ArgumentError("Action not defined. Valid arguments: \"list\", \"long\""))
        else
            data = Vector{DataFrame}(undef, mids.m)
            for i in 1:mids.m
                data[i] = mids.data
                for j in eachindex(mids.visitSequence)
                    var = mids.visitSequence[i]
                    data[i][ismissing.(data[:, var]), var] = data.imputations[j][:, i]
                end
            end
            if action == "list"
                return data
            else
                return vcat(data...)
            end
        end
    end
end

"""
    struct Mira(analyses::Vector)

A multiply imputed repeated analyses object.

The analyses are stored as a vector of analyses of individual imputations.
"""
struct Mira(analyses::Vector)

"""
    function with(
        mids::Mids,
        functionCall::String
    )

Conducts repeated analyses of a multiply imputed dataset (`Mids`).

The function takes two arguments: firstly the `Mids` object itself,
then a function call (`functionCall`) expressed as a string.

The function call should reflect that used to analyse an ordinary DataFrame.
For the data argument in the function substitute the term `data`.
For example: `with(mids, "lm(@formula(y ~ x1 + x2), data, Bernoulli())")`
"""
function with(
    mids::Mids, 
    functionCall::String
)
    analyses = Vector{Any}(undef, mids.m)

    for i in 1:mids.m
        data = complete(mids, i)
        analyses[i] = eval(parse(functionCall))
    end

    mira = Mira(analyses)

    return mira
end
            
