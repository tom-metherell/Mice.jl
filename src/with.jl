"""
    complete(
        mids::Mids,
        imputation::Int
        )

Produces a DataFrame with missings replaced with imputed values from a
multiply imputed dataset (`Mids`) object.

The Mids object must be supplied first.

The `imputation` argument is an integer identifying which specific imputation
is to be used to fill in the missing values.
"""
function complete(
    mids::Mids,
    imputation::Int
    )

    data = mids.data

    for i in eachindex(mids.visitSequence)
        var = mids.visitSequence[i]
        data[ismissing.(data[:, var]), var] = data.imputations[i][:, imputation]
    end
end

"""
    complete(
        mids::Mids,
        action::String
    )

Summarises the outputs of all imputations in a multiply imputed dataset (`Mids`).

The Mids object must be supplied first.

The `action` argument is a string identifying what format the output should take.
If specified as "long", the function will return a single DataFrame, containing
the results of each imputation in succession with an identifier (`imp`).
If specified as "list", the function will return a vector of individual
DataFrames.
"""
function complete(
    mids::Mids,
    action::String
)
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
            for i in eachindex(data)
                insertcols!(data[i], 1, imp => i)
            end
            return vcat(data...)
        end
    end
end

"""
    Mira(analyses::Vector)

A multiply imputed repeated analyses object.

The analyses are stored as a vector of analyses of individual imputations.
"""
struct Mira(analyses::Vector)

"""
    with(
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
            
