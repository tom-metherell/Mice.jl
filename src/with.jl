"""
    complete(
        mids::Mids,
        imputation::Int
        )

Produces a data table with missings replaced with imputed values from a
multiply imputed dataset (`Mids`) object.

The Mids object must be supplied first.

The `imputation` argument is an integer identifying which specific imputation
is to be used to fill in the missing values.
"""
function complete(
    mids::Mids,
    imputation::Int
    )

    # Get the observed data
    data = deepcopy(mids.data)

    # For each variable
    for i in eachindex(mids.visitSequence)
        # If it was imputed
        if isassigned(mids.imputations, i)
            # Replace missings with imputed values
            var = mids.visitSequence[i]
            data[ismissing.(data[:, var]), var] = mids.imputations[i][:, imputation]
        end
    end

    return data
end

"""
    complete(
        mids::Mids,
        action::String
        )

Summarises the outputs of all imputations in a multiply imputed dataset (`Mids`).

The Mids object must be supplied first.

The `action` argument is a string identifying what format the output should take.
If specified as "long", the function will return a single data table, containing
the results of each imputation in succession with an identifier (`imp`).
If specified as "list", the function will return a vector of individual
data tables.
"""
function complete(
    mids::Mids,
    action::String
    )

    # Wrong action specified
    if !(action ∈ ["list", "long"])
        throw(ArgumentError("Action not defined. Valid arguments: \"list\", \"long\""))
    else
        # Detect type of data object
        T = typeof(mids.data)

        # Initialise output
        data = Vector{T}(undef, mids.m)

        # For each imputation
        for i in 1:mids.m
            # Get the observed data
            data[i] = deepcopy(mids.data)
            # For each variable
            for j in eachindex(mids.visitSequence)
                # If it was imputed
                if isassigned(mids.imputations, j)
                    # Replace missings with imputed values
                    var = mids.visitSequence[j]
                    data[i][ismissing.(data[i][:, var]), var] = mids.imputations[j][:, i]
                end
            end
        end
        
        if action == "list"
            return data
        elseif action == "long"
            # Concatenate the data tables
            for i in eachindex(data)
                T = typeof(data[i])
                rt = rowtable(data[i])
                rt = [merge(r, (imp = i,)) for r in rt]
                data[i] = T(rt)
            end
            return vcat(data...)
        end
    end
end

"""
    Mira

A multiply imputed repeated analyses object.

The analyses are stored as a vector of analyses of individual imputations.
"""
struct Mira
    analyses::Vector
end

"""
    with(
        mids::Mids,
        func::Function
        )

Conducts repeated analyses of a multiply imputed dataset (`Mids`).

The function takes two arguments: firstly the `Mids` object itself,
then a function (`func`). The function should take the form
`data -> analysisFunction(arguments, data, moreArguments...)`,
where `data` represents the position of the data argument in the function.

For example: `with(mids, data -> lm(@formula(y ~ x1 + x2), data))`
"""
function with(
    mids::Mids, 
    func::Function
    )

    # Initialise output
    analyses = Vector{Any}(undef, mids.m)
    
    # Fill in the missing values
    datalist = complete(mids, "list")

    # For each imputation
    for i in eachindex(datalist)
        # Conduct the analysis on the completed data
        data = datalist[i]
        analyses[i] = func(data)
    end

    return Mira(analyses)
end

export complete, Mira, with