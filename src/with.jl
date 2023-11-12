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
    listComplete(
        mids::Mids
        )

Summarises the outputs of all imputations in a multiply imputed dataset (`Mids`) as a list
of completed datasets.
"""
function listComplete(
    mids::Mids
    )

    # Detect type of data object
    T = typeof(mids.data)

    # Initialise output
    data = Vector{T}(undef, mids.m)

    # For each imputation
    for i in 1:mids.m
        # Get the observed data
        theseData = deepcopy(columntable(mids.data))
        # For each variable
        for j in eachindex(mids.visitSequence)
            # If it was imputed
            if isassigned(mids.imputations, j)
                # Replace missings with imputed values
                var = Symbol(mids.visitSequence[j])
                theseData[var][mids.imputeWhere[string(var)]] = mids.imputations[j][:, i]
                data[i] = T(theseData)
            end
        end
    end
        
    return data
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
    datalist = listComplete(mids)

    # For each imputation
    for i in eachindex(datalist)
        # Conduct the analysis on the completed data
        data = datalist[i]
        analyses[i] = func(data)
    end

    return Mira(analyses)
end

export complete, listComplete, Mira, with