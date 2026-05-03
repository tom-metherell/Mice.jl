"""
    findmissings(data)

Returns an AxisVector of boolean vectors describing the locations of missing data in each
column of the provided data table.
"""
function findmissings(data::T)::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}} where{T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    imputewhere = AxisArray(
        [Vector{Bool}(ismissing.(getcolumn(data, i))) for i in columnnames(data)],
        collect(string.(columnnames(data)))
    )

    return imputewhere
end

function makemonotonesequence(imputewhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}})::Vector{String}
    # Sort the data frame names vector by missingness
    visitsequence = axes(imputewhere)[1][sortperm(sum.(imputewhere))][sum(sum.(imputewhere) .== 0)+1:end]

    return visitsequence
end

"""
    makemethods(data)

Returns an AxisVector of strings defining the method by which each variable in `data`
should be imputed in the `mice()` function. The default method is predictive mean matching
(pmm).
"""
function makemethods(data::T) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    names = collect(string.(columnnames(data)))
    no = length(names)

    # Use pmm for all variables by default
    methods = AxisArray(
        fill("pmm", no),
        names
    )

    return methods
end

"""
    makepredictormatrix(data)

Returns an AxisMatrix of integers defining the predictors for each variable in `data`.
The variables to be predicted are on the rows, and the predictors are on the columns.
The default is to use all variables as predictors for all other variables (i.e. all
1s except for the diagonal, which is 0).
"""
function makepredictormatrix(data::T) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    names = collect(string.(columnnames(data)))
    no = length(names)

    # Initialise the predictor matrix with 1s
    predictormatrix = AxisArray(
        fill(1, no, no),
        names,
        names
    )
    
    # Set the diagonal to 0
    for i in 1:no
        predictormatrix[i, i] = 0
    end

    return predictormatrix
end

function initialiseworkingdata(
    data::T,
    imputewhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}},
    m::Int,
    visitsequence::Vector{String},
    methods::AxisArray{String, 1, Vector{String}},
    predictormatrix::AxisArray{Int, 2, Matrix{Int}}
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Select only variables that will predict another or will be imputed themselves
    predictors = axes(predictormatrix)[2][findall(col -> any(x -> x != 0, col), eachcol(predictormatrix))]
    imputed = visitsequence[findall(x -> methods[x] ≠ "", visitsequence)]
    wdvars = collect(string.(columnnames(data)))[in.(collect(string.(columnnames(data))), Ref(vcat(imputed, predictors)))]

    # Initialise working data vectors
    workingdata = AxisArray(
        [[getcolumn(data, Symbol(var)) isa CategoricalArray ? CategoricalArray(getcolumn(data, Symbol(var))) : Vector{eltype(getcolumn(data, Symbol(var)))}(getcolumn(data, Symbol(var))) for i in 1:m] for var in wdvars],
        wdvars
    )

    for var ∈ wdvars
        # If the variable is to be imputed
        if var ∈ imputed
            # Get locations of data to be imputed in column
            where_y = imputewhere[var]

            # Count data to be imputed in column
            wherecount = sum(where_y)

            # For each imputation
            for j in 1:m
                # Initialise using a random sample from the observed data
                workingdata[var][j][where_y] = IMPUTERS["sample"].f(
                    workingdata[var][j],
                    nothing,
                    where_y,
                    wherecount,
                    var,
                    0,
                    j,
                    String[]
                )
            end

            # Convert to non-missing type
            if workingdata[var][1] isa CategoricalArray
                workingdata[var] = [CategoricalArray{nonmissingtype(eltype(workingdata[var][1]))}(workingdata[var][j]) for j in 1:m]
            else
                workingdata[var] = convert(Vector{Vector{nonmissingtype(eltype(workingdata[var][1]))}}, workingdata[var])
            end
        end

        # Fallback in case there is a variable with no missing data that allows missing values
        if Missing <: eltype(workingdata[var][1])
            workingdata[var] = convert(Vector{Vector{nonmissingtype(eltype(workingdata[var][1]))}}, workingdata[var])
        end
    end

    return workingdata
end

function initialiseworkingdata(
    data::T,
    imputations::Vector{Matrix},
    imputewhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}},
    m::Int,
    visitsequence::Vector{String},
    methods::AxisArray{String, 1, Vector{String}},
    predictormatrix::AxisArray{Int, 2, Matrix{Int}}
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Comments as above
    predictors = axes(predictormatrix)[2][findall(col -> any(x -> x != 0, col), eachcol(predictormatrix))]
    imputed = visitsequence[findall(x -> methods[x] ≠ "", visitsequence)]
    wdvars = collect(string.(columnnames(data)))[in.(collect(string.(columnnames(data))), Ref(vcat(imputed, predictors)))]

    # Initialise working data vectors
    workingdata = AxisArray(
        [[getcolumn(data, Symbol(var)) isa CategoricalArray ? CategoricalArray(getcolumn(data, Symbol(var))) : Vector{eltype(getcolumn(data, Symbol(var)))}(getcolumn(data, Symbol(var))) for i in 1:m] for var in wdvars],
        wdvars
    )

    for var ∈ wdvars
        # If the variable is to be imputed
        if var ∈ imputed
            # For each imputation
            for j in 1:m
                # Initialise using the provided imputations
                if workingdata[var][1] isa CategoricalArray || nonmissingtype(eltype(workingdata[var][1])) <: CategoricalValue
                    workingdata[var][j][imputewhere[var]] = CategoricalArray(imputations[findfirst(visitsequence .== var)][:, j], levels = levels(workingdata[var][j]))
                else
                    workingdata[var][j][imputewhere[var]] = imputations[findfirst(visitsequence .== var)][:, j]
                end
            end

            # Convert to non-missing type
            if workingdata[var][1] isa CategoricalArray
                for j in eachindex(workingdata[var])
                    workingdata[var][j] = convert(CategoricalArray{nonmissingtype(eltype(workingdata[var][1]))}, workingdata[var][j])
                end
            else
                workingdata[var] = convert(Vector{Vector{nonmissingtype(eltype(workingdata[var][1]))}}, workingdata[var])
            end
        end

        if Missing <: eltype(workingdata[var][1])
            workingdata[var] = convert(Vector{Vector{nonmissingtype(eltype(workingdata[var][1]))}}, workingdata[var])
        end
    end

    return workingdata
end

# Allow US spelling of initialise
const initializeworkingdata = initialiseworkingdata

function initialisetraces(
    visitsequence::Vector{String},
    iter::Int,
    m::Int
    )

    traces = [Matrix{Float64}(undef, iter, m) for _ = eachindex(visitsequence)]

    return traces
end

# Allow US spelling of initialise
const initializetraces = initialisetraces