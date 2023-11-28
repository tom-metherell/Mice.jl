"""
    findMissings(data)

Returns a named vector of boolean vectors describing the locations of missing data in each
column of the provided data table.
"""
function findMissings(data::T) where{T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    imputeWhere = AxisArray(
        [Vector{Bool}(ismissing.(getcolumn(data, i))) for i in columnnames(data)],
        collect(string.(columnnames(data)))
    )

    return imputeWhere
end

function makeMonotoneSequence(imputeWhere::AxisArray{Vector{Bool}})
    missingness = sum.(imputeWhere)

    numberOfCompletes = sum(sum.(imputeWhere) .== 0)

    missingness = sort(missingness)

    # Sort the data frame names vector by missingness
    visitSequence = names(missingness)[1][numberOfCompletes+1:end]

    return visitSequence
end

"""
    makeMethods(data)

Returns a named vector of strings defining the method by which each variable in `data`
should be imputed in the `mice()` function. The default method is predictive mean matching
(pmm).
"""
function makeMethods(data::T) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Use pmm for all variables by default
    methods = AxisArray(
        Vector{String}(fill("pmm", length(columns(data)))),
        collect(string.(columnnames(data)))    
    )

    return methods
end

"""
    makePredictorMatrix(data)

Returns a named matrix of booleans defining the predictors for each variable in `data`.
The variables to be predicted are on the rows, and the predictors are on the columns.
The default is to use all variables as predictors for all other variables (i.e. all
1s except for the diagonal, which is 0).
"""
function makePredictorMatrix(data::T) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Initialise the predictor matrix with 1s
    predictorMatrix = AxisArray(
        Matrix{Bool}(fill(1, length(Tables.columns(data)), length(Tables.columns(data)))),
        collect(string.(Tables.columnnames(data))),
        collect(string.(Tables.columnnames(data)))
    )
    
    # Set the diagonal to 0
    for i in 1:length(Tables.columns(data))
        predictorMatrix[i, i] = 0
    end

    return predictorMatrix
end

function initialiseImputations(
    data::T,
    imputeWhere::AxisVector{Bool},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisVector{String}
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Initialise vector of imputations matrices
    imputations = Vector{Matrix}(undef, length(visitSequence))

    # For each variable
    for i in eachindex(visitSequence)

        # Grab the variable name
        yVar = visitSequence[i]

        # If the variable is to be imputed
        if methods[yVar] != ""

            # Grab the variable data
            y = getcolumn(data, Symbol(yVar))
            
            # Get locations of data to be imputed in column
            whereY = imputeWhere[yVar]

            # Count data to be imputed in column
            whereCount = sum(whereY)

            # Initialise imputation matrix
            imputations[i] = Matrix{nonmissingtype(eltype(y))}(undef, whereCount, m)

            # For each imputation
            for j in 1:m
                # Initialise using a random sample from the observed data
                imputations[i][:, j] = sampleImpute!(y, whereY, whereCount)
            end
        end
    end

    return imputations
end

# Allow US spelling of initialise
const initializeImputations = initialiseImputations

function initialiseTraces(
    visitSequence::Vector{String},
    iter::Int,
    m::Int
    )

    traces = [Matrix{Float64}(undef, iter, m) for _ = eachindex(visitSequence)]

    return traces
end

# Allow US spelling of initialise
const initializeTraces = initialiseTraces

export findMissings, makeMethods, makePredictorMatrix