"""
    findMissings(data)

Returns an AxisVector of boolean vectors describing the locations of missing data in each
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

function makeMonotoneSequence(imputeWhere::AxisVector{Vector{Bool}})
    # Sort the data frame names vector by missingness
    visitSequence = axes(imputeWhere)[1][sortperm(sum.(imputeWhere))][sum(sum.(imputeWhere) .== 0)+1:end]

    return visitSequence
end

"""
    makeMethods(data)

Returns an AxisVector of strings defining the method by which each variable in `data`
should be imputed in the `mice()` function. The default method is predictive mean matching
(pmm).
"""
function makeMethods(data::T) where {T}
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
    makePredictorMatrix(data)

Returns an AxisMatrix of integers defining the predictors for each variable in `data`.
The variables to be predicted are on the rows, and the predictors are on the columns.
The default is to use all variables as predictors for all other variables (i.e. all
1s except for the diagonal, which is 0).
"""
function makePredictorMatrix(data::T) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    names = collect(string.(columnnames(data)))
    no = length(names)

    # Initialise the predictor matrix with 1s
    predictorMatrix = AxisArray(
        fill(1, no, no),
        names,
        names
    )
    
    # Set the diagonal to 0
    for i in 1:no
        predictorMatrix[i, i] = 0
    end

    return predictorMatrix
end

function initialiseWorkingData(
    data::T,
    imputeWhere::AxisVector{Vector{Bool}},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisVector{String},
    predictorMatrix::AxisMatrix{Int}
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Select only variables that will predict another or will be imputed themselves
    predictors = axes(predictorMatrix)[2][findall(col -> any(x -> x != 0, col), eachcol(predictorMatrix))]
    imputed = visitSequence[findall(x -> methods[x] ≠ "", visitSequence)]
    wdVars = collect(string.(columnnames(data)))[in.(collect(string.(columnnames(data))), Ref(vcat(imputed, predictors)))]

    # Initialise working data vectors
    workingData = AxisArray(
        [[getcolumn(data, Symbol(var)) isa CategoricalArray ? CategoricalArray(getcolumn(data, Symbol(var))) : Vector{eltype(getcolumn(data, Symbol(var)))}(getcolumn(data, Symbol(var))) for i in 1:m] for var in wdVars],
        wdVars
    )

    for var ∈ wdVars
        # If the variable is to be imputed
        if var ∈ imputed
            # Get locations of data to be imputed in column
            whereY = imputeWhere[var]

            # Count data to be imputed in column
            whereCount = sum(whereY)

            # For each imputation
            for j in 1:m
                # Initialise using a random sample from the observed data
                workingData[var][j][whereY] = sampleImpute!(workingData[var][j][.!whereY], whereCount)
            end

            # Convert to non-missing type
            if workingData[var][1] isa CategoricalArray
                for j in eachindex(workingData[var])
                    workingData[var][j] = convert(CategoricalArray{nonmissingtype(eltype(workingData[var][1]))}, workingData[var][j])
                end
            else
                workingData[var] = convert(Vector{Vector{nonmissingtype(eltype(workingData[var][1]))}}, workingData[var])
            end
        end
    end

    return workingData
end

function initialiseWorkingData(
    data::T,
    imputations::Vector{Matrix},
    imputeWhere::AxisVector{Vector{Bool}},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisVector{String},
    predictorMatrix::AxisMatrix{Int}
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Comments as above
    predictors = axes(predictorMatrix)[2][findall(col -> any(x -> x != 0, col), eachcol(predictorMatrix))]
    imputed = visitSequence[findall(x -> methods[x] ≠ "", visitSequence)]
    wdVars = collect(string.(columnnames(data)))[in.(collect(string.(columnnames(data))), Ref(vcat(imputed, predictors)))]

    # Initialise working data vectors
    workingData = AxisArray(
        [[getcolumn(data, Symbol(var)) isa CategoricalArray ? CategoricalArray(getcolumn(data, Symbol(var))) : Vector{eltype(getcolumn(data, Symbol(var)))}(getcolumn(data, Symbol(var))) for i in 1:m] for var in wdVars],
        wdVars
    )

    for var ∈ wdVars
        # If the variable is to be imputed
        if var ∈ imputed
            # For each imputation
            for j in 1:m
                # Initialise using the provided imputations
                workingData[var][j][imputeWhere[var]] = imputations[findfirst(visitSequence .== var)][:, j]
            end

            # Convert to non-missing type
            if workingData[var][1] isa CategoricalArray
                for j in eachindex(workingData[var])
                    workingData[var][j] = convert(CategoricalArray{nonmissingtype(eltype(workingData[var][1]))}, workingData[var][j])
                end
            else
                workingData[var] = convert(Vector{Vector{nonmissingtype(eltype(workingData[var][1]))}}, workingData[var])
            end
        end
    end

    return workingData
end

# Allow US spelling of initialise
const initializeWorkingData = initialiseWorkingData

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