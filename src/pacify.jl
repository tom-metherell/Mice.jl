# The fillXMissings! function includes a ! as it updates X in place
function fillXMissings!(
    X::T,
    imputeWhere::AxisVector{Vector{Bool}},
    predictors::Vector{String},
    visitSequence::Vector{String},
    imputations::Vector{Matrix},
    j::Int
    ) where {T}
    istable(X) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # For each predictor
    for k in predictors
        # Find its position in the visit sequence
        kVS = findfirst(visitSequence .== k)

        # Find the positions of data to be imputed in the predictor column
        whereX = imputeWhere[k]

        # If there are any missing data
        if any(whereX)
            # Fill them with the imputed values
            X[Symbol(k)][whereX] = imputations[kVS][:, j]
        end
    end
end

# The pacify! function includes a ! as it updates loggedEvents in place
function pacify!(
    X::U,
    predictors::Vector{String},
    loggedEvents::Vector{String},
    iterCounter::Int,
    yVar::AbstractString,
    j::Int
    ) where{U}
    istable(X) || throw(ArgumentError("Data not provided as a Tables.jl table."))

    # Initialise vector of categorical predictors
    categoricalPredictors = Vector{Symbol}([])

    # For each predictor
    for xVar in Symbol.(predictors)
        # Grab the predictor data
        x = getcolumn(X, xVar)

        # If the data are categorical (either CategoricalArray or a vector of strings)
        if x isa CategoricalArray || nonmissingtype(eltype(x)) <: AbstractString
            # If there is more than one level
            if length(levels(x)) > 1
                # Add this variable to the list of categorical predictors
                push!(categoricalPredictors, xVar)
            else
                # Otherwise, drop this variable
                X = columntable(NamedTuple{Tuple(setdiff(columnnames(X), [xVar]))}([X[c] for c in setdiff(columnnames(X), [xVar])]))
                predictors = predictors[predictors .≠ string(xVar)]

                # Log that this predictor has been dropped
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: predictor $xVar dropped because of zero variance.")
            end
        end
    end

    # Define the model frame
    mf = ModelFrame(term(0) ~ sum(term.(predictors)), X)

    # Set contrast coding for categorical predictors to orthogonal polynomial
    setcontrasts!(mf, Dict([xVar => PolynomialCoding() for xVar in categoricalPredictors]))

    # Produce the model matrix
    X = ModelMatrix(mf).m[:, 2:end]

    # Standardise everything
    for i in Base.axes(X, 2)
        X[:, i] = standardize(UnitRangeTransform, X[:, i])
    end

    return X
end

# The PolynomialCoding type is used to specify orthogonal polynomial contrast coding
mutable struct PolynomialCoding <: AbstractContrasts
end

# This extension of the contrasts_matrix function from StatsModels.jl allows orthogonal polynomial contrast coding
function contrasts_matrix(C::PolynomialCoding, _, n)
    X = reduce(hcat, [((1:n) .- mean(1:n)) .^ i for i in 0:n-1])
    qrX = qr(X)
    Z = qrX.Q * Diagonal(qrX.R)
    for i in Base.axes(Z, 2)
        Z[:, i] = Z[:, i] ./ sqrt(sum(Z[:, i].^2))
    end
    return Z[:, 2:end]
end

# This extension of the termnames function from StatsModels.jl names the new dummy variables correctly
function termnames(C::PolynomialCoding, levels::AbstractVector, _::Integer)
    return Vector{String}([".^$i" for i in 1:length(levels)])
end

function pacify(y::AbstractArray)
    # Grab the levels of the variable
    yLevels = levels(y)

    # Initialise the dummy matrix
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))

    # Convert the variable to dummy variables
    for q in eachindex(yLevels)
        yDummies[:, q] = y .== yLevels[q]
    end

    return yDummies
end

# The removeLinDeps! function includes a ! as it updates X in place
function removeLinDeps!(
    X::Matrix{Float64},
    y::AbstractArray,
    whereY::Vector{Bool},
    whereCount::Int
    )

    # If all y-values are missing, stop now
    if whereCount == length(whereY)
        return
    end

    # Grab observed predictor data
    Xₒ = Matrix{Float64}(X[.!whereY, :])
    
    # If y is categorical
    if y isa CategoricalArray || nonmissingtype(eltype(y)) <: AbstractString
        # Grab observed y-values
        yₒ = y[.!whereY]

        # Convert y to dummy variables (as floats)
        mapping = Dict(levels(yₒ)[i] => i-1 for i in eachindex(levels(yₒ)))
        yₒ = Vector{Float64}([mapping[v] for v in yₒ])
    else
        # Grab observed y-values (as floats)
        yₒ = Vector{Float64}(y[.!whereY])
    end

    # If the variance of observed y-values falls below the allowed threshold, delete all predictors and stop now
    if var(yₒ) < 1e-4
        X = X[:, []]
        return
    end

    # Define vector of predictors to keep as those with sufficiently high variance and sufficiently low correlation with the observed y-values
    keep = var.(eachcol(Xₒ)) .> 1e-4 .&& cor.(eachcol(Xₒ), [yₒ]) .< 0.99

    # Get the total number of predictors currently being kept
    keepSum = sum(keep)

    # If there are 0 or 1 remaining, stop now
    if keepSum < 2
        X = X[:, keep]
        return
    end

    # Otherwise, calculate the correlation matrix of the remaining predictorsa
    xCors = cor(Xₒ)
    nxCors = xCors[findall(keep), findall(keep)]

    # Calculate the eigenvalues and eigenvectors of the correlation matrix and sort them
    eigenCors = eigen(nxCors)
    eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
    sortedeigvals = eigenCors.values[eigvalsorder]
    sortedeigvecs = eigenCors.vectors[:, eigvalsorder]

    # While the largest eigenvalue is more than 10_000 times larger than the smallest 
    while sortedeigvals[keepSum] / sortedeigvals[1] < 1e-4
        # Remove the predictor contributing most to the eigenvector with the smallest eigenvalue
        w = sortperm(abs.(sortedeigvecs[:, keepSum]), rev = true)[1]
        keep[findall(keep)[w]] = false
        nxCors = xCors[findall(keep), findall(keep)]
        keepSum -= 1
        eigenCors = eigen(nxCors)
        eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
        sortedeigvals = eigenCors.values[eigvalsorder]
        sortedeigvecs = eigenCors.vectors[:, eigvalsorder]
    end

    # Remove the predictors that are not being kept
    X = X[:, keep]
end

function pacifyWorkingData(workingData::AxisVector{Vector})
    categoricalColumns = Vector{String}([])

    for i in eachindex(workingData)
        if workingData[i][1] isa CategoricalArray || nonmissingtype(eltype(workingData[i][1])) <: Union{AbstractString, CategoricalValue}
            push!(categoricalColumns, axes(workingData)[1][i])
        end
    end

    workingDataLevels = AxisArray(
        [levels(workingData[yVar][1]) for yVar in categoricalColumns],
        categoricalColumns
    )

    workingDataPacified = AxisArray(
        [[pacifyWorkingData(workingData[yVar][j], workingDataLevels[yVar]) for j in eachindex(workingData[yVar])] for yVar in categoricalColumns],
        categoricalColumns
    )

    return workingDataPacified, workingDataLevels
end

function pacifyWorkingData(workingData::AbstractVector, levels::Vector)
    contrastsMatrix = contrasts_matrix(PolynomialCoding(), 1, length(levels))

    workingDataPacified = Matrix{Float64}(undef, length(workingData), size(contrastsMatrix, 2))

    for i in Base.axes(workingDataPacified, 2)
        for j in eachindex(levels)
            workingDataPacified[workingData .== levels[j], i] .= contrastsMatrix[j, i]
        end
    end

    # Standardise everything
    for i in Base.axes(workingDataPacified, 2)
        workingDataPacified[:, i] = standardize(UnitRangeTransform, workingDataPacified[:, i])
    end

    return workingDataPacified
end