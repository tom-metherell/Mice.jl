# The PolynomialCoding type is used to specify orthogonal polynomial contrast coding
mutable struct PolynomialCoding <: AbstractContrasts
end

# This extension of the contrasts_matrix function from StatsModels.jl allows orthogonal polynomial contrast coding
function contrasts_matrix(C::PolynomialCoding, _, n)
    X = reduce(hcat, [((1:n) .- mean(1:n)) .^ i for i ∈ 0:n-1])
    qrX = qr(X)
    Z = qrX.Q * Diagonal(qrX.R)
    for i ∈ Base.axes(Z, 2)
        Z[:, i] = Z[:, i] ./ sqrt(sum(Z[:, i].^2))
    end
    return Z[:, 2:end]
end

# This extension of the termnames function from StatsModels.jl names the new dummy variables correctly
function termnames(C::PolynomialCoding, levels::AbstractVector, _::Integer)
    return Vector{String}([".^$i" for i ∈ 1:length(levels)])
end

function pacify(y::AbstractArray)
    # Grab the levels of the variable
    yLevels = levels(y)

    # Initialise the dummy matrix
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))

    # Convert the variable to dummy variables
    for q ∈ eachindex(yLevels)
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
    if y isa CategoricalArray || nonmissingtype(eltype(y)) <: Union{AbstractString, CategoricalValue}
        # Grab observed y-values
        yₒ = y[.!whereY]

        # Convert y to dummy variables (as floats)
        mapping = Dict(levels(yₒ)[i] => i-1 for i ∈ eachindex(levels(yₒ)))
        yₒ = Vector{Float64}([mapping[v] for v ∈ yₒ])
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

    for i ∈ eachindex(workingData)
        if workingData[i][1] isa CategoricalArray || nonmissingtype(eltype(workingData[i][1])) <: Union{AbstractString, CategoricalValue}
            push!(categoricalColumns, axes(workingData)[1][i])
        end
    end

    workingDataLevels = AxisArray(
        [levels(workingData[yVar][1]) for yVar ∈ categoricalColumns],
        categoricalColumns
    )

    workingDataPacified = AxisArray(
        [[pacifyWorkingData(workingData[yVar][j], workingDataLevels[yVar]) for j ∈ eachindex(workingData[yVar])] for yVar ∈ categoricalColumns],
        categoricalColumns
    )

    return workingDataPacified, workingDataLevels
end

function pacifyWorkingData(workingData::AbstractVector, levels::Vector)
    contrastsMatrix = contrasts_matrix(PolynomialCoding(), 1, length(levels))

    workingDataPacified = Matrix{Float64}(undef, length(workingData), size(contrastsMatrix, 2))

    for i ∈ Base.axes(workingDataPacified, 2)
        for j ∈ eachindex(levels)
            workingDataPacified[workingData .== levels[j], i] .= contrastsMatrix[j, i]
        end
    end

    # Standardise everything
    for i ∈ Base.axes(workingDataPacified, 2)
        workingDataPacified[:, i] = standardize(UnitRangeTransform, workingDataPacified[:, i])
    end

    return workingDataPacified
end