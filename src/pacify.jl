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
    ydummies = Matrix{Float64}(undef, length(y), length(yLevels))

    # Convert the variable to dummy variables
    for q in eachindex(yLevels)
        ydummies[:, q] = y .== yLevels[q]
    end

    return ydummies
end

# The removelindeps! function includes a ! as it updates X in place
function removelindeps!(
    X::Matrix{Float64},
    y::AbstractArray,
    where_y::Vector{Bool},
    wherecount::Int
    )

    # If all y-values are missing, stop now
    if wherecount == length(where_y)
        return
    end

    # Grab observed predictor data
    Xₒ = Matrix{Float64}(X[.!where_y, :])
    
    # If y is categorical
    if y isa CategoricalArray || nonmissingtype(eltype(y)) <: Union{AbstractString, CategoricalValue}
        # Grab observed y-values
        yₒ = y[.!where_y]

        # Convert y to dummy variables (as floats)
        mapping = Dict(levels(yₒ)[i] => i-1 for i in eachindex(levels(yₒ)))
        yₒ = Vector{Float64}([mapping[v] for v in yₒ])
    else
        # Grab observed y-values (as floats)
        yₒ = Vector{Float64}(y[.!where_y])
    end

    # If the variance of observed y-values falls below the allowed threshold, delete all predictors and stop now
    if var(yₒ) < 1e-4
        X = X[:, []]
        return
    end

    # Define vector of predictors to keep as those with sufficiently high variance and sufficiently low correlation with the observed y-values
    keep = var.(eachcol(Xₒ)) .> 1e-4 .&& cor.(eachcol(Xₒ), [yₒ]) .< 0.99

    # Get the total number of predictors currently being kept
    keepsum = sum(keep)

    # If there are 0 or 1 remaining, stop now
    if keepsum < 2
        X = X[:, keep]
        return
    end

    # Otherwise, calculate the correlation matrix of the remaining predictorsa
    xcors = cor(Xₒ)
    nxcors = xcors[findall(keep), findall(keep)]

    # Calculate the eigenvalues and eigenvectors of the correlation matrix and sort them
    eigencors = eigen(nxcors)
    eigvalsorder = sortperm(abs.(eigencors.values), rev = true)
    sortedeigvals = eigencors.values[eigvalsorder]
    sortedeigvecs = eigencors.vectors[:, eigvalsorder]

    # While the largest eigenvalue is more than 10_000 times larger than the smallest 
    while sortedeigvals[keepsum] / sortedeigvals[1] < 1e-4
        # Remove the predictor contributing most to the eigenvector with the smallest eigenvalue
        w = sortperm(abs.(sortedeigvecs[:, keepsum]), rev = true)[1]
        keep[findall(keep)[w]] = false
        nxcors = xcors[findall(keep), findall(keep)]
        keepsum -= 1
        eigencors = eigen(nxcors)
        eigvalsorder = sortperm(abs.(eigencors.values), rev = true)
        sortedeigvals = eigencors.values[eigvalsorder]
        sortedeigvecs = eigencors.vectors[:, eigvalsorder]
    end

    # Remove the predictors that are not being kept
    X = X[:, keep]
end

function pacifyworkingdata(workingdata::AxisVector)
    categoricalcolumns = Vector{String}([])

    for i in eachindex(workingdata)
        if workingdata[i][1] isa CategoricalArray || nonmissingtype(eltype(workingdata[i][1])) <: Union{AbstractString, CategoricalValue}
            push!(categoricalcolumns, axes(workingdata)[1][i])
        end
    end

    workingdatalevels = AxisArray(
        [collect(levels(workingdata[yvar][1])) for yvar in categoricalcolumns],
        categoricalcolumns
    )

    workingdatapacified = AxisArray(
        [[pacifyworkingdata(workingdata[yvar][j], workingdatalevels[yvar]) for j in eachindex(workingdata[yvar])] for yvar in categoricalcolumns],
        categoricalcolumns
    )

    return workingdatapacified, workingdatalevels
end

function pacifyworkingdata(workingdata::AbstractVector, levels::AbstractVector)
    contrastsmatrix = contrasts_matrix(PolynomialCoding(), 1, length(levels))

    workingdatapacified = Matrix{Float64}(undef, length(workingdata), size(contrastsmatrix, 2))

    for i in Base.axes(workingdatapacified, 2)
        for j in eachindex(levels)
            workingdatapacified[workingdata .== levels[j], i] .= contrastsmatrix[j, i]
        end
    end

    # Standardise everything
    for i in Base.axes(workingdatapacified, 2)
        workingdatapacified[:, i] = standardize(UnitRangeTransform, workingdatapacified[:, i])
    end

    return workingdatapacified
end