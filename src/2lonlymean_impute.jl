function secondlevelonlymean_impute(
    ydata::Vector{Float64},
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    types::Vector{Int},
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String};
    unusedkwargs...
    )

    if wherecount == 0
        return Float64[]
    end

    classcols = findall(types .== -2)

    if isempty(classcols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end

    classmatrix = Matrix{Float64}(X[:, classcols])
    classkeys = [Tuple(classmatrix[r, c] for c in axes(classmatrix, 2)) for r in axes(classmatrix, 1)]
    classlevels = unique(classkeys)
    nclasses = length(classlevels)
    classmap = Dict(level => idx for (idx, level) in enumerate(classlevels))
    gf_full = [classmap[key] for key in classkeys]
    gf = gf_full[.!where_y]

    yₒ = ydata[.!where_y]

    if length(unique(gf)) < nclasses
        throw(ArgumentError("Two-level imputation requires at least one observed outcome per class."))
    end

    # Get class assignments for missing observations
    gfₘ = gf_full[where_y]

    # Calculate class means
    classstats = [mean(yₒ[gf .== class]) for class in 1:nclasses]

    # Impute with class means
    return [classstats[gfₘ[i]] for i in 1:wherecount]
end

const SECOND_LEVEL_ONLY_MEAN_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; kwargs...) -> begin
    secondlevelonlymean_impute(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; kwargs...)
end; twolevel = true)

registerimputer!("2lonly.mean", SECOND_LEVEL_ONLY_MEAN_IMPUTER)
