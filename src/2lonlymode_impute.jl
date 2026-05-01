function secondlevelonlymode_impute(
    ydata::AbstractArray,
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
        return Vector{eltype(ydata)}(undef, 0)
    end

    classcols = findall(types .== -2)
    if isempty(classcols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end

    # Extract class information
    classmatrix = Matrix{Float64}(X[:, classcols])
    classkeys = [Tuple(classmatrix[r, c] for c in axes(classmatrix, 2)) for r in axes(classmatrix, 1)]
    classlevels = unique(classkeys)
    nclasses = length(classlevels)
    classmap = Dict(level => idx for (idx, level) in enumerate(classlevels))
    gf_full = [classmap[key] for key in classkeys]

    # Calculate mode for each class
    classmodes = Vector{eltype(ydata)}(undef, nclasses)
    for i in 1:nclasses
        classindex = findall(gf_full .== i)
        obsindex = classindex[.!where_y[classindex]]
        if !isempty(obsindex)
            classvalues = ydata[obsindex]
            classmodes[i] = mode(classvalues)
        end
    end

    # Map modes back to original missing rows
    gfₘ = gf_full[where_y]
    return [classmodes[gfₘ[i]] for i in 1:wherecount]
end

const SECOND_LEVEL_ONLY_MODE_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; kwargs...) -> begin
    secondlevelonlymode_impute(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; kwargs...)
end; twolevel = true)

registerimputer!("2lonly.mode", SECOND_LEVEL_ONLY_MODE_IMPUTER)
