function secondlevelonlynorm_impute!(
    ydata::Vector{Float64},
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    types::Vector{Int},
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String};
    ridge::Float64 = 1e-4,
    unusedkwargs...
    )

    if wherecount == 0
        return Float64[]
    end

    return _imputationlevel2(
        ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents,
        (yₒ, X2l, wy2l, wherecount2l, yvar, itercounter, j, loggedevents; ridge=ridge, kwargs...) -> 
            norm_impute!(yₒ, X2l, wy2l, wherecount2l, yvar, itercounter, j, loggedevents; ridge=ridge, kwargs...);
        ridge=ridge,
        unusedkwargs...
    )
end

const SECOND_LEVEL_ONLY_NORM_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; ridge::Float64 = 1e-4, kwargs...) -> begin
    secondlevelonlynorm_impute!(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; ridge = ridge, kwargs...)
end; twolevel = true)

registerimputer!("2lonly.norm", SECOND_LEVEL_ONLY_NORM_IMPUTER)

# Helper function for level-2 only imputation methods
# Aggregates level-1 data to level-2, calls the specified imputation method, 
# and maps results back to original rows

function _imputationlevel2(
    ydata::AbstractArray,
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    types::Vector{Int},
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String},
    imputationmethod::Function;
    ridge::Float64 = 1e-4,
    kwargs...
    )

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

    # Check for partial missing level-2 data
    # (level-2 data should be constant within each class)
    for class in 1:nclasses
        classindex = findall(gf_full .== class)
        obsindex = classindex[.!where_y[classindex]]
        misindex = classindex[where_y[classindex]]
        
        if !isempty(obsindex) && !isempty(misindex)
            clusterids = join(classlevels[class], ", ")
            throw(ArgumentError("Two-level imputation found partially missing level-2 data in cluster $clusterids. Use 2lonly.mean or 2lonly.mode to fix inconsistencies."))
        end
    end

    # Aggregate level-1 predictors to level-2 by class means
    randomcols = findall(types .== 2)
    X2lagg = Matrix{Float64}(undef, nclasses, length(randomcols))
    for i in 1:nclasses
        classindex = findall(gf_full .== i)
        if !isempty(classindex)
            X2lagg[i, :] = vec(mean(X[classindex, randomcols], dims=1))
        end
    end

    # Create level-2 missing indicator
    # where_y2l[i] = true if all values in class i are missing
    where_y2l = Vector{Bool}(undef, nclasses)
    for i in 1:nclasses
        where_y2l[i] = all(where_y[findall(gf_full .== i)])
    end
    
    # Get observed y values at level-2.
    # Numeric data are averaged; non-numeric data use the class mode.
    ytype = nonmissingtype(eltype(ydata))
    ynumeric = ytype <: Real
    y2l = ynumeric ? Vector{Float64}(undef, nclasses) : Vector{ytype}(undef, nclasses)
    for i in 1:nclasses
        classindex = findall(gf_full .== i)
        obsindex = classindex[.!where_y[classindex]]
        if !isempty(obsindex)
            classvalues = ydata[obsindex]
            y2l[i] = ynumeric ? Float64(mean(classvalues)) : mode(classvalues)
        end
    end

    # Call the specified imputation method at the aggregated level-2
    imps2l = imputationmethod(y2l, X2lagg, where_y2l, sum(where_y2l), yvar, itercounter, j, loggedevents; ridge=ridge, kwargs...)

    # Map the missing class ids to positions in imps2l.
    missingclasspositions = findall(where_y2l)
    missingclassmap = Dict(missingclass => idx for (idx, missingclass) in enumerate(missingclasspositions))

    # Map level-2 imputations back to original missing rows
    gfₘ = gf_full[where_y]
    return [imps2l[missingclassmap[gfₘ[i]]] for i in 1:wherecount]
end