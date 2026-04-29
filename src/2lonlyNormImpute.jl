# The secondLevelOnlyNormImpute! function includes a ! as it updates loggedEvents in place
function secondLevelOnlyNormImpute!(
    yData::Vector{Float64},
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    types::Vector{Int},
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    ridge::Float64 = 1e-4,
    unusedKwargs...
    )

    if whereCount == 0
        return Float64[]
    end

    return _imputationLevel2!(
        yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents,
        (yₒ, X2l, wy2l, whereCount2l, yVar, iterCounter, j, loggedEvents; ridge=ridge, kwargs...) -> 
            normImpute!(yₒ, X2l, wy2l, whereCount2l, yVar, iterCounter, j, loggedEvents; ridge=ridge, kwargs...);
        ridge=ridge,
        unusedKwargs...
    )
end

const SECOND_LEVEL_ONLY_NORM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge::Float64 = 1e-4, kwargs...) -> begin
    secondLevelOnlyNormImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge = ridge, kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.norm", SECOND_LEVEL_ONLY_NORM_IMPUTER)

# Helper function for level-2 only imputation methods
# Aggregates level-1 data to level-2, calls the specified imputation method, 
# and maps results back to original rows

function _imputationLevel2!(
    yData::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    types::Vector{Int},
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String},
    imputationMethod::Function;
    ridge::Float64 = 1e-4,
    kwargs...
    )

    classCols = findall(types .== -2)
    if isempty(classCols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end

    # Extract class information
    classMatrix = Matrix{Float64}(X[:, classCols])
    classKeys = [Tuple(classMatrix[r, c] for c in axes(classMatrix, 2)) for r in axes(classMatrix, 1)]
    classLevels = unique(classKeys)
    nClasses = length(classLevels)
    classMap = Dict(level => idx for (idx, level) in enumerate(classLevels))
    gfFull = [classMap[key] for key in classKeys]

    # Check for partial missing level-2 data
    # (level-2 data should be constant within each class)
    for class in 1:nClasses
        classIdx = findall(gfFull .== class)
        obsIdx = classIdx[.!whereY[classIdx]]
        misIdx = classIdx[whereY[classIdx]]
        
        if !isempty(obsIdx) && !isempty(misIdx)
            clusterIds = join(classLevels[class], ", ")
            throw(ArgumentError("Two-level imputation found partially missing level-2 data in cluster $clusterIds. Use 2lonly.mean or 2lonly.mode to fix inconsistencies."))
        end
    end

    # Aggregate level-1 predictors to level-2 by class means
    randomCols = findall(types .== 2)
    X2lAgg = Matrix{Float64}(undef, nClasses, length(randomCols))
    for i in 1:nClasses
        classIdx = findall(gfFull .== i)
        if !isempty(classIdx)
            X2lAgg[i, :] = vec(mean(X[classIdx, randomCols], dims=1))
        end
    end

    # Create level-2 missing indicator
    # whereY2l[i] = true if all values in class i are missing
    whereY2l = Vector{Bool}(undef, nClasses)
    for i in 1:nClasses
        whereY2l[i] = all(whereY[findall(gfFull .== i)])
    end
    
    # Get observed y values at level-2.
    # Numeric data are averaged; non-numeric data use the class mode.
    yType = nonmissingtype(eltype(yData))
    isNumeric = yType <: Real
    y2l = isNumeric ? Vector{Float64}(undef, nClasses) : Vector{yType}(undef, nClasses)
    for i in 1:nClasses
        classIdx = findall(gfFull .== i)
        obsIdx = classIdx[.!whereY[classIdx]]
        if !isempty(obsIdx)
            classValues = yData[obsIdx]
            y2l[i] = isNumeric ? Float64(mean(classValues)) : mode(classValues)
        end
    end

    # Call the specified imputation method at the aggregated level-2
    imps2l = imputationMethod(y2l, X2lAgg, whereY2l, sum(whereY2l), yVar, iterCounter, j, loggedEvents; ridge=ridge, kwargs...)

    # Map the missing class ids to positions in imps2l.
    missingClassPositions = findall(whereY2l)
    missingClassMap = Dict(missingClass => idx for (idx, missingClass) in enumerate(missingClassPositions))

    # Map level-2 imputations back to original missing rows
    gfₘ = gfFull[whereY]
    return [imps2l[missingClassMap[gfₘ[i]]] for i in 1:whereCount]
end