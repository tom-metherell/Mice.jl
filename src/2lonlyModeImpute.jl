# The secondLevelOnlyModeImpute! function for categorical data
function secondLevelOnlyModeImpute!(
    yData::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    types::Vector{Int},
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    unusedKwargs...
    )

    if whereCount == 0
        return Vector{eltype(yData)}(undef, 0)
    end

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

    # Calculate mode for each class
    classModes = Vector{eltype(yData)}(undef, nClasses)
    for i in 1:nClasses
        classIdx = findall(gfFull .== i)
        obsIdx = classIdx[.!whereY[classIdx]]
        if !isempty(obsIdx)
            classValues = yData[obsIdx]
            classModes[i] = mode(classValues)
        end
    end

    # Map modes back to original missing rows
    gfₘ = gfFull[whereY]
    return [classModes[gfₘ[i]] for i in 1:whereCount]
end

const SECOND_LEVEL_ONLY_MODE_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...) -> begin
    secondLevelOnlyModeImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.mode", SECOND_LEVEL_ONLY_MODE_IMPUTER)
