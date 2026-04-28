# The secondLevelOnlyMeanImpute! function includes a ! as it updates loggedEvents in place
function secondLevelOnlyMeanImpute!(
    yData::Vector{Float64},
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
        return Float64[]
    end

    classCols = findall(types .== -2)

    if isempty(classCols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end

    classMatrix = Matrix{Float64}(X[:, classCols])
    classKeys = [Tuple(classMatrix[r, c] for c in axes(classMatrix, 2)) for r in axes(classMatrix, 1)]
    classLevels = unique(classKeys)
    nClasses = length(classLevels)
    classMap = Dict(level => idx for (idx, level) in enumerate(classLevels))
    gfFull = [classMap[key] for key in classKeys]
    gf = gfFull[.!whereY]

    yₒ = yData[.!whereY]

    if length(unique(gf)) < nClasses
        throw(ArgumentError("Two-level imputation requires at least one observed outcome per class."))
    end

    # Get class assignments for missing observations
    gfₘ = gfFull[whereY]

    # Calculate class means
    classStats = [mean(yₒ[gf .== class]) for class in 1:nClasses]

    # Impute with class means
    return [classStats[gfₘ[i]] for i in 1:whereCount]
end

const SECOND_LEVEL_ONLY_MEAN_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...) -> begin
    secondLevelOnlyMeanImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.mean", SECOND_LEVEL_ONLY_MEAN_IMPUTER)
