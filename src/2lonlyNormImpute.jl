# The secondLevelOnlyNormImpute! function includes a ! as it updates loggedEvents in place
function secondLevelOnlyNormImpute!(
    yData::AbstractArray,
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
        return Vector{eltype(yData)}(undef, 0)
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

    yₒRaw = yData[.!whereY]
    yₒ = Vector{Float64}(yₒRaw)

    if length(unique(gf)) < nClasses
        throw(ArgumentError("Two-level imputation requires at least one observed outcome per class."))
    end

    # Calculate mean and variance for each group
    groupMeans = Float64[mean(yₒ[gf .== class]) for class in 1:nClasses]
    groupVars = Float64[var(yₒ[gf .== class]) for class in 1:nClasses]

    # Handle zero variance groups
    groupVars = [v > 0 ? v : 1e-6 for v in groupVars]

    # Get class assignments for missing observations
    gfMissing = gfFull[whereY]

    # Impute from normal distribution for each group
    imputedValues = [randn() * sqrt(groupVars[gfMissing[i]]) + groupMeans[gfMissing[i]] for i in 1:whereCount]

    return imputedValues
end

const SECOND_LEVEL_ONLY_NORM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge::Float64 = 1e-4, kwargs...) -> begin
    secondLevelOnlyNormImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge = ridge, kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.norm", SECOND_LEVEL_ONLY_NORM_IMPUTER)
