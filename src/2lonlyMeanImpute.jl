# The secondLevelOnlyMeanImpute! function includes a ! as it updates loggedEvents in place
function secondLevelOnlyMeanImpute!(
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

    # Calculate mean for each group
    groupMeans = [mean(yₒ[gf .== class]) for class in 1:nClasses]

    # Get class assignments for missing observations
    gfMissing = gfFull[whereY]

    # Impute with group means
    imputedValues = [groupMeans[gfMissing[i]] for i in 1:whereCount]

    return eltype(yData) <: Union{AbstractString, CategoricalValue} ? yₒRaw[[argmin(abs.(yₒ .- val)) for val in imputedValues]] : imputedValues
end

const SECOND_LEVEL_ONLY_MEAN_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...) -> begin
    secondLevelOnlyMeanImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.mean", SECOND_LEVEL_ONLY_MEAN_IMPUTER)
