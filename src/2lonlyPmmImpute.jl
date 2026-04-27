# The secondLevelOnlyPmmImpute! function includes a ! as it updates loggedEvents in place
function secondLevelOnlyPmmImpute!(
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
    donors::Int = 5,
    unusedKwargs...
    )

    if whereCount == 0
        return Vector{eltype(yData)}(undef, 0)
    end

    classCols = findall(types .== -2)
    randomCols = findall(types .== 2)

    if isempty(classCols)
        throw(ArgumentError("Two-level PMM method specified, but no class variable (coded -2) found."))
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
        throw(ArgumentError("Two-level PMM requires at least one observed outcome per class."))
    end

    # Get observed and missing predictor values
    Xₒ = Matrix{Float64}(X[.!whereY, randomCols])
    Xₘ = Matrix{Float64}(X[whereY, randomCols])
    
    # Add intercept
    Xₒ = hcat(repeat([1], size(Xₒ, 1)), Xₒ)
    Xₘ = hcat(repeat([1], size(Xₘ, 1)), Xₘ)

    # Convert y to numeric if categorical
    if nonmissingtype(eltype(yₒRaw)) <: Union{AbstractString, CategoricalValue}
        mapping = Dict(levels(yₒRaw)[i] => i-1 for i in eachindex(levels(yₒRaw)))
        yNum = Vector{Float64}([mapping[v] for v in yₒRaw])
        yNum = quantify(yNum, Xₒ)
    else
        yNum = yₒ
    end

    # Perform group-level regression
    gfMissing = gfFull[whereY]
    imputedIndices = Vector{Int}(undef, whereCount)

    for missIdx in 1:whereCount
        classIdx = gfMissing[missIdx]
        # Get data for this class
        classObsMask = gf .== classIdx
        Xₒ_class = Xₒ[classObsMask, :]
        yNum_class = yNum[classObsMask]
        yₒ_class = yₒ[classObsMask]

        # Fit regression on class data
        β̂ = Xₒ_class \ yNum_class

        # Get prediction for missing observation
        ŷₘ = Xₘ[missIdx, :]' * β̂

        # Calculate residuals and predictions for donors
        ŷₒ = Xₒ_class * β̂
        
        # Find nearest donors within the class (PMM)
        distances = abs.(ŷₒ .- ŷₘ)
        nDonors = min(donors, length(distances))
        _, nearestIndices = findmin(distances), sortperm(distances)[1:nDonors]
        
        # Sample from donors
        selectedIdx = rand(nearestIndices)
        # Map back to original indices
        classIndices = findall(classObsMask)
        imputedIndices[missIdx] = classIndices[selectedIdx]
    end

    return yₒRaw[imputedIndices]
end

const SECOND_LEVEL_ONLY_PMM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge::Float64 = 1e-4, donors::Int = 5, kwargs...) -> begin
    secondLevelOnlyPmmImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge = ridge, donors = donors, kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.pmm", SECOND_LEVEL_ONLY_PMM_IMPUTER)
