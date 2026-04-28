# The twoLevelPmmImpute! function includes a ! as it updates loggedEvents in place
function twoLevelPmmImpute!(
    yData::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    types::Vector{Int},
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    intercept::Bool = true,
    nIterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    noise::Float64 = 1e5,
    unusedKwargs...
    )

    if whereCount == 0
        return Vector{eltype(yData)}(undef, 0)
    end

    prep = prepareTwoLevelImputationInputs(X, whereY, types; intercept = intercept)
    gfFull = prep.gfFull
    gf = prep.gf
    nClasses = prep.nClasses

    yₒRaw = yData[.!whereY]

    if length(unique(gf)) < nClasses
        throw(ArgumentError("Two-level PMM requires at least one observed outcome per class."))
    end

    Xₒ = prep.Xₒ
    Xₘ = prep.Xₘ

    # If y is categorical
    if yₒRaw isa CategoricalArray
        # Convert to dummy variables (as floats) via CCA
        yNum = quantify(yₒRaw, Xₒ)
    elseif nonmissingtype(eltype(yₒRaw)) <: Union{AbstractString, CategoricalValue}
        mapping = Dict(levels(yₒRaw)[i] => i-1 for i in eachindex(levels(yₒRaw)))
        yNum = Vector{Float64}([mapping[v] for v in yₒRaw])
        yNum = quantify(yNum, Xₒ)
    else
        yNum = Vector{Float64}(yₒRaw)
    end

    XG = [Matrix{Float64}(Xₒ[gf .== class, :]) for class in 1:nClasses]
    Xₛₛ = [transpose(A) * A for A in XG]
    yG = [yNum[gf .== class] for class in 1:nClasses]
    nG = [sum(gf .== class) for class in 1:nClasses]

    β, invσ² = drawTwoLevelParams(XG, Xₛₛ, yG, nG, nIterations, ridge)

    gfₘ = gfFull[whereY]

    # Predicted values for missing units based on sampled parameters.
    ŷₘ = [dot(Xₘ[r, :], β[gfₘ[r], :]) for r in eachindex(gfₘ)]

    βList = [vec(β[c, :]) for c in 1:nClasses]
    ŷₒ = [dot(Xₒ[r, :], βList[gf[r]]) for r in eachindex(gf)]

    indices = matchIndex(ŷₒ, ŷₘ, donors)
    return yₒRaw[indices]
end

const TWO_LEVEL_PMM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents;
    intercept::Bool = true,
    nIterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    matchSampledPars::Bool = true,
    noise::Float64 = 1e5,
    kwargs...) -> begin
    twoLevelPmmImpute!(
        yData,
        X,
        whereY,
        whereCount,
        types,
        yVar,
        iterCounter,
        j,
        loggedEvents;
        intercept = intercept,
        nIterations = nIterations,
        ridge = ridge,
        donors = donors,
        matchSampledPars = matchSampledPars,
        noise = noise,
        kwargs...
    )
end; twoLevel = true)

registerImputer!("2l.pmm", TWO_LEVEL_PMM_IMPUTER)
