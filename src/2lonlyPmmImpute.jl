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
        return Float64[]
    end

    return _imputationLevel2!(
        yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents,
        (yₒ, X2l, wy2l, whereCount2l, yVar, iterCounter, j, loggedEvents; ridge=ridge, donors=donors, kwargs...) -> 
            pmmImpute!(yₒ, X2l, wy2l, whereCount2l, yVar, iterCounter, j, loggedEvents; donors=donors, ridge=ridge, kwargs...);
        ridge=ridge,
        donors=donors,
        unusedKwargs...
    )
end

const SECOND_LEVEL_ONLY_PMM_IMPUTER = Imputer((yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge::Float64 = 1e-4, donors::Int = 5, kwargs...) -> begin
    secondLevelOnlyPmmImpute!(yData, X, whereY, whereCount, types, yVar, iterCounter, j, loggedEvents; ridge = ridge, donors = donors, kwargs...)
end; twoLevel = true)

registerImputer!("2lonly.pmm", SECOND_LEVEL_ONLY_PMM_IMPUTER)
