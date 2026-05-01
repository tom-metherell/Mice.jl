function secondlevelonlypmm_impute!(
    ydata::AbstractArray,
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    types::Vector{Int},
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String};
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    unusedkwargs...
    )

    if wherecount == 0
        return Vector{eltype(ydata)}([])
    end

    return _imputationlevel2(
        ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents,
        (yₒ, X2l, wy2l, wherecount2l, yvar, itercounter, j, loggedevents; ridge=ridge, donors=donors, kwargs...) -> 
            pmm_impute!(yₒ, X2l, wy2l, wherecount2l, yvar, itercounter, j, loggedevents; donors=donors, ridge=ridge, kwargs...);
        ridge=ridge,
        donors=donors,
        unusedkwargs...
    )
end

const SECOND_LEVEL_ONLY_PMM_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; ridge::Float64 = 1e-4, donors::Int = 5, kwargs...) -> begin
    secondlevelonlypmm_impute!(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; ridge = ridge, donors = donors, kwargs...)
end; twolevel = true)

registerimputer!("2lonly.pmm", SECOND_LEVEL_ONLY_PMM_IMPUTER)
