function sample_impute(
    ydata::AbstractArray,
    where_y::Vector{Bool},
    wherecount::Int
    )

    yₒ = ydata[.!where_y]

    # If there are at least some non-missing data
    if length(yₒ) > 0
        imputeddata = sample(yₒ, wherecount)
    elseif yₒ isa CategoricalArray
        # Sample from the levels of the categorical variable
        imputeddata = CategoricalArray{nonmissingtype(eltype(yₒ))}(sample(levels(yₒ), wherecount))
    else
        # Sample from a standard normal distribution
        imputeddata = randn(wherecount)
    end

    return imputeddata
end

const SAMPLE_IMPUTER = Imputer((ydata, X, where_y, wherecount, yvar, itercounter, j, loggedevents; kwargs...) -> begin
    sample_impute(ydata, where_y, wherecount)
end; requirespredictors = false)

registerimputer!("sample", SAMPLE_IMPUTER)