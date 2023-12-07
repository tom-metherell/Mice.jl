function sampleImpute!(
    yₒ::AbstractArray,
    whereCount::Int
    )

    # If there are at least some non-missing data
    if length(yₒ) > 0
        imputedData = sample(yₒ, whereCount)
    elseif y isa CategoricalArray
        # Sample from the levels of the categorical variable
        imputedData = CategoricalArray{nonmissingtype(eltype(yₒ))}(sample(levels(yₒ), whereCount))
    else
        # Sample from a standard normal distribution
        imputedData = randn(whereCount)
    end

    return imputedData
end