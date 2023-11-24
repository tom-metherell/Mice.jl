function sampleImpute!(
    y::AbstractArray,
    whereY::Vector{Bool},
    whereCount::Int
    )

    # If there are at least some non-missing data
    if whereCount < length(y)
        imputedData = sample(y[.!whereY], whereCount)
    elseif y isa CategoricalArray
        # Sample from the levels of the categorical variable
        imputedData = CategoricalArray{nonmissingtype(eltype(y))}(sample(levels(y), length(y)))
    else
        # Sample from a standard normal distribution
        imputedData = randn(length(y))
    end

    return imputedData
end