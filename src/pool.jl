"""
    Mipo

A type for storing the pooled results of multiply imputed repeated analyses (`Mira`).
"""
struct Mipo
    coeftable::CoefTable
    coefnames::Vector{String}
    coefs::Vector{Float64}
    stderrors::Vector{Float64}
    tvalues::Vector{Float64}
    pvalues::Vector{PValue}
end

"""
    pool(mira::Mira)

Pools the results of multiply imputed repeated analyses (`Mira`).
The function will work on any `Mira` object containing model outputs which are
receptive to the `coef`, `stderror` and `nobs` functions from StatsAPI.jl.
"""
function pool(mira::Mira)
    # Grab coefficients and standard errors from each analysis
    coefs = transpose(reduce(hcat, coef.(mira.analyses)))
    stderrors = transpose(reduce(hcat, stderror.(mira.analyses)))

    # Calculate pooled coefficients and standard errors
    pooledcoefs = mean.(eachcol(coefs))
    V_W = mean.(eachcol(stderrors .^ 2))
    V_B = var.(eachcol(coefs))
    V_T = V_W + V_B + V_B/length(mira.analyses)
    pooledstderrors = sqrt.(V_T)

    # Calculate degrees of freedom, t-values and p-values
    λ = (V_B .+ V_B ./ length(mira.analyses)) ./ V_T 
    df_old = (length(mira.analyses) - 1) ./ λ.^2
    n = nobs(mira.analyses[1])
    k = length(pooledcoefs)
    df_observed = (n - k + 1)/(n - k + 3) * (n - k) .* (1 .- λ)
    df_adjusted = (df_old .* df_observed) ./ (df_old .+ df_observed)
    tvalues = pooledcoefs ./ pooledstderrors
    pvalues = PValue.(ccdf.(FDist.(1, df_adjusted), abs2.(tvalues)))

    # Producing tidy table of coefficients (to mirror outputs from StatsModels.jl)
    pooledcoefficients = CoefTable(
        [
            pooledcoefs,
            pooledstderrors,
            tvalues,
            df_adjusted,
            pvalues,
            pooledcoefs .+ quantile.(TDist.(df_adjusted), 0.025) .* pooledstderrors,
            pooledcoefs .+ quantile.(TDist.(df_adjusted), 0.975) .* pooledstderrors
        ],
        ["Coef.", "Std. Error", "t", "df", "Pr(>|t|)", "Lower 95%", "Upper 95%"],
        coefnames(mira.analyses[1]),
        5,
        3
    )

    return Mipo(
        pooledcoefficients,
        coefnames(mira.analyses[1]),
        pooledcoefs,
        pooledstderrors,
        tvalues,
        pvalues
    )
end