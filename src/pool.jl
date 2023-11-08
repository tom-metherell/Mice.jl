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
    pooledCoefs = mean.(eachcol(coefs))
    V_W = mean.(eachcol(stderrors .^ 2))
    V_B = var.(eachcol(coefs))
    V_T = V_W + V_B + V_B/length(mira.analyses)
    pooledStderrors = sqrt.(V_T)

    # Calculate degrees of freedom, t-values and p-values
    λ = (V_B .+ V_B ./ length(mira.analyses)) ./ V_T 
    df_Old = (length(mira.analyses) - 1) ./ λ.^2
    n = nobs(mira.analyses[1])
    k = length(pooledCoefs)
    df_Observed = (n - k + 1)/(n - k + 3) * (n - k) .* (1 .- λ)
    df_Adjusted = (df_Old .* df_Observed) ./ (df_Old .+ df_Observed)
    tvalues = pooledCoefs ./ pooledStderrors
    pvalues = PValue.(1 .- cdf.(TDist.(df_Adjusted), abs.(tvalues)))

    # Producing tidy table of coefficients (to mirror outputs from StatsModels.jl)
    pooledCoefficients = CoefTable(
        [
            pooledCoefs,
            pooledStderrors,
            tvalues,
            df_Adjusted,
            pvalues,
            pooledCoefs .+ quantile.(TDist.(df_Adjusted), 0.025) .* pooledStderrors,
            pooledCoefs .+ quantile.(TDist.(df_Adjusted), 0.975) .* pooledStderrors
        ],
        ["Coef.", "Std. Error", "t", "df", "Pr(>|t|)", "Lower 95%", "Upper 95%"],
        coefnames(mira.analyses[1]),
        5,
        3
    )

    return Mipo(
        pooledCoefficients,
        coefnames(mira.analyses[1]),
        pooledCoefs,
        pooledStderrors,
        tvalues,
        pvalues
    )
end

export Mipo, pool