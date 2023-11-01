struct Mipo
    coeftable::CoefTable
    coefnames::Vector{String}
    coefs::Vector{Float64}
    stderrors::Vector{Float64}
    tvalues::Vector{Float64}
    pvalues::Vector{PValue}
end

function pool(
    mira::Mira
    )

    coefs = transpose(reduce(hcat, coef.(mira.analyses)))
    stderrors = transpose(reduce(hcat, stderror.(mira.analyses)))

    pooledCoefs = mean.(eachcol(coefs))
    V_W = mean.(eachcol(stderrors .^ 2))
    V_B = var.(eachcol(coefs))
    V_T = V_W + V_B + V_B/length(mira.analyses)

    pooledStderrors = sqrt.(V_T)
    λ = (V_B .+ V_B ./ length(mira.analyses)) ./ V_T 
    df_Old = (length(mira.analyses) - 1) ./ λ.^2
    n = nobs(mira.analyses[1])
    k = length(pooledCoefs)
    df_Observed = (n - k + 1)/(n - k + 3) * (n - k) .* (1 .- λ)
    df_Adjusted = (df_Old .* df_Observed) ./ (df_Old .+ df_Observed)
    tvalues = pooledCoefs ./ pooledStderrors
    pvalues = PValue.(1 .- cdf.(TDist.(df_Adjusted), abs.(tvalues)))

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