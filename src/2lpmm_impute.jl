function twolevelpmm_impute!(
    ydata::AbstractArray,
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    types::Vector{Int},
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String};
    intercept::Bool = true,
    n_iterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    unusedkwargs...
    )

    if wherecount == 0
        return Vector{eltype(ydata)}(undef, 0)
    end

    prep = preparetwolevelimputationinputs(X, where_y, types; intercept = intercept)
    gf_full = prep.gf_full
    gf = prep.gf
    nclasses = prep.nclasses

    yₒ_raw = ydata[.!where_y]

    if length(unique(gf)) < nclasses
        throw(ArgumentError("Two-level PMM requires at least one observed outcome per class."))
    end

    Xₒ = prep.Xₒ
    Xₘ = prep.Xₘ

    # If y is categorical
    if yₒ_raw isa CategoricalArray
        # Convert to dummy variables (as floats) via CCA
        ynum = quantify(yₒ_raw, Xₒ)
    elseif nonmissingtype(eltype(yₒ_raw)) <: Union{AbstractString, CategoricalValue}
        mapping = Dict(levels(yₒ_raw)[i] => i-1 for i in eachindex(levels(yₒ_raw)))
        ynum = Vector{Float64}([mapping[v] for v in yₒ_raw])
        ynum = quantify(ynum, Xₒ)
    else
        ynum = Vector{Float64}(yₒ_raw)
    end

    XG = [Matrix{Float64}(Xₒ[gf .== class, :]) for class in 1:nclasses]
    Xₛₛ = [transpose(A) * A for A in XG]
    yg = [ynum[gf .== class] for class in 1:nclasses]
    ng = [sum(gf .== class) for class in 1:nclasses]

    β, invσ² = drawtwolevelparams(XG, Xₛₛ, yg, ng, n_iterations, ridge)

    gfₘ = gf_full[where_y]

    # Predicted values for missing units based on sampled parameters.
    ŷₘ = [dot(Xₘ[r, :], β[gfₘ[r], :]) for r in eachindex(gfₘ)]

    βlist = [vec(β[c, :]) for c in 1:nclasses]
    ŷₒ = [dot(Xₒ[r, :], βlist[gf[r]]) for r in eachindex(gf)]

    indices = matchindex(ŷₒ, ŷₘ, donors)
    return yₒ_raw[indices]
end

const TWO_LEVEL_PMM_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents;
    intercept::Bool = true,
    n_iterations::Int = 100,
    ridge::Float64 = 1e-4,
    donors::Int = 5,
    kwargs...) -> begin
    twolevelpmm_impute!(
        ydata,
        X,
        where_y,
        wherecount,
        types,
        yvar,
        itercounter,
        j,
        loggedevents;
        intercept = intercept,
        n_iterations = n_iterations,
        ridge = ridge,
        donors = donors,
        kwargs...
    )
end; twolevel = true)

registerimputer!("2l.pmm", TWO_LEVEL_PMM_IMPUTER)
