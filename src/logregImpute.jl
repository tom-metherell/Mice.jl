# The logregImpute! function includes a ! as it updates loggedEvents in place
function logregImpute!(
    y::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    unusedKwargs...
    )

    # Augment data to avoid perfect prediction (à la White et al., 2010)
    if sum(.!whereY) > 1
        XAug, yAug, whereYAug, weights = augment(y, X, whereY)
    end

    Xₒ = hcat(ones(length(whereYAug) - whereCount), XAug[.!whereYAug, :])
    Xₘ = hcat(ones(whereCount), XAug[whereYAug, :])
    yₒ = yAug[.!whereYAug]

    modelFit = glm(Xₒ, yₒ, Binomial(), LogitLink(), wts = weights[.!whereYAug])

    β̂ = coef(modelFit)
    V = cholesky(Hermitian(inv(cholesky(modelFit.pp)))).L
    β̇ = β̂ + V * randn(size(V, 2))

    p = 1 ./ (1 .+ exp.(-Xₘ * β̇))
    return rand.(Bernoulli.(p))
end

function augment(
    y::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool}
    )

    cats = sort(unique(skipmissing(y)))
    numCats = length(cats)
    ncolX = size(X, 2)

    means = mean.(skipmissing(eachcol(X)))
    sds = std.(skipmissing(eachcol(X)))
    minX = minimum.(skipmissing(eachcol(X)))
    maxX = maximum.(skipmissing(eachcol(X)))

    nrow = 2 * ncolX * numCats

    A = repeat(means', nrow)
    B = reshape(vcat(repeat(vcat(repeat([0.5, -0.5], numCats), repeat([0], nrow)), ncolX - 1), repeat([0.5, -0.5], numCats)), nrow, ncolX)
    C = repeat(sds', nrow)
    D = A + B .* C
    D = [maximum(skipmissing([minX[j], D[i, j]])) for i ∈ Base.axes(D, 1), j ∈ Base.axes(D, 2)]
    D = [minimum(skipmissing([maxX[j], D[i, j]])) for i ∈ Base.axes(D, 1), j ∈ Base.axes(D, 2)]
    e = repeat(repeat(cats, inner = 2), ncolX)

    XAug = vcat(X, D)
    yAug = vcat(y, e)
    whereYAug = vcat(whereY, fill(false, nrow))
    weights = vcat(ones(length(y)), fill((ncolX + 1) / nrow, nrow))

    return XAug, yAug, whereYAug, weights
end
