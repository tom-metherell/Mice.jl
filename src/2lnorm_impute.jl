# The twolevelnorm_impute! function includes a ! as it updates loggedevents in place
function twolevelnorm_impute!(
    ydata::Vector{Float64},
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
    unusedkwargs...
    )

    prep = preparetwolevelimputationinputs(X, where_y, types; intercept = intercept)
    gf_full = prep.gf_full
    gf = prep.gf
    nclasses = prep.nclasses
    yₒ = Vector{Float64}(ydata[.!where_y])

    if length(unique(gf)) < nclasses
        throw(ArgumentError("Two-level imputation requires at least one observed outcome per class."))
    end

    Xₒ = prep.Xₒ
    Xₘ = prep.Xₘ
    XG = [Matrix{Float64}(Xₒ[gf .== class, :]) for class in 1:nclasses]
    Xₛₛ = [transpose(A) * A for A in XG]
    yg = [yₒ[gf .== class] for class in 1:nclasses]
    ng = [sum(gf .== class) for class in 1:nclasses]

    β, invσ² = drawtwolevelparams(XG, Xₛₛ, yg, ng, n_iterations, ridge)

    if wherecount == 0
        return Vector{Float64}(undef, 0)
    end

    gfₘ = gf_full[where_y]
    imps = Vector{Float64}(undef, wherecount)

    for r in eachindex(gfₘ)
        class = gfₘ[r]
        imps[r] = dot(Xₘ[r, :], β[class, :]) + randn() * sqrt(1 / invσ²[class])
    end

    return imps
end

function preparetwolevelimputationinputs(
    X::Matrix{Float64},
    where_y::Vector{Bool},
    types::Vector{Int};
    intercept::Bool = true
    )

    Xwork = X
    typeswork = types
    if intercept
        Xwork = hcat(ones(Float64, size(Xwork, 1)), Xwork)
        typeswork = vcat(2, typeswork)
    end

    classcols = findall(typeswork .== -2)
    randomcols = findall(typeswork .== 2)
    fixedcols = findall(typeswork .> 0)

    if isempty(classcols)
        throw(ArgumentError("Two-level imputation method specified, but no class variable (coded -2) found."))
    end
    if isempty(randomcols)
        throw(ArgumentError("Two-level imputation requires at least one random predictor (coded 2)."))
    end

    classmatrix = Matrix{Float64}(Xwork[:, classcols])
    classkeys = [Tuple(classmatrix[r, c] for c in axes(classmatrix, 2)) for r in axes(classmatrix, 1)]
    classlevels = unique(classkeys)
    nclasses = length(classlevels)
    classmap = Dict(level => idx for (idx, level) in enumerate(classlevels))
    gf_full = [classmap[key] for key in classkeys]
    gf = gf_full[.!where_y]

    Xₒ = Matrix{Float64}(Xwork[.!where_y, randomcols])
    Xₘ = Matrix{Float64}(Xwork[where_y, randomcols])

    return (
        X = Xwork,
        types = typeswork,
        classcols = classcols,
        randomcols = randomcols,
        fixedcols = fixedcols,
        gf_full = gf_full,
        gf = gf,
        nclasses = nclasses,
        Xₒ = Xₒ,
        Xₘ = Xₘ,
    )
end

function drawtwolevelparams(
    XG::Vector{Matrix{Float64}},
    Xₛₛ::Vector{Matrix{Float64}},
    yg::Vector{Vector{Float64}},
    ng::Vector{Int},
    n_iterations::Int,
    ridge::Float64
    )

    nclasses = length(XG)
    nrc = size(XG[1], 2)

    β = zeros(nclasses, nrc)
    ss = zeros(nclasses)
    μ = zeros(nrc)
    invψ = Matrix{Float64}(Diagonal(ones(nrc)))
    invσ² = ones(nclasses)
    σ²₀ = 1.0
    θ = 1.0

    for _ in 1:n_iterations
        for class in 1:nclasses
            vv = symridge(invσ²[class] * Xₛₛ[class] + invψ, ridge)
            βvar = chol2inv(cholesky(vv))

            meanpart = βvar * (
                transpose(invσ²[class] * XG[class]) * yg[class] + invψ * μ
            )

            noisepart = cholesky(symridge(βvar, ridge)).U' * randn(nrc)
            β[class, :] = meanpart + noisepart

            resid = yg[class] - XG[class] * β[class, :]
            ss[class] = dot(resid, resid)
        end

        μVar = chol2inv(cholesky(symridge(invψ, ridge))) / nclasses
        μ = vec(mean(β, dims = 1)) + cholesky(μVar).U' * randn(nrc)

        centred = β .- transpose(μ)
        ψscale = chol2inv(cholesky(symridge(centred' * centred, ridge)))
        invψ = randwishart(nclasses - nrc - 1, cholesky(ψscale).U)

        for class in 1:nclasses
            shape = ng[class] / 2 + 1 / (2 * θ)
            scale = 2 * θ / (ss[class] * θ + σ²₀)
            invσ²[class] = rand(Gamma(shape, scale))
        end

        H = 1 / mean(invσ²)
        σ²₀ = rand(Gamma(nclasses / (2 * θ) + 1, 2 * θ * H / nclasses))

        G = exp(mean(log.(1 ./ invσ²)))
        θshape = nclasses / 2 - 1
        θscale = 2 / (nclasses * (σ²₀ / H - log(σ²₀) + log(G) - 1))
        θ = 1 / rand(Gamma(θshape, θscale))
    end

    return β, invσ²
end

function symridge(A::AbstractMatrix{Float64}, ridge::Float64 = 1e-4)
    Asym = (A + transpose(A)) ./ 2
    if size(Asym, 1) == 1
        return Matrix{Float64}(Asym)
    end
    return Matrix{Float64}(Asym + Diagonal(diag(Asym) .* ridge))
end

function chol2inv(C)
    return inv(C.U' * C.U)
end

function randwishart(df::Int, sqrtΣ::AbstractMatrix{Float64})
    p = size(sqrtΣ, 1)
    Z = zeros(p, p)

    for i in 1:p
        Z[i, i] = sqrt(rand(Chisq(df - i + 1)))
    end
    for i in 2:p
        for j in 1:i-1
            Z[i, j] = randn()
        end
    end

    B = Z * sqrtΣ
    return transpose(B) * B
end

const TWO_LEVEL_NORM_IMPUTER = Imputer((ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; intercept::Bool = true, n_iterations::Int = 100, ridge::Float64 = 1e-4, kwargs...) -> begin
    twolevelnorm_impute!(ydata, X, where_y, wherecount, types, yvar, itercounter, j, loggedevents; intercept = intercept, n_iterations = n_iterations, ridge = ridge, kwargs...)
end; twolevel = true)

registerimputer!("2l.norm", TWO_LEVEL_NORM_IMPUTER)