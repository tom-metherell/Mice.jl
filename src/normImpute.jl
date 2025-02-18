function normImpute!(
    yₒ::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    ridge::Float64 = 1e-5,
    unusedKwargs...
    )

    Xₒ = Matrix{Float64}(hcat(ones(length(whereY) - whereCount), X[.!whereY, :]))
    Xₘ = Matrix{Float64}(hcat(ones(whereCount), X[whereY, :]))

    β̂, β̇, σ̇ = blrDraw!(yₒ, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    return Xₘ * β̇ + randn(whereCount) .* σ̇
end

function blrDraw!(
    yₒ::Vector{Float64},
    Xₒ::Matrix{Float64}, 
    κ::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    β̂ = Xₒ \ yₒ 
    R = qr(Xₒ).R

    V = try
        inv(transpose(R) * R)
    catch
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: ridge penalty applied (unstable results) - predictors are highly multicollinear.")
        S = transpose(R) * R;
        inv(S + Diagonal(S) * κ)
    end

    σ̇ = sqrt(sum((yₒ - Xₒ * β̂).^2)) / rand(Chisq(max(length(yₒ) - size(Xₒ, 2), 1)))
    β̇ = β̂ + σ̇ * cholesky(Hermitian(V)).L * randn(size(Xₒ, 2))

    return β̂, β̇, σ̇
end
