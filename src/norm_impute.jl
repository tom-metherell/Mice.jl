function norm_impute!(
    ydata::Vector{Float64},
    X::Matrix{Float64},
    where_y::Vector{Bool},
    wherecount::Int,
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String};
    ridge::Float64 = 1e-5,
    unusedKwargs...
    )

    yₒ = ydata[.!where_y]

    Xₒ = Matrix{Float64}(hcat(repeat([1], length(where_y) - wherecount), X[.!where_y, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], wherecount), X[where_y, :]))

    β̂, β̇, σ̇ = blrdraw!(yₒ, Xₒ, ridge, yvar, itercounter, j, loggedevents)

    return Xₘ * β̇ + randn(wherecount) .* σ̇
end

function blrdraw!(
    yₒ::Vector{Float64},
    Xₒ::Matrix{Float64}, 
    κ::Float64,
    yvar::String,
    itercounter::Int,
    j::Int,
    loggedevents::Vector{String}
    )

    β̂ = Xₒ \ yₒ 
    R = qr(Xₒ).R

    V = try
        inv(transpose(R) * R)
    catch
        push!(loggedevents, "Iteration $itercounter, variable $yvar, imputation $j: ridge penalty applied (unstable results) - predictors are highly multicollinear.")
        S = transpose(R) * R;
        inv(S + Diagonal(S) * κ)
    end

    σ̇ = sqrt(sum((yₒ - Xₒ * β̂).^2)) / rand(Chisq(max(length(yₒ) - size(Xₒ, 2), 1)))
    β̇ = β̂ + σ̇ * cholesky(Hermitian(V)).L * randn(size(Xₒ, 2))

    return β̂, β̇, σ̇
end

const NORM_IMPUTER = Imputer((ydata, X, where_y, wherecount, yvar, itercounter, j, loggedevents; ridge::Float64 = 1e-5, kwargs...) -> begin
    norm_impute!(ydata, X, where_y, wherecount, yvar, itercounter, j, loggedevents; ridge = ridge, kwargs...)
end)

registerimputer!("norm", NORM_IMPUTER)
