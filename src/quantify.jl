function quantify(
    yₒ::AbstractArray,
    Xₒ::Matrix{Float64}
    )

    # Use CCA to convert categorical variables to dummy variables
    yDummies = pacify(yₒ)
    Ycoef = miceCCA(Xₒ, yDummies)
    yₒ = Vector{Float64}(zscore(yDummies * Ycoef[:, 2]))

    return yₒ
end

function miceCCA(
    X::Matrix{Float64},
    Y::Matrix{Float64}
    )

    qrX = qr(X)
    qrY = qr(Y)
    nr = size(X, 1)
    dx = rank(X)
    dy = rank(Y)

    Z = svd((transpose(qrY.Q) * (qrX.Q * diagm(nr, dx, repeat([1], min(nr, dx)))))[1:dy, :])

    Ycoef = qrY.R \ Z.U

    return Ycoef
end