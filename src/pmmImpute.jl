# The pmmImpute! function includes a ! as it updates loggedEvents in place
function pmmImpute!(
    yₒ::AbstractArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    donors::Int = 5,
    ridge::Float64 = 1e-5,
    unusedKwargs...
    )

    # Get the X-values for the rows with observed and missing y-values, respectively
    Xₒ = Matrix{Float64}(hcat(repeat([1], length(whereY) - whereCount), X[.!whereY, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], whereCount), X[whereY, :]))

    # If y is categorical
    if nonmissingtype(eltype(yₒ)) <: Union{AbstractString, CategoricalValue}
        # Convert to dummy variables (as floats) via CCA
        mapping = Dict(levels(yₒ)[i] => i-1 for i in eachindex(levels(yₒ)))
        yNum = Vector{Float64}([mapping[v] for v in yₒ])
        yNum = quantify(yNum, Xₒ)
    else
        yNum = Vector{Float64}(yₒ)
    end

    # Draw from Bayesian linear regression
    β̂, β̇, σ̇ = blrDraw!(yNum, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    # Calculate predicted y-values (for type-1 matching)
    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    # Match predicted y-values with donors
    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return yₒ[indices]    
end

# Comments are as above
function pmmImpute!(
    yₒ::CategoricalArray,
    X::Matrix{Float64},
    whereY::Vector{Bool},
    whereCount::Int,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String};
    donors::Int = 5,
    ridge::Float64 = 1e-5,
    unusedKwargs...
    )

    Xₒ = Matrix{Float64}(hcat(repeat([1], sum(.!whereY)), X[.!whereY, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], whereCount), X[whereY, :]))

    yNum = quantify(yₒ, Xₒ)

    β̂, β̇, σ̇ = blrDraw!(yNum, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return yₒ[indices]
end

function matchIndex(
    ŷₒ::Vector{Float64}, 
    ẏₘ::Vector{Float64},
    donors::Int
    )

    # Shuffle records to remove effects of ties
    nₒ = length(ŷₒ)
    ishuf = randperm(nₒ)
    yshuf = ŷₒ[ishuf]

    # Obtain sorting order on shuffled data
    isort = sortperm(yshuf)

    # Calculate index on input data and sort
    id = ishuf[isort]
    ysort = ŷₒ[id]

    # Pre-sample n_m values between 1 and the number of donors
    nₘ = length(ẏₘ)
    donors = min(donors, nₒ)
    donors = max(donors, 1)
    selections = sample(1:donors, nₘ, replace = true)

    indices = similar(ẏₘ, Int)

    # Loop over the target units
    for i in eachindex(ẏₘ)
        value = ẏₘ[i]
        donorID = selections[i]
        count = 0

        # Find the two adjacent neighbours
        r = searchsortedfirst(ysort, value)
        l = r - 1

        # Find the donorID'th nearest neighbour
        # Store the index of that neighbour
        while count < donorID && l >= 1 && r <= nₒ
            if value - ysort[l] < ysort[r] - value
                indices[i] = id[l]
                l -= 1
            else
                indices[i] = id[r]
                r += 1
            end
            count += 1
        end

        # If right side is exhausted, take left elements
        while count < donorID && l >= 1
            indices[i] = id[l]
            l -= 1
            count += 1
        end

        # If left side is exhausted, take right elements
        while count < donorID && r <= nₒ
            indices[i] = id[r]
            r += 1
            count += 1
        end
    end

    return indices
end
