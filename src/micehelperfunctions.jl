"""
    makeMonotoneSequence(data::DataFrame)

Returns a vector of the column names in a DataFrame in ascending order of missingness.
This facilitates convergence in cases where missingness follows a "monotone" pattern.
It is the default visit sequence for the `mice()` function.
"""
function makeMonotoneSequence(data::DataFrame)
    # Initialise missingness vector
    missingness = Vector{Int}(undef, size(data, 2))

    # Count missing data in each column
    for i in axes(data, 2)
        missingness[i] = sum(ismissing.(data[:, i]))
    end

    # Sort the missingness vector in ascending order
    missingness = sortperm(missingness)

    # Sort the data frame names vector by missingness
    visitSequence = names(data)[missingness]

    return visitSequence
end

"""
    makeMethods(data::DataFrame)

Returns a named vector of strings defining the method by which each variable in `data`
should be imputed in the `mice()` function. The default (and only supported) method is
predictive mean matching (pmm).
"""
function makeMethods(data::DataFrame)
    # Use pmm for all variables by default
    methods = NamedArray(Vector{String}(fill("pmm", ncol(data))))

    # Grab the names of the variables
    setnames!(methods, names(data), 1)

    return methods
end

"""
    makePredictorMatrix(data::DataFrame)

Returns a named matrix of booleans defining the predictors for each variable in `data`.
The variables to be predicted are on the rows, and the predictors are on the columns.
The default is to use all variables as predictors for all other variables (i.e. all
1s except for the diagonal, which is 0).
"""
function makePredictorMatrix(data::DataFrame)
    # Initialise the predictor matrix with 1s
    predictorMatrix = NamedArray(Matrix{Bool}(fill(1, ncol(data), ncol(data))))
    
    # Set the diagonal to 0
    for i in 1:ncol(data)
        predictorMatrix[i, i] = 0
    end

    # Grab the names of the variables
    setnames!(predictorMatrix, names(data), 1)
    setnames!(predictorMatrix, names(data), 2)

    return predictorMatrix
end

function initialiseImputations(
    data::DataFrame,
    m::Int,
    visitSequence::Vector{String},
    methods::NamedVector{String}
    )

    # Initialise vector of imputations matrices
    imputations = Vector{Matrix}(undef, ncol(data))

    # For each variable
    for i in eachindex(visitSequence)

        # Grab the variable name
        yVar = visitSequence[i]

        # If the variable is to be imputed
        if methods[yVar] != ""

            # Grab the variable data
            y = data[:, yVar]
            
            # Get locations of missing data in column
            yMissings = ismissing.(data[:, yVar])

            # Count missing data in column
            missingDataCount = sum(yMissings)

            # Initialise imputation matrix
            imputations[i] = Matrix{nonmissingtype(eltype(y))}(undef, missingDataCount, m)

            # If there are at least some non-missing data
            if !all(yMissings)
                # For each imputation
                for j in 1:m
                    # Sample observed data at random to serve as initial values
                    imputations[i][:, j] = sample(y[.!yMissings], missingDataCount)
                end
            else
                if y isa CategoricalArray
                    # For each imputation
                    for j in 1:m
                        # Sample from the levels of the categorical variable
                        imputations[i][:, j] = CategoricalArray{nonmissingtype(eltype(y))}(sample(levels(y), length(y)))
                    end
                else
                    # For each imputation
                    for j in 1:m
                        # Sample from a standard normal distribution
                        imputations[i][:, j] .= randn(length(y))
                    end
                end
            end
        end
    end

    return imputations
end

# Allow US spelling of initialise
const initializeImputations = initialiseImputations

function initialiseTraces(
    visitSequence::Vector{String},
    iter::Int,
    m::Int
    )

    traces = [Matrix{Float64}(undef, iter, m) for _ = eachindex(visitSequence)]

    return traces
end

# Allow US spelling of initialise
const initializeTraces = initialiseTraces

# The sampler! function includes a ! as it updates imputations, meanTraces, varTraces and loggedEvents in place
function sampler!(
    imputations::Vector{Matrix},
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    data::DataFrame,
    m::Int,
    visitSequence::Vector{String},
    methods::NamedVector{String},
    predictorMatrix::NamedArray{Bool},
    iter::Int,
    iterCounter::Int,
    i::Int,
    progressReports::Bool,
    loggedEvents::Vector{String},
    threads::Bool;
    kwargs...
    )
    
    # Grab name of variable to be imputed
    yVar = visitSequence[i]

    # Grab the variable data
    y = data[:, yVar]

    # Grab the variable's predictors
    predictorVector = predictorMatrix[yVar, :]

    # Grab the names of the predictors
    predictors = names(predictorVector)[1][predictorVector]

    # If there is at least one predictor
    if length(predictors) > 0

        # For variables using pmm (in v0.0.0 this is the only supported method)
        # Also check that there is actually missing data in the column
        if methods[yVar] == "pmm" && any(ismissing.(y))

            # For multithreaded execution
            if threads
                # Create a lock to prevent data race conditions
                lk = ReentrantLock()

                # For each imputation
                Threads.@threads for j in 1:m
                    # Initialise a temporary events log
                    tempLog = Vector{String}([])

                    # Grab the predictors' data
                    X = data[:, predictors]

                    # Fill missings in the predictors with their imputed values
                    fillXMissings!(X, predictors, visitSequence, imputations, j)

                    # Convert categorical variables to dummy equivalents
                    # NB: "pacifier" = "dummy" in British English
                    X = pacify!(X, predictors, tempLog, iterCounter, yVar, j)

                    # Store the original number of (dummy) columns
                    origNCol = size(X, 2)

                    # Remove linear dependencies
                    removeLinDeps!(X, y)

                    # If there are still some predictors
                    if size(X, 2) > 0
                        # If some (dummy) predictors were removed
                        if size(X, 2) < origNCol
                            # Calculate the difference
                            diff = OrigNCol - size(X, 2)

                            # Log an event explaining why the predictors where removed
                            push!(tempLog, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end

                        # Impute the missing data by predictive mean matching
                        imputedData = pmmImpute!(y, X, 5, 1e-5, yVar, iterCounter, j, tempLog)
                    else
                        # Log an event explaining why the imputation was skipped
                        push!(tempLog, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                        
                        # Use the imputed data from the previous iteration
                        imputedData = imputations[i][:, j]
                    end

                    # Activate reentrant lock
                    lock(lk)
                    try
                        # Update mean and variance traces
                        updateTraces!(meanTraces, varTraces, imputedData, i, iterCounter, j)
                        
                        # Append imputed data to imputations matrix
                        imputations[i][:, j] = imputedData

                        # Append temporary events log to main log
                        append!(loggedEvents, tempLog)
                    finally
                        unlock(lk)
                    end                    
                end

                # Print progress indicator
                if progressReports
                    progress = ((iterCounter - 1)/iter + (i/length(visitSequence))/iter) * 100
                    progressRound = floor(Int8, progress / 10)
                    miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                    @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar miceEmojis progress length(loggedEvents)
                end
            else
                # Comments are as above
                for j in 1:m
                    X = data[:, predictors]
                    fillXMissings!(X, predictors, visitSequence, imputations, j)
                    X = pacify!(X, predictors, loggedEvents, iterCounter, yVar, j)
                    origNCol = size(X, 2)
                    removeLinDeps!(X, y)

                    if size(X, 2) > 0
                        if size(X, 2) < origNCol
                            diff = origNCol - size(X, 2)
                            push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end
                        imputedData = pmmImpute!(y, X, 5, 1e-5, yVar, iterCounter, j, loggedEvents)
                    else
                        push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                        imputedData = imputations[i][:, j]
                    end

                    updateTraces!(meanTraces, varTraces, imputedData, i, iterCounter, j)

                    imputations[i][:, j] = imputedData

                    if progressReports
                        progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                        progressRound = floor(Int8, progress / 10)
                        miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                        @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress length(loggedEvents)
                    end
                end
            end
        else
            # No method specified
            if methods[yVar] == ""
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no method specified.")
            # Method other than pmm specified
            elseif methods[yVar] != "pmm"
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - method not supported.")
            # Neither of these => there is no missing data
            else
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no missing data.")
            end
        end
    # No predictors remain
    else
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no predictors.")
    end
end

# The fillXMissings! function includes a ! as it updates X in place
function fillXMissings!(
    X::DataFrame,
    predictors::Vector{String},
    visitSequence::Vector{String},
    imputations::Vector{Matrix},
    j::Int
    )

    # For each predictor
    for k in predictors
        # Find its position in the visit sequence
        kVS = findfirst(visitSequence .== k)

        # Find the positions of missing data in the predictor column
        xMissings = ismissing.(X[:, k])

        # If there are any missing data
        if any(xMissings)
            # Fill them with the imputed values
            X[xMissings, k] = imputations[kVS][:, j]
        end
    end
end

# The pacify! function includes a ! as it updates loggedEvents in place
function pacify!(
    X::DataFrame,
    predictors::Vector{String},
    loggedEvents::Vector{String},
    iterCounter::Int,
    yVar::AbstractString,
    j::Int
    )

    # Initialise vector of categorical predictors
    categoricalPredictors = Vector{String}([])

    # For each predictor
    for xVar in predictors
        # Grab the predictor data
        x = X[:, xVar]

        # If the data are categorical (either CategoricalArray or a vector of strings)
        if x isa CategoricalArray || nonmissingtype(eltype(x)) <: AbstractString
            # If there is more than one level
            if length(levels(x)) > 1
                # Add this variable to the list of categorical predictors
                push!(categoricalPredictors, xVar)
            else
                # Otherwise, drop this variable
                select!(X, Not(xVar))
                predictors = predictors[predictors .!= xVar]

                # Log that this predictor has been dropped
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: predictor $xVar dropped because of zero variance.")
            end
        # Otherwise, convert the data to float and standardise
        elseif nonmissingtype(eltype(x)) <: Real
            X[!, xVar] = convert.(Float64, X[:, xVar])
            X[:, xVar] = standardize(UnitRangeTransform, X[:, xVar])
        end
    end

    # Define the model frame
    mf = ModelFrame(term(0) ~ sum(term.(predictors)), X)

    # Set contrast coding for categorical predictors to orthogonal polynomial
    setcontrasts!(mf, Dict([Symbol(xVar) => PolynomialCoding() for xVar in categoricalPredictors]))

    # Produce the model matrix
    X = ModelMatrix(mf).m[:, 2:end]

    return X
end

# The PolynomialCoding type is used to specify orthogonal polynomial contrast coding
mutable struct PolynomialCoding <: AbstractContrasts
end

# This extension of the contrasts_matrix function from StatsModels.jl allows orthogonal polynomial contrast coding
function contrasts_matrix(C::PolynomialCoding, _, n)
    X = reduce(hcat, [((1:n) .- mean(1:n)) .^ i for i in 0:n-1])
    qrX = qr(X)
    Z = qrX.Q * Diagonal(qrX.R)
    for i in axes(Z, 2)
        Z[:, i] = Z[:, i] ./ sqrt(sum(Z[:, i].^2))
    end
    return Z[:, 2:end]
end

# This extension of the termnames function from StatsModels.jl names the new dummy variables correctly
function termnames(C::PolynomialCoding, levels::AbstractVector, _::Integer)
    return Vector{String}([".^$i" for i in 1:length(levels)])
end

function pacify(y::AbstractArray)
    # Grab the levels of the variable
    yLevels = levels(y)

    # Initialise the dummy matrix
    yDummies = Matrix{Float64}(undef, length(y), length(yLevels))

    # Convert the variable to dummy variables
    for q in eachindex(yLevels)
        yDummies[:, q] = y .== yLevels[q]
    end

    return yDummies
end

# The removeLinDeps! function includes a ! as it updates X in place
function removeLinDeps!(
    X::Matrix{Float64},
    y::AbstractArray
    )

    # If all y-values are missing, stop now
    if all(ismissing.(y))
        return
    end

    # Grab observed predictor data
    Xₒ = Matrix{Float64}(X[.!ismissing.(y), :])
    
    # If y is categorical
    if y isa CategoricalArray || nonmissingtype(eltype(y)) <: AbstractString
        # Grab observed y-values
        yₒ = y[.!ismissing.(y)]

        # Convert y to dummy variables (as floats)
        mapping = Dict(levels(yₒ)[i] => i-1 for i in eachindex(levels(yₒ)))
        yₒ = Vector{Float64}([mapping[v] for v in yₒ])
    else
        # Grab observed y-values (as floats)
        yₒ = Vector{Float64}(y[.!ismissing.(y)])
    end

    # If the variance of observed y-values falls below the allowed threshold, delete all predictors and stop now
    if var(yₒ) < 1e-4
        select!(X, [])
        return
    end

    # Define vector of predictors to keep as those with sufficiently high variance and sufficiently low correlation with the observed y-values
    keep = var.(eachcol(Xₒ)) .> 1e-4 .&& cor.(eachcol(Xₒ), [yₒ]) .< 0.99

    # Get the total number of predictors currently being kept
    keepSum = sum(keep)

    # If there are 0 or 1 remaining, stop now
    if keepSum < 2
        X = X[:, keep]
        return
    end

    # Otherwise, calculate the correlation matrix of the remaining predictorsa
    xCors = cor(Xₒ)
    nxCors = xCors[findall(keep), findall(keep)]

    # Calculate the eigenvalues and eigenvectors of the correlation matrix and sort them
    eigenCors = eigen(nxCors)
    eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
    sortedeigvals = eigenCors.values[eigvalsorder]
    sortedeigvecs = eigenCors.vectors[:, eigvalsorder]

    # While the largest eigenvalue is more than 10_000 times larger than the smallest 
    while sortedeigvals[keepSum] / sortedeigvals[1] < 1e-4
        # Remove the predictor contributing most to the eigenvector with the smallest eigenvalue
        w = sortperm(abs.(sortedeigvecs[:, keepSum]), rev = true)[1]
        keep[findall(keep)[w]] = false
        nxCors = xCors[findall(keep), findall(keep)]
        keepSum -= 1
        eigenCors = eigen(nxCors)
        eigvalsorder = sortperm(abs.(eigenCors.values), rev = true)
        sortedeigvals = eigenCors.values[eigvalsorder]
        sortedeigvecs = eigenCors.vectors[:, eigvalsorder]
    end

    # Remove the predictors that are not being kept
    X = X[:, keep]
end

# The updateTraces! function includes a ! as it updates meanTraces and varTraces in place
function updateTraces!(
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    imputedData::AbstractArray,
    i::Int,
    iterCounter::Int,
    j::Int
    )

    # If the imputed data are categorical
    if imputedData isa CategoricalArray || nonmissingtype(eltype(imputedData)) <: AbstractString
        # Convert the imputed data to integers
        mapping = Dict(levels(imputedData)[i] => i-1 for i in eachindex(levels(imputedData)))
        imputedData = [mapping[v] for v in imputedData]
    end

    # Find the mean and variance and append these to the traces
    meanTraces[i][iterCounter, j] = mean(imputedData)
    varTraces[i][iterCounter, j] = var(imputedData)
end

# The pmmImpute! function includes a ! as it updates loggedEvents in place
function pmmImpute!(
    y::AbstractArray,
    X::Matrix{Float64},
    donors::Int,
    ridge::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    # Get the positions of missing y-values
    yMissings = ismissing.(y)

    # Get the X-values for the rows with observed and missing y-values, respectively
    Xₒ = Matrix{Float64}(hcat(repeat([1], sum(.!yMissings)), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], sum(yMissings)), X[yMissings, :]))

    # If y is categorical
    if nonmissingtype(eltype(y)) <: AbstractString
        # Convert to dummy variables (as floats) via CCA
        mapping = Dict(levels(y[.!yMissings])[i] => i-1 for i in eachindex(levels(y[.!yMissings])))
        yₒ = Vector{Float64}([mapping[v] for v in y[.!yMissings]])
        yₒ = quantify(y[.!yMissings], Xₒ)
    else
        # Grab observed y-values (as floats)
        yₒ = Vector{Float64}(y[.!yMissings])
    end

    # Draw from Bayesian linear regression
    β̂, β̇ = blrDraw!(yₒ, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    # Calculate predicted y-values (for type-1 matching)
    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    # Match predicted y-values with donors
    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]    
end

# Comments are as above
function pmmImpute!(
    y::CategoricalArray,
    X::Matrix{Float64},
    donors::Int,
    ridge::Float64,
    yVar::String,
    iterCounter::Int,
    j::Int,
    loggedEvents::Vector{String}
    )

    yMissings = ismissing.(y)

    Xₒ = Matrix{Float64}(hcat(repeat([1], sum(.!yMissings)), X[.!yMissings, :]))
    Xₘ = Matrix{Float64}(hcat(repeat([1], sum(yMissings)), X[yMissings, :]))

    yₒ = quantify(y[.!yMissings], Xₒ)

    β̂, β̇ = blrDraw!(yₒ, Xₒ, ridge, yVar, iterCounter, j, loggedEvents)

    ŷₒ = Xₒ * β̂
    ẏₘ = Xₘ * β̇

    indices = matchIndex(ŷₒ, ẏₘ, donors)

    return y[.!yMissings][indices]
end

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
    β̇ = β̂ + σ̇ * cholesky((V + transpose(V)) / 2).factors * randn(size(Xₒ, 2))

    return β̂, β̇
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

export makeMonotoneSequence, makeMethods, makePredictorMatrix