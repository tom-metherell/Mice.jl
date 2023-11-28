# The sampler! function includes a ! as it updates imputations, meanTraces, varTraces and loggedEvents in place
function sampler!(
    imputations::Vector{Matrix},
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    data::T,
    imputeWhere::AxisVector{Vector{Bool}},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisVector{String},
    predictorMatrix::AxisMatrix{Bool},
    iter::Int,
    iterCounter::Int,
    i::Int,
    progressReports::Bool,
    loggedEvents::Vector{String},
    threads::Bool;
    kwargs...
    ) where {T}
    istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))
    
    # Grab name of variable to be imputed
    yVar = visitSequence[i]

    # Grab the variable data
    y = deepcopy(getcolumn(data, Symbol(yVar)))

    # Grab locations of data to be imputed, and set these values to missing (in case of over-imputation)
    whereY = imputeWhere[yVar]
    whereCount = sum(whereY)
    
    if whereCount > 0
        if !(Missing <: eltype(y))
            y = similar(y, Union{Missing, eltype(y)}) .= y
        end
        y[whereY] .= missing
    end

    # Grab the variable's predictors
    predictorVector = predictorMatrix[yVar, :]

    # Grab the names of the predictors
    predictors = names(predictorVector)[1][predictorVector]

    # If there is at least one predictor
    if length(predictors) > 0
        # For variables using an unconditional method
        # Also check that there are actually data to be imputed in the column
        if methods[yVar] ∈ ["mean", "sample"] && any(whereY)
            # For each imputation
            for j in 1:m
                # Impute the data by the specified method
                if methods[yVar] == "mean"
                    # Impute the missing data with the mean of the observed data
                    imputedData = repeat([mean(y[.!whereY])], whereCount)
                elseif methods[yVar] == "sample"
                    # Impute the missing data by sampling from the observed data
                    imputedData = sampleImpute!(y, whereY, whereCount)
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

        # For variables using a conditional method
        # Also check that there are actually data to be imputed in the column
        elseif methods[yVar] ∈ ["norm", "pmm"] && any(whereY)

            # For multithreaded execution
            if threads
                # Create a lock to prevent data race conditions
                lk = ReentrantLock()

                # For each imputation
                Threads.@threads for j in 1:m
                    # Initialise a temporary events log
                    tempLog = Vector{String}([])

                    # Grab the predictors' data
                    X = deepcopy(columntable(NamedTuple{Tuple(Symbol.(predictors))}(Tuple(columntable(data)[c] for c in Symbol.(predictors)))))

                    # Fill missings in the predictors with their imputed values
                    fillXMissings!(X, imputeWhere, predictors, visitSequence, imputations, j)

                    # Convert categorical variables to dummy equivalents
                    # NB: "pacifier" = "dummy" in British English
                    X = pacify!(X, predictors, tempLog, iterCounter, yVar, j)

                    # Store the original number of (dummy) columns
                    origNCol = size(X, 2)

                    # Remove linear dependencies
                    removeLinDeps!(X, y, whereY, whereCount)

                    # If there are still some predictors
                    if size(X, 2) > 0
                        # If some (dummy) predictors were removed
                        if size(X, 2) < origNCol
                            # Calculate the difference
                            diff = OrigNCol - size(X, 2)

                            # Log an event explaining why the predictors where removed
                            push!(tempLog, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end

                        # Impute the missing data by the specified method
                        if methods[yVar] == "pmm"
                            imputedData = pmmImpute!(y, X, whereY, whereCount, 5, 1e-5, yVar, iterCounter, j, tempLog)
                        elseif methods[yVar] == "norm"
                            imputedData = normImpute!(y, X, whereY, whereCount, 1e-5, yVar, iterCounter, j, tempLog)
                        end
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
                    X = deepcopy(columntable(NamedTuple{Tuple(Symbol.(predictors))}(Tuple(columntable(data)[c] for c in Symbol.(predictors)))))
                    fillXMissings!(X, imputeWhere, predictors, visitSequence, imputations, j)
                    X = pacify!(X, predictors, loggedEvents, iterCounter, yVar, j)
                    origNCol = size(X, 2)
                    removeLinDeps!(X, y, whereY, whereCount)

                    if size(X, 2) > 0
                        if size(X, 2) < origNCol
                            diff = origNCol - size(X, 2)
                            push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                        end
                        if methods[yVar] == "pmm"
                            imputedData = pmmImpute!(y, X, whereY, whereCount, 5, 1e-5, yVar, iterCounter, j, loggedEvents)
                        elseif methods[yVar] == "norm"
                            imputedData = normImpute!(y, X, whereY, whereCount, 1e-5, yVar, iterCounter, j, loggedEvents)
                        end
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
            # Invalid method specified
            elseif !(methods[yVar] ∈ ["mean", "norm", "pmm", "sample"])
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
    if imputedData isa CategoricalArray || nonmissingtype(eltype(imputedData)) <: Union{AbstractString, CategoricalValue}
        # Convert the imputed data to integers
        mapping = Dict(levels(imputedData)[i] => i-1 for i in eachindex(levels(imputedData)))
        imputedData = [mapping[v] for v in imputedData]
    end

    # Find the mean and variance and append these to the traces
    meanTraces[i][iterCounter, j] = mean(imputedData)
    varTraces[i][iterCounter, j] = var(imputedData)
end