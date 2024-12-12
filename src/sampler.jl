# The sampler! function includes a ! as it updates workingData, workingDataPacified, meanTraces, varTraces and loggedEvents in place
function sampler!(
    workingData::AxisVector{Vector},
    workingDataPacified::AxisVector,
    workingDataLevels::AxisVector,
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    imputeWhere::AxisVector{Vector{Bool}},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisVector{String},
    predictorMatrix::AxisMatrix{Int},
    iter::Int,
    iterCounter::Int,
    i::Int,
    progressReports::Bool,
    loggedEvents::Vector{String};
    kwargs...
    )
    
    # Grab name of variable to be imputed
    yVar = visitSequence[i]

    # Grab locations of data to be imputed, and set these values to missing (in case of over-imputation)
    whereY = imputeWhere[yVar]
    whereCount = sum(whereY)

    # Grab the names of the predictors
    predictors = axes(predictorMatrix[yVar, :])[1][predictorMatrix[yVar, :] .== 1]

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
                    workingData[yVar][j][whereY] = repeat([mean(workingData[yVar][j][.!whereY])], whereCount)
                elseif methods[yVar] == "sample"
                    # Impute the missing data by sampling from the observed data
                    workingData[yVar][j][whereY] = sampleImpute!(workingData[yVar][j][.!whereY], whereCount)
                end

                updateTraces!(meanTraces, varTraces, workingData[yVar][j][whereY], i, iterCounter, j)

                if workingData[yVar][j] isa CategoricalArray || nonmissingtype(eltype(workingData[yVar][j])) <: Union{AbstractString, CategoricalValue}
                    workingDataPacified[yVar][j] = pacifyWorkingData(workingData[yVar][j], workingDataLevels[yVar])
                end

                if progressReports
                    progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                    progressRound = floor(Int8, progress / 10)
                    miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                    @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress length(loggedEvents)
                end
            end

        # For variables using a conditional method
        # Also check that there are actually data to be imputed in the column
        elseif methods[yVar] ∈ ["norm", "pmm", "rf"] && any(whereY)

            # For each imputation
            for j in 1:m
                # Grab the predictors' data
                X = Matrix{Float64}(reduce(hcat, [predictor ∈ axes(workingDataPacified)[1] ? workingDataPacified[predictor][j] : workingData[predictor][j] for predictor in predictors]))
                    
                # Store the original number of (dummy) columns
                origNCol = size(X, 2)

                # Remove linear dependencies
                removeLinDeps!(X, workingData[yVar][j], whereY, whereCount)

                # If there are still some predictors
                if size(X, 2) > 0
                    # If some (dummy) predictors were removed
                    if size(X, 2) < origNCol
                        # Calculate the difference
                        diff = origNCol - size(X, 2)

                        # Log an event explaining why the predictors where removed
                        push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                    end

                    # Impute the missing data by the specified method
                    if methods[yVar] == "pmm"
                        workingData[yVar][j][whereY] = pmmImpute!(workingData[yVar][j][.!whereY], X, whereY, whereCount, yVar, iterCounter, j, loggedEvents; kwargs...)
                    elseif methods[yVar] == "norm"
                        workingData[yVar][j][whereY] = normImpute!(workingData[yVar][j][.!whereY], X, whereY, whereCount, yVar, iterCounter, j, loggedEvents; kwargs...)
                    elseif methods[yVar] == "rf"
                        workingData[yVar][j][whereY] = rfImpute!(workingData[yVar][j], X, whereY; kwargs...)
                    end
                else
                    # Log an event explaining why the imputation was skipped
                    push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                end

                # Update mean and variance traces
                updateTraces!(meanTraces, varTraces, workingData[yVar][j][whereY], i, iterCounter, j)

                # Convert to dummy variables if necessary
                if workingData[yVar][j] isa CategoricalArray || nonmissingtype(eltype(workingData[yVar][j])) <: Union{AbstractString, CategoricalValue}
                    workingDataPacified[yVar][j] = pacifyWorkingData(workingData[yVar][j], workingDataLevels[yVar])
                end

                if progressReports
                    progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                    progressRound = floor(Int8, progress / 10)
                    miceEmojis = string(repeat("🐁", progressRound), repeat("🐭", 10 - progressRound))
                    @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress length(loggedEvents)
                end
            end
        else
            # No method specified
            if methods[yVar] == ""
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no method specified.")
            # Invalid method specified
            elseif !(methods[yVar] ∈ ["mean", "norm", "pmm", "rf", "sample"])
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
