# The sampler! function includes a ! as it updates workingData, workingDataPacified, meanTraces, varTraces and loggedEvents in place
function sampler!(
    workingData::AxisVector,
    workingDataPacified::AxisVector,
    workingDataLevels::AxisVector,
    meanTraces::Vector{Matrix{Float64}},
    varTraces::Vector{Matrix{Float64}},
    imputeWhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}},
    m::Int,
    visitSequence::Vector{String},
    methods::AxisArray{String, 1, Vector{String}},
    predictorMatrix::AxisArray{Int, 2, Matrix{Int}},
    iter::Int,
    iterCounter::Int,
    i::Int,
    progressReports::Bool,
    loggedEvents::Vector{String};
    imputers::AbstractDict{String, <:Imputer} = IMPUTERS,
    kwargs...
    )
    
    # Grab name of variable to be imputed
    yVar = visitSequence[i]

    # Grab locations of data to be imputed, and set these values to missing (in case of over-imputation)
    whereY = imputeWhere[yVar]
    whereCount = sum(whereY)

    # Grab the names of the predictors
    predictors = axes(predictorMatrix[yVar, :])[1][predictorMatrix[yVar, :] .== 1]

    methodName = methods[yVar]

    if methodName == ""
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no method specified.")
        return
    end

    if !isSupportedMethod(methodName, imputers)
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - method not supported.")
        return
    end

    if !any(whereY)
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no missing data.")
        return
    end

    methodImputer = imputers[methodName]

    if methodImputer.requiresPredictors && isempty(predictors)
        push!(loggedEvents, "Iteration $iterCounter, variable $yVar: imputation skipped - no predictors.")
        return
    end

    for j in 1:m
        X = nothing

        if methodImputer.requiresPredictors
            X = Matrix{Float64}(reduce(hcat, [predictor ∈ axes(workingDataPacified)[1] ? workingDataPacified[predictor][j] : workingData[predictor][j] for predictor in predictors]))
            origNCol = size(X, 2)
            removeLinDeps!(X, workingData[yVar][j], whereY, whereCount)

            if size(X, 2) == 0
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                continue
            end

            if size(X, 2) < origNCol
                diff = origNCol - size(X, 2)
                push!(loggedEvents, "Iteration $iterCounter, variable $yVar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
            end
        end

        workingData[yVar][j][whereY] = methodImputer.f(
            workingData[yVar][j],
            X,
            whereY,
            whereCount,
            yVar,
            iterCounter,
            j,
            loggedEvents;
            kwargs...
        )

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
