# The sampler! function includes a ! as it updates workingdata, workingdatapacified, meantraces, vartraces and loggedevents in place
function sampler!(
    workingdata::AxisVector,
    workingdatapacified::AxisVector,
    workingdatalevels::AxisVector,
    meantraces::Vector{Matrix{Float64}},
    vartraces::Vector{Matrix{Float64}},
    imputewhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}},
    m::Int,
    visitsequence::Vector{String},
    methods::AxisArray{String, 1, Vector{String}},
    predictormatrix::AxisArray{Int, 2, Matrix{Int}},
    iter::Int,
    itercounter::Int,
    i::Int,
    progressreports::Bool,
    loggedevents::Vector{String};
    imputers::AbstractDict{String, <:Imputer} = IMPUTERS,
    kwargs...
    )
    
    # Grab name of variable to be imputed
    yvar = visitsequence[i]

    # Grab locations of data to be imputed, and set these values to missing (in case of over-imputation)
    where_y = imputewhere[yvar]
    wherecount = sum(where_y)

    # Grab the names of the predictors
    predictors = axes(predictormatrix[yvar, :])[1][predictormatrix[yvar, :] .≠ 0]

    methodname = methods[yvar]

    if methodname == ""
        push!(loggedevents, "Iteration $itercounter, variable $yvar: imputation skipped - no method specified.")
        return
    end

    if !isSupportedMethod(methodname, imputers)
        push!(loggedevents, "Iteration $itercounter, variable $yvar: imputation skipped - method not supported.")
        return
    end

    methodimputer = imputers[methodname]
    twolevel = methodimputer.twolevel

    if !any(where_y) && !methodimputer.ispassive
        push!(loggedevents, "Iteration $itercounter, variable $yvar: imputation skipped - no missing data.")
        return
    end

    if methodimputer.requirespredictors && isempty(predictors)
        push!(loggedevents, "Iteration $itercounter, variable $yvar: imputation skipped - no predictors.")
        return
    end

    for j in 1:m
        X = nothing
        types = Int[]

        if methodimputer.ispassive
            workingdata[yvar][j] = methodimputer.f(
                workingdata,
                where_y,
                wherecount,
                yvar,
                itercounter,
                j,
                loggedevents;
                kwargs...
            )
        else
            if methodimputer.requirespredictors
                predictordata = Vector{Any}(undef, length(predictors))

                for p in eachindex(predictors)
                    predictor = predictors[p]
                    predictortype = predictormatrix[yvar, predictor]

                    # For two-level methods, classing variables (coded -2) should remain in their original form.
                    if twolevel && predictortype == -2
                        predictordata[p] = workingdata[predictor][j]
                    elseif predictor ∈ axes(workingdatapacified)[1]
                        predictordata[p] = workingdatapacified[predictor][j]
                    else
                        predictordata[p] = workingdata[predictor][j]
                    end
                end

                X = Matrix{Float64}(reduce(hcat, predictordata))
                orig_ncol = size(X, 2)
                removelindeps!(X, workingdata[yvar][j], where_y, wherecount)

                types = vcat([
                    repeat([predictormatrix[yvar, predictor]], size(predictordata[p], 2))
                    for (p, predictor) in enumerate(predictors)
                ]...)

                if size(X, 2) == 0
                    push!(loggedevents, "Iteration $itercounter, variable $yvar, imputation $j: imputation skipped - all predictors dropped because of high multicollinearity.")
                    continue
                end

                if size(X, 2) < orig_ncol
                    diff = orig_ncol - size(X, 2)
                    push!(loggedevents, "Iteration $itercounter, variable $yvar, imputation $j: $diff (dummy) predictors were dropped because of high multicollinearity.")
                end
            end

            if twolevel
                workingdata[yvar][j][where_y] = methodimputer.f(
                    workingdata[yvar][j],
                    X,
                    where_y,
                    wherecount,
                    types,
                    yvar,
                    itercounter,
                    j,
                    loggedevents;
                    kwargs...
                )
            elseif !methodimputer.ispassive
                if methodimputer.requirespredictors && any(types .≠ 1)
                    push!(loggedevents, "Iteration $itercounter, variable $yvar, imputation $j: imputation skipped - predictor matrix contains unsupported values.")
                    continue
                end

                workingdata[yvar][j][where_y] = methodimputer.f(
                    workingdata[yvar][j],
                    X,
                    where_y,
                    wherecount,
                    yvar,
                    itercounter,
                    j,
                    loggedevents;
                    kwargs...
                )
            end
        end

        updatetraces!(meantraces, vartraces, workingdata[yvar][j][where_y], i, itercounter, j)

        if workingdata[yvar][j] isa CategoricalArray || nonmissingtype(eltype(workingdata[yvar][j])) <: Union{AbstractString, CategoricalValue}
            workingdatapacified[yvar][j] = pacifyworkingdata(workingdata[yvar][j], workingdatalevels[yvar])
        end

        if progressreports
            progress = ((itercounter - 1)/iter + ((i-1)/length(visitsequence))/iter + (j/m)/length(visitsequence)/iter) * 100
            progressround = floor(Int8, progress / 10)
            miceemojis = string(repeat("🐁", progressround), repeat("🐭", 10 - progressround))
            @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n\33[2KLogged events: %u\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r" itercounter iter i length(visitsequence) yvar j m miceemojis progress length(loggedevents)
        end
    end
end

# The updatetraces! function includes a ! as it updates meantraces and vartraces in place
function updatetraces!(
    meantraces::Vector{Matrix{Float64}},
    vartraces::Vector{Matrix{Float64}},
    imputeddata::AbstractArray,
    i::Int,
    itercounter::Int,
    j::Int
    )

    # If the imputed data are categorical
    if imputeddata isa CategoricalArray || nonmissingtype(eltype(imputeddata)) <: Union{AbstractString, CategoricalValue}
        # Convert the imputed data to integers
        mapping = Dict(levels(imputeddata)[i] => i-1 for i in eachindex(levels(imputeddata)))
        imputeddata = [mapping[v] for v in imputeddata]
    end

    # Find the mean and variance and append these to the traces
    meantraces[i][itercounter, j] = mean(imputeddata)
    vartraces[i][itercounter, j] = var(imputeddata)
end
