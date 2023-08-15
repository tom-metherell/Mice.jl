module Mice

    # Dependencies
    using CategoricalArrays, DataFrames, Distributions, LinearAlgebra, NamedArrays, Plots, Printf, Random, StatsBase, Statistics

    # All functions except the main 'mice' function are defined in this file
    include("micehelperfunctions.jl")

"""
    Mids

A multiply imputed dataset object.

The data originally supplied are stored as `data`.

The imputed data are stored as `imputations` (one column per imputation).

The number of imputations is stored as `m`.

The imputation method for each variable is stored as `methods`.

The predictor matrix is stored as `predictorMatrix`.

The order in which the variables are imputed is stored as `visitSequence`.

The number of iterations is stored as `iter`.

The mean of each variable across the imputations is stored as `meanTraces`.

The variance of each variable across the imputations is stored as `varTraces`.
"""
    struct Mids
        data::DataFrame
        imputations::Vector{Matrix}
        m::Int
        methods::NamedArray
        predictorMatrix::NamedArray
        visitSequence::Vector{String}
        iter::Int
        meanTraces::Vector{Matrix}
        varTraces::Vector{Matrix}
    end

"""
    mice(
        data::DataFrame;
        m::Int = 5,
        visitSequence = "monotone",
        methods = nothing,
        predictorMatrix = nothing,
        iter::Int = 10,
        args...
    )

Imputes missing values in a dataset using the MICE algorithm. 
Heavily based on the R package `mice` (van Buuren & Groothuis-Oudshoorn, 2011).

The data containing missing values (`data`) must be supplied as a `DataFrame`.

The number of imputations created is specified by `m`.

The variables will be imputed in the order specified by `visitSequence`. 
The default is sorted by proportion of missing data in descending order ("monotone"); 
the order can be customised using a vector of variable names in the desired order.

The imputation method for each variable is specified by the `NamedArray` `methods`. 
The default is to use predictive mean matching (`pmm`) for all variables. 
Currently only `pmm` is supported. 
Any variable not to be imputed can be marked as such using an empty string ("").

The predictor matrix is specified by the `NamedArray` `predictorMatrix`. 
The default is to use all other variables as predictors for each variable. 
Any variable not predicting another variable can be marked as such in the matrix using a 0.

The number of iterations is specified by `iter`.

If `progressReports` is `true`, a progress indicator will be displayed in the console.
"""
    function mice(
        data::DataFrame;
        m::Int = 5,
        visitSequence = "monotone",
        methods = nothing,
        predictorMatrix = nothing,
        iter::Int = 10,
        progressReports::Bool = true,
        args...
        )

        if visitSequence === "monotone"
            visitSequence = makeMonotoneSequence(data)
        end

        if methods === nothing
            methods = makeMethods(data)
        end

        if predictorMatrix === nothing
            predictorMatrix = makePredictorMatrix(data)
        end

        imputations = initialiseImputations(data, m, visitSequence, methods)

        meanTraces = initialiseTraces(visitSequence, iter, m)
        varTraces = initialiseTraces(visitSequence, iter, m)

        if(progressReports)
            @printf "======= MICE progress =======\n"
        end

        for iterCounter in 1:iter
            for i in eachindex(visitSequence)
                yVar = visitSequence[i]
                y = data[:, yVar]
                X = data[:, predictorMatrix[:, yVar]]

                X = removeLinDeps(X, y)

                if methods[yVar] == "pmm" && any(ismissing.(y)) && !isnothing(X) && size(X, 2) > 0
                    for j in 1:m
                        for k in findall(predictorMatrix[:, yVar])
                            xVar = names(predictorMatrix)[2][k]
                            kVS = findfirst(visitSequence .== xVar)
                            if any(ismissing.(X[:, xVar]))
                                X[ismissing.(X[:, xVar]) .== 1, xVar] = imputations[kVS][:, j]
                            end
                        end
                        imputations[i][:, j] = pmmImpute(y, X, 5, 1e-5)
                    
                        plottingData = deepcopy(data[:, yVar])
                        plottingData[ismissing.(plottingData) .== 1] = imputations[i][:, j]
                    
                        if plottingData isa CategoricalArray || nonmissingtype(eltype(plottingData)) <: AbstractString
                            mapping = Dict(levels(plottingData)[i] => i-1 for i in eachindex(levels(plottingData)))
                            plottingData = [mapping[v] for v in plottingData]
                        end
                    
                        meanTraces[i][iterCounter, j] = mean(plottingData)
                        varTraces[i][iterCounter, j] = var(plottingData)
                        if(progressReports)
                            progress = ((iterCounter - 1)/iter + ((i-1)/length(visitSequence))/iter + (j/m)/length(visitSequence)/iter) * 100
                            miceEmojis = string(repeat("🐁", floor(Int8, progress/10)), repeat("🐭", ceil(Int8, (100 - progress)/10)))
                            @printf "\33[2KIteration:  %u / %u\n\33[2KVariable:   %u / %u (%s)\n\33[2KImputation: %u / %u\n\33[2K%s   %.1f %%\n=============================\u1b[A\u1b[A\u1b[A\u1b[A\r" iterCounter iter i length(visitSequence) yVar j m miceEmojis progress
                        end
                    end
                end
            end
        end

        if(progressReports)
            @printf "\u1b[A\33[2K"
        end

        midsObj = Mids(
            data,
            imputations,
            m,
            methods,
            predictorMatrix,
            visitSequence,
            iter,
            meanTraces,
            varTraces
        )

        return midsObj
    end

    import Plots.plot

    function plot(
        mids::Mids,
        var::String
        )

        var_no = findfirst(mids.visitSequence .== var)

        a = plot(mids.meanTraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(mids.varTraces[var_no], xlabel = "Iteration", ylabel = "Variance")

        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    function plot(
        mids::Mids,
        var_no::Int
        )

        var = mids.visitSequence[var_no]

        a = plot(mids.meanTraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(mids.varTraces[var_no], xlabel = "Iteration", ylabel = "Variance")

        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    export Mids, mice, plot
end