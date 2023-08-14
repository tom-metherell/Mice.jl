module Mice

    # Dependencies
    using CategoricalArrays, DataFrames, Distributions, LinearAlgebra, NamedArrays, Plots, Random, StatsBase, Statistics

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
"""
    function mice(
        data::DataFrame;
        m::Int = 5,
        visitSequence = "monotone",
        methods = nothing,
        predictorMatrix = nothing,
        iter::Int = 10,
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

        select!(data, visitSequence)

        imputations = initialiseImputations(data, m, visitSequence, methods)

        meanTraces = initialiseTraces(visitSequence, iter, m)
        varTraces = initialiseTraces(visitSequence, iter, m)

        for iterCounter in 1:iter
            for i in eachindex(visitSequence)
                yVar = visitSequence[i]
                y = data[:, yVar]
                X = data[:, predictorMatrix[yVar, :]]

                if methods[yVar] == "pmm" && any(ismissing.(y))
                    for j in 1:m
                        for k in findall(predictorMatrix[yVar, :])      
                            xVar = visitSequence[k]
                            replacements = Vector{Union{Missing, nonmissingtype(eltype(X[:, xVar]))}}(missing, size(X, 1))
                            counter = 1
                            for z in axes(X, 1)
                                if ismissing(X[z, xVar])
                                    replacements[z] = imputations[k][counter, j]
                                    counter += 1
                                end
                            end
                            X[:, xVar] = coalesce.(X[:, xVar], replacements)
                        end
                                        
                        imputations[i][:, j] = pmmImpute(y, X, 5)
                    
                        plottingData = deepcopy(data[:, yVar])
                        plottingData[ismissing.(plottingData) .== 1] = imputations[i][:, j]
                    
                        if plottingData isa CategoricalArray
                            mapping = Dict(levels(plottingData)[i] => i for i in eachindex(levels(plottingData)))
                            plottingData = [mapping[v] for v in plottingData]
                        end
                    
                        meanTraces[i][iterCounter, j] = mean(plottingData)
                        varTraces[i][iterCounter, j] = var(plottingData)
                    end
                end
            end
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
        mids::Mids
        )

        plot_grid = plot(layout = (length(mids.meanTraces), 2), legend = false, size = (1200, 400*length(mids.meanTraces)))

        for i in eachindex(mids.meanTraces)
            plot!(plot_grid[i, 1], mids.meanTraces[i], xlabel = "Iteration", ylabel = "Mean")
            plot!(plot_grid[i, 2], mids.varTraces[i], xlabel = "Iteration", ylabel = "Variance")
        end

        return plot_grid
    end

    export Mids, mice, plot
end