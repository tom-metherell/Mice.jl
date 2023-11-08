module Mice
    # Dependencies
    using CategoricalArrays: CategoricalArray, levels
    using DataFrames: DataFrame, insertcols!, ncol, Not, select!
    using Distributions: cdf, Chisq, Normal, TDist
    using LinearAlgebra: cholesky, Diagonal, diagm, eigen, inv, qr, rank, svd
    using NamedArrays: NamedArray, NamedMatrix, NamedVector, setnames!
    import Plots: plot
    using Printf: @printf
    using Random: rand, randn, randperm
    using Statistics: cor, mean, quantile, var
    using StatsAPI: coef, coefnames, nobs, stderror
    using StatsBase: CoefTable, PValue, sample, standardize, UnitRangeTransform, zscore
    import StatsModels: contrasts_matrix, termnames
    using StatsModels: AbstractContrasts, ModelFrame, ModelMatrix, setcontrasts!, term

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
        meanTraces::Vector{Matrix{Float64}}
        varTraces::Vector{Matrix{Float64}}
        loggedEvents::Vector{String}
    end

    include("micehelperfunctions.jl")
    include("with.jl")
    include("pool.jl")

    """
        mice(
            data::DataFrame;
            m::Int = 5,
            visitSequence::Union{Vector{String}, Nothing} = nothing,
            methods::Union{NamedVector{String}, Nothing} = nothing,
            predictorMatrix::Union{NamedMatrix{Bool}, Nothing} = nothing,
            iter::Int = 10,
            progressReports::Bool = true,
            gcSchedule::Float64 = 1.0,
            threads::Bool = true,
            kwargs...
            )

    Imputes missing values in a dataset using the MICE algorithm. 
    Heavily based on the R package `mice` (van Buuren & Groothuis-Oudshoorn, 2011).

    The data containing missing values (`data`) must be supplied as a `DataFrame`.

    The number of imputations created is specified by `m`.

    The variables will be imputed in the order specified by `visitSequence`. 
    The default is sorted by proportion of missing data in ascending order; 
    the order can be customised using a vector of variable names in the desired order.

    The imputation method for each variable is specified by the `NamedArray` `methods`. 
    The default is to use predictive mean matching (`pmm`) for all variables. 
    Currently only `pmm` is supported. 
    Any variable not to be imputed can be marked as such using an empty string ("").

    The predictor matrix is specified by the `NamedArray` `predictorMatrix`. 
    The default is to use all other variables as predictors for each variable. 
    Any variable not predicting another variable can be marked as such in the matrix
    using a 0.

    The number of iterations is specified by `iter`.

    If `progressReports` is `true`, a progress indicator will be displayed in the console.

    `gcSchedule` dictates when the garbage collector will be (additionally) invoked. The 
    number provided is the fraction of your RAM remaining at which the GC will be called.
    For small datasets, you may get away with a value of `0.0` (never called), but for larger
    datasets, it may be worthwhile to call it more frequently. The default is to call it 
    after each iteration of each variable (`1.0`), but this may negatively affect
    performance if it is not necessary for your dataset.

    `threads` dictates whether multi-threading will be used. This will improve performance
    for larger jobs if and only if Julia has been launched with multiple threads (which you
    can verify by calling `Threads.nthreads()`). The default is `true`.
    """
    function mice(
        data::DataFrame;
        m::Int = 5,
        visitSequence::Union{Vector{String}, Nothing} = nothing,
        methods::Union{NamedVector{String}, Nothing} = nothing,
        predictorMatrix::Union{NamedMatrix{Bool}, Nothing} = nothing,
        iter::Int = 10,
        progressReports::Bool = true,
        gcSchedule::Float64 = 1.0,
        threads::Bool = true,
        kwargs...
        )

        # If no visit sequence specified: sort by proportion of missing data
        if visitSequence === nothing
            visitSequence = makeMonotoneSequence(data)
        end

        # If no methods specified: use pmm for all variables
        if methods === nothing
            methods = makeMethods(data)
        end

        # If no predictor matrix specified: everything predicts everything
        if predictorMatrix === nothing
            predictorMatrix = makePredictorMatrix(data)
        end

        # Initialise imputations with random draws from the observed data
        imputations = initialiseImputations(data, m, visitSequence, methods)

        # Initialise mean and variance traces (for plotting)
        meanTraces = initialiseTraces(visitSequence, iter, m)
        varTraces = initialiseTraces(visitSequence, iter, m)

        # Initialise log of events
        loggedEvents = Vector{String}([])

        # Print header of progress indicator
        if progressReports
            @printf "======= MICE progress =======\n"
        end

        # For each iteration, for each variable
        for iterCounter in 1:iter, i in eachindex(visitSequence)
            # Run the Gibbs sampler
            sampler!(imputations, meanTraces, varTraces, data, m, visitSequence, methods, predictorMatrix, iter, iterCounter, i, progressReports, loggedEvents, threads)
            
            # If free RAM falls below specified threshold, invoke the garbage collector
            if Sys.free_memory()/Sys.total_memory() < gcSchedule
                GC.gc()
            end
        end

        # Clear the progress indicator
        if progressReports
            if threads
                @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
            else
                @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
            end
        end

        # Define Mids output
        midsObj = Mids(
            data,
            imputations,
            m,
            methods,
            predictorMatrix,
            visitSequence,
            iter,
            meanTraces,
            varTraces,
            loggedEvents
        )

        return midsObj
    end

    """
        mice(
            mids::Mids;
            iter::Int = 10,
            progressReports::Bool = true,
            gcSchedule::Float64 = 1.0,
            threads::Bool = true,
            kwargs...
            )

    Adds additional iterations to an existing `Mids` object.

    The number of *additional* iterations is specified by `iter`.

    If `progressReports` is `true`, a progress indicator will be displayed in the console.

    `gcSchedule` dictates when the garbage collector will be (additionally) invoked. The 
    number provided is the fraction of your RAM remaining at which the GC will be called.
    For small datasets, you may get away with a value of `0.0` (never called), but for larger
    datasets, it may be worthwhile to call it more frequently. The default is to call it 
    after each iteration of each variable (`1.0`), but this may negatively affect
    performance if it is not necessary for your dataset.

    `threads` dictates whether multi-threading will be used. This will improve performance
    for larger jobs if and only if Julia has been launched with multiple threads (which you
    can verify by calling `Threads.nthreads()`). The default is `true`.
    """
    function mice(
        mids::Mids;
        iter::Int = 10,
        progressReports::Bool = true,
        gcSchedule::Float64 = 1.0,
        threads::Bool = true,
        kwargs...
        )

        # Grab existing parameters from the input Mids
        data = mids.data
        imputations = mids.imputations
        m = mids.m
        methods = mids.methods
        predictorMatrix = mids.predictorMatrix
        visitSequence = mids.visitSequence
        prevIter = mids.iter
        prevMeanTraces = mids.meanTraces
        prevVarTraces = mids.varTraces
        loggedEvents = mids.loggedEvents

        # Initialise new mean and variance traces
        meanTraces = initialiseTraces(visitSequence, iter+prevIter, m)
        for w in eachindex(meanTraces)
            meanTraces[w][1:prevIter, :] = prevMeanTraces[w]
        end

        varTraces = initialiseTraces(visitSequence, iter+prevIter, m)
        for w in eachindex(varTraces)
            varTraces[w][1:prevIter, :] = prevVarTraces[w]
        end

        # Print header of progress indicator
        if progressReports
            @printf "======= MICE progress =======\n"
        end
 
        # For each new iteration, for each variable
        for iterCounter in prevIter+1:prevIter+iter, i in eachindex(visitSequence)
            # Run the Gibbs sampler
            sampler!(imputations, meanTraces, varTraces, data, m, visitSequence, methods, predictorMatrix, prevIter+iter, iterCounter, i, progressReports, loggedEvents, threads)
            
            # If free RAM falls below specified threshold, invoke the garbage collector
            if Sys.free_memory()/Sys.total_memory() < gcSchedule
                GC.gc()
            end
        end

        # Clear the progress indicator
        if progressReports
            @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
        end

        # Define the new Mids output
        midsObj = Mids(
            data,
            imputations,
            m,
            methods,
            predictorMatrix,
            visitSequence,
            prevIter+iter,
            meanTraces,
            varTraces,
            loggedEvents
        )
        
        return midsObj
    end

    """
        plot(
            mids::Mids,
            var::String
            )

    Plots the mean and standard deviation of the imputed values for a given variable.
    Here `var` is given as a string (the name of the variable).
    """
    function plot(
        mids::Mids,
        var::String
        )

        # Find index of variable in the visit sequence
        var_no = findfirst(mids.visitSequence .== var)

        # Plot the means and standard deviations across iterations
        a = plot(mids.meanTraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(sqrt.(mids.varTraces[var_no]), xlabel = "Iteration", ylabel = "Standard deviation")

        # Combine plots in a 1x2 grid
        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    """
        plot(
            mids::Mids,
            var_no::Int
            )
        
    Plots the mean and standard deviation of the imputed values for a given variable.
    Here `var_no` is given as an integer (the index of the variable in the `visitSequence`).
    """
    function plot(
        mids::Mids,
        var_no::Int
        )
        
        # Find the variable name in the visit sequence
        var = mids.visitSequence[var_no]

        # Plot the means and standard deviations across iterations
        a = plot(mids.meanTraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(sqrt.(mids.varTraces[var_no]), xlabel = "Iteration", ylabel = "Standard deviation")

        # Combine plots in a 1x2 grid
        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    export makeMethods, makePredictorMatrix, makeVisitSequence, Mids, mice, plot

    include("precompile.jl")
end