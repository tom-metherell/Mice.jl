module Mice
    # Dependencies
    using AxisArrays: axes, AxisArray, AxisMatrix, AxisVector
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Distributions: ccdf, Chisq, FDist, Normal, TDist
    using LinearAlgebra: cholesky, Diagonal, diagm, eigen, Hermitian, inv, qr, rank, svd
    using PrecompileTools: @compile_workload
    using Printf: @printf
    using Random: rand, randn, randperm
    import RecipesBase: plot
    using Statistics: cor, mean, quantile, var
    using StatsAPI: coef, coefnames, nobs, stderror
    using StatsBase: CoefTable, PValue, sample, standardize, UnitRangeTransform, zscore
    import StatsModels: contrasts_matrix, termnames
    using StatsModels: AbstractContrasts, ModelFrame, ModelMatrix, setcontrasts!, term
    using Tables: columnnames, columns, columntable, getcolumn, istable

    """
        Mids

    A multiply imputed dataset object.

    The data originally supplied are stored as `data`.

    The imputed data are stored as `imputations` (one column per imputation).

    The locations at which data have been imputed are stored as `imputeWhere`.

    The number of imputations is stored as `m`.

    The imputation method for each variable is stored as `methods`.

    The predictor matrix is stored as `predictorMatrix`.

    The order in which the variables are imputed is stored as `visitSequence`.

    The number of iterations is stored as `iter`.

    The mean of each variable across the imputations is stored as `meanTraces`.

    The variance of each variable across the imputations is stored as `varTraces`.
    """
    struct Mids
        data
        imputations::Vector{Matrix}
        imputeWhere::AxisVector{Vector{Bool}}
        m::Int
        methods::AxisVector{String}
        predictorMatrix::AxisMatrix{Int}
        visitSequence::Vector{String}
        iter::Int
        meanTraces::Vector{Matrix{Union{Missing, Float64}}}
        varTraces::Vector{Matrix{Union{Missing, Float64}}}
        loggedEvents::Vector{String}

        function Mids(data, imputations, imputeWhere, m, methods, predictorMatrix, visitSequence, iter, meanTraces, varTraces, loggedEvents)
            istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))
            new(data, imputations, imputeWhere, m, methods, predictorMatrix, visitSequence, iter, meanTraces, varTraces, loggedEvents)
        end
    end

    include("makeFunctions.jl")
    include("with.jl")
    include("pool.jl")
    include("sampler.jl")
    include("normImpute.jl")
    include("pacify.jl")
    include("pmmImpute.jl")
    include("quantify.jl")
    include("rfImpute.jl")
    include("sampleImpute.jl")

    """
        mice(
            data;
            m::Int = 5,
            imputeWhere::AxisVector{Vector{Bool}} = findMissings(data),
            visitSequence::Vector{String} = makeMonotoneSequence(imputeWhere),
            methods::AxisVector{String} = makeMethods(data),
            predictorMatrix::AxisMatrix{Int} = makePredictorMatrix(data),
            iter::Int = 10,
            progressReports::Bool = true,
            gcSchedule::Float64 = 0.3,
            kwargs...
            )

    Imputes missing values in a dataset using the MICE algorithm. 
    The output is a `Mids` object.

    The data containing missing values (`data`) must be supplied as a `Tables.jl` table.

    The number of imputations created is specified by `m`.

    `imputeWhere` is an `AxisVector` of boolean vectors specifying where data are to be
    imputed. The default is to impute all missing data.

    The variables will be imputed in the order specified by `visitSequence`. 
    The default is sorted by proportion of missing data in ascending order; 
    the order can be customised using a vector of variable names in the desired order.
    Any column not to be imputed at all can be left out of the visit sequence.

    The imputation method for each variable is specified by the `AxisVector` `methods`. 
    The default is to use predictive mean matching (`pmm`) for all variables.
    Any variable not to be imputed can be marked as such using an empty string ("").

    The predictor matrix is specified by the `AxisMatrix` `predictorMatrix`. 
    The default is to use all other variables as predictors for each variable. 
    Any variable not predicting another variable can be marked as such in the matrix
    using a 0.

    The number of iterations is specified by `iter`.

    If `progressReports` is `true`, a progress indicator will be displayed in the console.

    `gcSchedule` dictates when the garbage collector will be (additionally) invoked. The 
    number provided is the fraction of your RAM remaining at which the GC will be called.
    For small datasets, you may get away with a value of `0.0` (never called), but for larger
    datasets, it may be worthwhile to call it more frequently. The default is `0.3`, but for
    really large jobs you may want to increase this value.
    """
    function mice(
        data::T;
        m::Int = 5,
        imputeWhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}} = findMissings(data),
        visitSequence::Vector{String} = makeMonotoneSequence(imputeWhere),
        methods::AxisArray{String, 1, Vector{String}} = makeMethods(data),
        predictorMatrix::AxisArray{Int, 2, Matrix{Int}} = makePredictorMatrix(data),
        iter::Int = 10,
        progressReports::Bool = true,
        gcSchedule::Float64 = 0.0,
        kwargs...
        ) where {T}
        istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

        # If nothing to be imputed: throw error
        if sum(sum.(imputeWhere)) == 0
            throw(ArgumentError("Provided dataset contains no missing data to be imputed."))
        end

        # Initialise working data, with imputed locations replaced with random draws from the observed data
        workingData = initialiseWorkingData(data, imputeWhere, m, visitSequence, methods, predictorMatrix)

        # Replacing categorical values with dummies where necessary
        workingDataPacified, workingDataLevels = pacifyWorkingData(workingData)

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
            sampler!(workingData, workingDataPacified, workingDataLevels, meanTraces, varTraces, imputeWhere, m, visitSequence, methods, predictorMatrix, iter, iterCounter, i, progressReports, 
loggedEvents; kwargs...)
            
            # If free RAM falls below specified threshold, invoke the garbage collector
            if Sys.free_memory()/Sys.total_memory() < gcSchedule
                GC.gc()
            end
        end

        # Clear the progress indicator
        if progressReports
            @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
        end

        imputations = [reduce(hcat, [workingData[yVar][j][imputeWhere[yVar]] for j in 1:m]) for yVar in visitSequence]

        # Define Mids output
        midsObj = Mids(
            data,
            imputations,
            imputeWhere,
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
            gcSchedule::Float64 = 0.3;
            kwargs...
            )

    Adds additional iterations to an existing `Mids` object.

    The number of *additional* iterations is specified by `iter`.

    `progressReports` and `gcSchedule` can also be specified: all other arguments will be
    ignored.
    """
    function mice(
        mids::Mids;
        iter::Int = 10,
        progressReports::Bool = true,
        gcSchedule::Float64 = 0.3,
        kwargs...
        )

        # Grab existing parameters from the input Mids
        data = mids.data
        imputations = mids.imputations
        imputeWhere = mids.imputeWhere
        m = mids.m
        methods = mids.methods
        predictorMatrix = mids.predictorMatrix
        visitSequence = mids.visitSequence
        prevIter = mids.iter
        prevMeanTraces = mids.meanTraces
        prevVarTraces = mids.varTraces
        loggedEvents = mids.loggedEvents

        # Initialise working data & version with dummy variables
        workingData = initialiseWorkingData(data, imputations, imputeWhere, m, visitSequence, methods, predictorMatrix)

        # Replacing categorical values with dummies where necessary
        workingDataPacified, workingDataLevels = pacifyWorkingData(workingData)

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
            sampler!(workingData, workingDataPacified, workingDataLevels, meanTraces, varTraces, imputeWhere, m, visitSequence, methods, predictorMatrix, prevIter+iter, iterCounter, i, progressReports, 
loggedEvents; kwargs...)
            
            # If free RAM falls below specified threshold, invoke the garbage collector
            if Sys.free_memory()/Sys.total_memory() < gcSchedule
                GC.gc()
            end
        end

        # Clear the progress indicator
        if progressReports
            @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
        end

        imputations = [reduce(hcat, [workingData[yVar][j][imputeWhere[yVar]] for j in 1:m]) for yVar in visitSequence]

        # Define the new Mids output
        midsObj = Mids(
            data,
            imputations,
            imputeWhere,
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

    """
        bindImputations(
            mids1::Mids,
            mids2::Mids
            )

    Combines two `Mids` objects into one. The two objects must have been created from the
    same dataset, with the same imputation methods, predictor matrix, visit sequence and
    number of iterations. The numbers of imputations can be different.
    """
    function bindImputations(
        mids1::Mids,
        mids2::Mids
        )

        data = mids1.data
        if data !== mids2.data
            throw(ArgumentError("Cannot bind these Mids objects: they appear to result from different datasets."))
        end

        imputeWhere = mids1.imputeWhere
        if imputeWhere ≠ mids2.imputeWhere
            throw(ArgumentError("Cannot bind these Mids objects: the locations of imputed data are different."))
        end

        methods = mids1.methods
        if methods ≠ mids2.methods
            throw(ArgumentError("Cannot bind these Mids objects: the imputation methods are different."))
        end

        predictorMatrix = mids1.predictorMatrix
        if predictorMatrix ≠ mids2.predictorMatrix
            throw(ArgumentError("Cannot bind these Mids objects: the predictor matrices are different."))
        end

        visitSequence = mids1.visitSequence
        if visitSequence ≠ mids2.visitSequence
            throw(ArgumentError("Cannot bind these Mids objects: the visit sequences are different."))
        end

        iter = mids1.iter
        if iter ≠ mids2.iter
            throw(ArgumentError("Cannot bind these Mids objects: the numbers of iterations are different."))
        end

        m = mids1.m + mids2.m
        loggedEvents = vcat(mids1.loggedEvents, mids2.loggedEvents)

        # Initialise new imputations object
        imputations = Vector{Matrix}(undef, length(mids1.imputations))
        # Concatenate imputations
        for i in eachindex(imputations)
            if isassigned(mids1.imputations, i)
                imputations[i] = hcat(mids1.imputations[i], mids2.imputations[i])
            end
        end

        # Initialise new mean and variance traces
        meanTraces = Vector{Matrix}(undef, length(mids1.meanTraces))
        varTraces = Vector{Matrix}(undef, length(mids1.varTraces))
        # Concatenate traces
        for i in eachindex(meanTraces)
            meanTraces[i] = hcat(mids1.meanTraces[i], mids2.meanTraces[i])
        end
        for i in eachindex(varTraces)
            varTraces[i] = hcat(mids1.varTraces[i], mids2.varTraces[i])
        end    
        
        # Define the new Mids output
        midsObj = Mids(
            data,
            imputations,
            imputeWhere,
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
        bindImputations(
            midsVector::Vector{Mids}
            )

    Combines a vector of `Mids` objects into one `Mids` object. They must all have been
    created from the same dataset with the same imputation methods, predictor matrix,
    visit sequence and number of iterations. The number of imputations can be different.
    """
    function bindImputations(
        midsVector::Vector{Mids}
        )

        midsObj = midsVector[1]
    
        for i in eachindex(midsVector)[2:end]
            midsObj = bindImputations(midsObj, midsVector[i])
        end

        return midsObj
    end

    """
        bindImputations(
            mids...
            )

    Combines any number of `Mids` objects into one `Mids` object. They must all have been
    created from the same dataset with the same imputation methods, predictor matrix, visit
    sequence and number of iterations. The number of imputations can be different.
    """
    function bindImputations(
        mids...
        )

        midsObj = reduce(bindImputations, mids)

        return(midsObj)
    end

    export bindImputations, complete, findMissings, listComplete, makeMethods, makePredictorMatrix, mice, Mids, Mipo, Mira, pool, plot, with

    include("precompile.jl")
end
