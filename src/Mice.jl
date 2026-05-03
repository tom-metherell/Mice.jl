module Mice
    # Dependencies
    using AxisArrays: axes, AxisArray, AxisMatrix, AxisVector
    using CategoricalArrays: CategoricalArray, CategoricalPool, CategoricalValue, levels
    using Distributions: ccdf, Chisq, FDist, Gamma, Normal, TDist
    using LinearAlgebra: cholesky, diag, Diagonal, diagm, dot, eigen, Hermitian, inv, qr, rank, svd
    using PrecompileTools: @compile_workload
    using Printf: @printf
    using Random: rand, randn, randperm
    import RecipesBase: plot
    using Statistics: cor, mean, quantile, var
    using StatsAPI: coef, coefnames, nobs, stderror
    using StatsBase: CoefTable, PValue, mode, sample, standardize, UnitRangeTransform, zscore
    import StatsModels: contrasts_matrix, termnames
    using StatsModels: AbstractContrasts, ModelFrame, ModelMatrix, setcontrasts!, term
    using Tables: columnnames, columns, columntable, getcolumn, istable

    """
        Mids

    A multiply imputed dataset object.

    The data originally supplied are stored as `data`.

    The imputed data are stored as `imputations` (one column per imputation).

    The locations at which data have been imputed are stored as `imputewhere`.

    The number of imputations is stored as `m`.

    The imputation method for each variable is stored as `methods`.

    The predictor matrix is stored as `predictormatrix`.

    The order in which the variables are imputed is stored as `visitsequence`.

    The number of iterations is stored as `iter`.

    The mean of each variable across the imputations is stored as `meantraces`.

    The variance of each variable across the imputations is stored as `vartraces`.
    """
    struct Mids
        data
        imputations::Vector{Matrix}
        imputewhere::AxisVector{Vector{Bool}}
        m::Int
        methods::AxisVector{String}
        predictormatrix::AxisMatrix{Int}
        visitsequence::Vector{String}
        iter::Int
        meantraces::Vector{Matrix{Union{Missing, Float64}}}
        vartraces::Vector{Matrix{Union{Missing, Float64}}}
        loggedevents::Vector{String}

        function Mids(data, imputations, imputewhere, m, methods, predictormatrix, visitsequence, iter, meantraces, vartraces, loggedevents)
            istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))
            new(data, imputations, imputewhere, m, methods, predictormatrix, visitsequence, iter, meantraces, vartraces, loggedevents)
        end
    end

    include("imputers.jl")
    include("makefunctions.jl")
    include("with.jl")
    include("pool.jl")
    include("sampler.jl")
    include("mean_impute.jl")
    include("norm_impute.jl")
    include("2lnorm_impute.jl")
    include("2lpmm_impute.jl")
    include("2lonlymean_impute.jl")
    include("2lonlynorm_impute.jl")
    include("2lonlypmm_impute.jl")
    include("2lonlymode_impute.jl")
    include("pacify.jl")
    include("pmm_impute.jl")
    include("quantify.jl")
    include("sample_impute.jl")

    """
        mice(
            data;
            m::Int = 5,
            imputewhere::AxisVector{Vector{Bool}} = findmissings(data),
            visitsequence::Vector{String} = makemonotonesequence(imputewhere),
            methods::AxisVector{String} = makemethods(data),
            predictormatrix::AxisMatrix{Int} = makepredictormatrix(data),
            iter::Int = 10,
            progressreports::Bool = true,
            kwargs...
            )

    Imputes missing values in a dataset using the MICE algorithm. 
    The output is a `Mids` object.

    The data containing missing values (`data`) must be supplied as a `Tables.jl` table.

    The number of imputations created is specified by `m`.

    `imputewhere` is an `AxisVector` of boolean vectors specifying where data are to be
    imputed. The default is to impute all missing data.

    The variables will be imputed in the order specified by `visitsequence`. 
    The default is sorted by proportion of missing data in ascending order; 
    the order can be customised using a vector of variable names in the desired order.
    Any column not to be imputed at all can be left out of the visit sequence.

    The imputation method for each variable is specified by the `AxisVector` `methods`. 
    The default is to use predictive mean matching (`pmm`) for all variables.
    Any variable not to be imputed can be marked as such using an empty string ("").

    The predictor matrix is specified by the `AxisMatrix` `predictormatrix`. 
    The default is to use all other variables as predictors for each variable. 
    Any variable not predicting another variable can be marked as such in the matrix
    using a 0.

    The number of iterations is specified by `iter`.

    If `progressreports` is `true`, a progress indicator will be displayed in the console.
    """
    function mice(
        data::T;
        m::Int = 5,
        imputewhere::AxisArray{Vector{Bool}, 1, Vector{Vector{Bool}}} = findmissings(data),
        visitsequence::Vector{String} = makemonotonesequence(imputewhere),
        methods::AxisArray{String, 1, Vector{String}} = makemethods(data),
        predictormatrix::AxisArray{Int, 2, Matrix{Int}} = makepredictormatrix(data),
        iter::Int = 10,
        progressreports::Bool = true,
        imputers::AbstractDict{String, <:Imputer} = IMPUTERS,
        kwargs...
        ) where {T}
        istable(data) || throw(ArgumentError("Data not provided as a Tables.jl table."))

        # If nothing to be imputed: throw error
        if sum(sum.(imputewhere)) == 0
            throw(ArgumentError("Provided dataset contains no missing data to be imputed."))
        end

        # Initialise working data, with imputed locations replaced with random draws from the observed data
        workingdata = initialiseworkingdata(data, imputewhere, m, visitsequence, methods, predictormatrix)

        # Replacing categorical values with dummies where necessary
        workingdatapacified, workingdatalevels = pacifyworkingdata(workingdata)

        # Initialise mean and variance traces (for plotting)
        meantraces = initialisetraces(visitsequence, iter, m)
        vartraces = initialisetraces(visitsequence, iter, m)

        # Initialise log of events
        loggedevents = Vector{String}([])

        # Print header of progress indicator
        if progressreports
            @printf "======= MICE progress =======\n"
        end

        # For each iteration, for each variable
        for itercounter in 1:iter, i in eachindex(visitsequence)
            # Run the Gibbs sampler
            sampler!(workingdata, workingdatapacified, workingdatalevels, meantraces, vartraces, imputewhere, m, visitsequence, methods, predictormatrix, iter, itercounter, i, progressreports, loggedevents; imputers = imputers, kwargs...)
        end

        # Clear the progress indicator
        if progressreports
            @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
        end

        imputations = [reduce(hcat, [workingdata[yvar][j][imputewhere[yvar]] for j in 1:m]) for yvar in visitsequence]

        # Define Mids output
        midsobj = Mids(
            data,
            imputations,
            imputewhere,
            m,
            methods,
            predictormatrix,
            visitsequence,
            iter,
            meantraces,
            vartraces,
            loggedevents
        )

        return midsobj
    end

    """
        mice(
            mids::Mids;
            iter::Int = 10,
            progressreports::Bool = true,
            kwargs...
            )

    Adds additional iterations to an existing `Mids` object.

    The number of *additional* iterations is specified by `iter`.

    `progressreports` can also be specified: all other arguments will be
    ignored or passed to inner functions.
    """
    function mice(
        mids::Mids;
        iter::Int = 10,
        progressreports::Bool = true,
        imputers::AbstractDict{String, <:Imputer} = IMPUTERS,
        kwargs...
        )

        # Grab existing parameters from the input Mids
        data = mids.data
        imputations = mids.imputations
        imputewhere = mids.imputewhere
        m = mids.m
        methods = mids.methods
        predictormatrix = mids.predictormatrix
        visitsequence = mids.visitsequence
        previter = mids.iter
        prevmeantraces = mids.meantraces
        prevvartraces = mids.vartraces
        loggedevents = mids.loggedevents

        # Initialise working data & version with dummy variables
        workingdata = initialiseworkingdata(data, imputations, imputewhere, m, visitsequence, methods, predictormatrix)

        # Replacing categorical values with dummies where necessary
        workingdatapacified, workingdatalevels = pacifyworkingdata(workingdata)

        # Initialise new mean and variance traces
        meantraces = initialisetraces(visitsequence, iter+previter, m)
        for w in eachindex(meantraces)
            meantraces[w][1:previter, :] = prevmeantraces[w]
        end

        vartraces = initialisetraces(visitsequence, iter+previter, m)
        for w in eachindex(vartraces)
            vartraces[w][1:previter, :] = prevvartraces[w]
        end

        # Print header of progress indicator
        if progressreports
            @printf "======= MICE progress =======\n"
        end

        # For each new iteration, for each variable
        for itercounter in previter+1:previter+iter, i in eachindex(visitsequence)
            # Run the Gibbs sampler
            sampler!(workingdata, workingdatapacified, workingdatalevels, meantraces, vartraces, imputewhere, m, visitsequence, methods, predictormatrix, previter+iter, itercounter, i, progressreports, loggedevents; imputers = imputers, kwargs...)
        end

        # Clear the progress indicator
        if progressreports
            @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
        end

        imputations = [reduce(hcat, [workingdata[yvar][j][imputewhere[yvar]] for j in 1:m]) for yvar in visitsequence]

        # Define the new Mids output
        midsobj = Mids(
            data,
            imputations,
            imputewhere,
            m,
            methods,
            predictormatrix,
            visitsequence,
            previter+iter,
            meantraces,
            vartraces,
            loggedevents
        )
        
        return midsobj
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
        var_no = findfirst(mids.visitsequence .== var)

        # Plot the means and standard deviations across iterations
        a = plot(mids.meantraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(sqrt.(mids.vartraces[var_no]), xlabel = "Iteration", ylabel = "Standard deviation")

        # Combine plots in a 1x2 grid
        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    """
        plot(
            mids::Mids,
            var_no::Int
            )
        
    Plots the mean and standard deviation of the imputed values for a given variable.
    Here `var_no` is given as an integer (the index of the variable in the `visitsequence`).
    """
    function plot(
        mids::Mids,
        var_no::Int
        )
        
        # Find the variable name in the visit sequence
        var = mids.visitsequence[var_no]

        # Plot the means and standard deviations across iterations
        a = plot(mids.meantraces[var_no], xlabel = "Iteration", ylabel = "Mean")
        b = plot(sqrt.(mids.vartraces[var_no]), xlabel = "Iteration", ylabel = "Standard deviation")

        # Combine plots in a 1x2 grid
        plot(a, b, layout = (1, 2), legend = false, title = var)
    end

    """
        bindimputations(
            mids1::Mids,
            mids2::Mids
            )

    Combines two `Mids` objects into one. The two objects must have been created from the
    same dataset, with the same imputation methods, predictor matrix, visit sequence and
    number of iterations. The numbers of imputations can be different.
    """
    function bindimputations(
        mids1::Mids,
        mids2::Mids
        )

        data = mids1.data
        if data !== mids2.data
            throw(ArgumentError("Cannot bind these Mids objects: they appear to result from different datasets."))
        end

        imputewhere = mids1.imputewhere
        if imputewhere ≠ mids2.imputewhere
            throw(ArgumentError("Cannot bind these Mids objects: the locations of imputed data are different."))
        end

        methods = mids1.methods
        if methods ≠ mids2.methods
            throw(ArgumentError("Cannot bind these Mids objects: the imputation methods are different."))
        end

        predictormatrix = mids1.predictormatrix
        if predictormatrix ≠ mids2.predictormatrix
            throw(ArgumentError("Cannot bind these Mids objects: the predictor matrices are different."))
        end

        visitsequence = mids1.visitsequence
        if visitsequence ≠ mids2.visitsequence
            throw(ArgumentError("Cannot bind these Mids objects: the visit sequences are different."))
        end

        iter = mids1.iter
        if iter ≠ mids2.iter
            throw(ArgumentError("Cannot bind these Mids objects: the numbers of iterations are different."))
        end

        m = mids1.m + mids2.m
        loggedevents = vcat(mids1.loggedevents, mids2.loggedevents)

        # Initialise new imputations object
        imputations = Vector{Matrix}(undef, length(mids1.imputations))
        # Concatenate imputations
        for i in eachindex(imputations)
            if isassigned(mids1.imputations, i)
                imputations[i] = hcat(mids1.imputations[i], mids2.imputations[i])
            end
        end

        # Initialise new mean and variance traces
        meantraces = Vector{Matrix}(undef, length(mids1.meantraces))
        vartraces = Vector{Matrix}(undef, length(mids1.vartraces))
        # Concatenate traces
        for i in eachindex(meantraces)
            meantraces[i] = hcat(mids1.meantraces[i], mids2.meantraces[i])
        end
        for i in eachindex(vartraces)
            vartraces[i] = hcat(mids1.vartraces[i], mids2.vartraces[i])
        end    
        
        # Define the new Mids output
        midsobj = Mids(
            data,
            imputations,
            imputewhere,
            m,
            methods,
            predictormatrix,
            visitsequence,
            iter,
            meantraces,
            vartraces,
            loggedevents
        )

        return midsobj
    end

    """
        bindimputations(
            midsvector::Vector{Mids}
            )

    Combines a vector of `Mids` objects into one `Mids` object. They must all have been
    created from the same dataset with the same imputation methods, predictor matrix,
    visit sequence and number of iterations. The number of imputations can be different.
    """
    function bindimputations(
        midsvector::Vector{Mids}
        )

        midsobj = midsvector[1]
    
        for i in eachindex(midsvector)[2:end]
            midsobj = bindimputations(midsobj, midsvector[i])
        end

        return midsobj
    end

    """
        bindimputations(
            mids...
            )

    Combines any number of `Mids` objects into one `Mids` object. They must all have been
    created from the same dataset with the same imputation methods, predictor matrix, visit
    sequence and number of iterations. The number of imputations can be different.
    """
    function bindimputations(
        mids...
        )

        midsobj = reduce(bindimputations, mids)

        return(midsobj)
    end

    export bindimputations, complete, findmissings, Imputer, listcomplete, makemethods, makepredictormatrix, mice, Mids, Mipo, Mira, pool, plot, registerimputer!, with

    include("precompile.jl")
end
