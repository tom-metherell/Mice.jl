function fastMice(
    data::DataFrame;
    m::Int = 5,
    visitSequence::Union{Vector{String}, Nothing} = nothing,
    methods::Union{NamedVector{String}, Nothing} = nothing,
    predictorMatrix::Union{NamedMatrix{Bool}, Nothing} = nothing,
    iter::Int = 10,
    progressReports::Bool = true,
    gcSchedule::Float64 = 1.0,
    kwargs...
    )

    if visitSequence === nothing
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

    loggedEvents = Vector{String}([])

    if(progressReports)
        @printf "======= MICE progress =======\n"
    end

    for iterCounter in 1:iter, i in eachindex(visitSequence)
        fastSampler!(imputations, meanTraces, varTraces, data, m, visitSequence, methods, predictorMatrix, iter, iterCounter, i, progressReports, loggedEvents)
        if(Sys.free_memory()/Sys.total_memory() < gcSchedule)
            GC.gc()
        end
    end

    if(progressReports)
        @printf "\u1b[A\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\n\33[2K\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\u1b[A\r"
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
        varTraces,
        loggedEvents
    )

    return midsObj
end

function fastMice(
    mids::Mids,
    iter::Int = 10,
    progressReports::Bool = true,
    gcschedule::Float64 = 1.0,
    kwargs...
    )

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

    meanTraces = initialiseTraces(visitSequence, iter+prevIter, m)
    for w in eachindex(meanTraces)
        meanTraces[w][1:prevIter, :] = prevMeanTraces[w]
    end

    varTraces = initialiseTraces(visitSequence, iter+prevIter, m)
    for w in eachindex(varTraces)
        varTraces[w][1:prevIter, :] = prevVarTraces[w]
    end
