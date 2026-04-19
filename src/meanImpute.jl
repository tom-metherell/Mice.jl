const MEAN_IMPUTER = Imputer((yData, X, whereY, whereCount, yVar, iterCounter, j, loggedEvents; kwargs...) -> begin
    repeat([mean(yData[.!whereY])], whereCount)
end; requiresPredictors = false)

registerImputer!("mean", MEAN_IMPUTER)