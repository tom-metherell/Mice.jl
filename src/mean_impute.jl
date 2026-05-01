const MEAN_IMPUTER = Imputer((ydata, X, where_y, wherecount, yvar, itercounter, j, loggedevents; kwargs...) -> begin
    repeat([mean(ydata[.!where_y])], wherecount)
end; requirespredictors = false)

registerimputer!("mean", MEAN_IMPUTER)