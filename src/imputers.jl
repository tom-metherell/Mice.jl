"""
    Imputer

Generic wrapper for an imputation function.

The wrapped function should accept:
`(yData, X, whereY, whereCount, yVar, iterCounter, j, loggedEvents; kwargs...)`
and return the imputed values for `yData[whereY]`.
"""
struct Imputer
    f::Function
    requiresPredictors::Bool
end

Imputer(f::Function; requiresPredictors::Bool = true) = Imputer(f, requiresPredictors)

const IMPUTERS = Dict{String, Imputer}()

function registerImputer!(name::String, imputer::Imputer; overwrite::Bool = true)
    if !overwrite && haskey(IMPUTERS, name)
        throw(ArgumentError("Imputer '$name' is already registered."))
    end
    IMPUTERS[name] = imputer
    return IMPUTERS
end

isSupportedMethod(method::String, imputers::AbstractDict{String, <:Imputer}) = haskey(imputers, method)