"""
    Imputer

Generic wrapper for an imputation function.

The wrapped function should accept:
`(yData, X, whereY, whereCount, (types,) yVar, iterCounter, j, loggedEvents; kwargs...)`
and return the imputed values for `yData[whereY]`. `types` is only passed for two-level imputation methods.
"""
struct Imputer
    f::Function
    requiresPredictors::Bool
    twoLevel::Bool
end

Imputer(f::Function; requiresPredictors::Bool = true, twoLevel::Bool = false) = Imputer(f, requiresPredictors, twoLevel)

const IMPUTERS = Dict{String, Imputer}()

function registerImputer!(name::String, imputer::Imputer; overwrite::Bool = true)
    if !overwrite && haskey(IMPUTERS, name)
        throw(ArgumentError("Imputer '$name' is already registered."))
    end
    IMPUTERS[name] = imputer
    return IMPUTERS
end

isSupportedMethod(method::String, imputers::AbstractDict{String, <:Imputer}) = haskey(imputers, method)