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

"""
    registerImputer!(name::String, imputer::Imputer; overwrite::Bool = true)

Registers a new imputation method.
- `name`: The name of the imputation method (e.g., "norm", "logreg", "2l.norm").
- `imputer`: An instance of `Imputer` wrapping the imputation function.
- `overwrite`: If `false`, will throw an error if an imputer with the same name already exists. If `true`, will overwrite any existing imputer with the same name.
"""
function registerImputer!(name::String, imputer::Imputer; overwrite::Bool = true)
    if !overwrite && haskey(IMPUTERS, name)
        throw(ArgumentError("Imputer '$name' is already registered."))
    end
    IMPUTERS[name] = imputer
    return IMPUTERS
end

isSupportedMethod(method::String, imputers::AbstractDict{String, <:Imputer}) = haskey(imputers, method)