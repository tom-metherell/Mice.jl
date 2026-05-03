"""
    Imputer

Generic wrapper for an imputation function.

If `ispassive` is `false`, the wrapped function should accept:
`(ydata, X, where_y, wherecount, (types,) yvar, itercounter, j, loggedevents; kwargs...)`
and return only the imputed values for `y`. `types` is only passed for two-level imputation methods.

If `ispassive` is `true`, the wrapped function should accept `(workingdata, where_y, wherecount, yvar, itercounter, j, loggedevents; kwargs...)` and return all values for `y`.
"""
struct Imputer
    f::Function
    requirespredictors::Bool
    twolevel::Bool
    ispassive::Bool
end

Imputer(f::Function; requirespredictors::Bool = true, twolevel::Bool = false, ispassive::Bool = false) = Imputer(f, requirespredictors, twolevel, ispassive)

const IMPUTERS = Dict{String, Imputer}()

"""
    registerimputer!(name::String, imputer::Imputer; overwrite::Bool = true)

Registers a new imputation method.
- `name`: The name of the imputation method (e.g., "norm", "logreg", "2l.norm").
- `imputer`: An instance of `Imputer` wrapping the imputation function.
- `overwrite`: If `false`, will throw an error if an imputer with the same name already exists. If `true`, will overwrite any existing imputer with the same name.
"""
function registerimputer!(name::String, imputer::Imputer; overwrite::Bool = true)
    if !overwrite && haskey(IMPUTERS, name)
        throw(ArgumentError("Imputer '$name' is already registered."))
    end
    IMPUTERS[name] = imputer
    return IMPUTERS
end

isSupportedMethod(method::String, imputers::AbstractDict{String, <:Imputer}) = haskey(imputers, method)