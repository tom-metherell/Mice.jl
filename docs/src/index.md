# What is `Mice.jl`?

`Mice.jl` is the Julia equivalent of the R package `mice` by Stef van Buuren, Karin Groothuis-Oudshoorn and collaborators [van_buuren_mice_2011](@cite). It allows you to impute missing values in a dataset using multiple imputation by chained equations (MICE).

Currently, only predictive mean matching (PMM) and Bayesian linear regression are supported as methods. `Mice.jl` also currently does not support hybrid imputation models.

If you want to learn more about multiple imputation, this is not the guide for you. Instead, I recommend consulting ["Flexible Imputation of Missing Data"](https://stefvanbuuren.name/fimd/) by Stef van Buuren (ed.) [van_buuren_flexible_2018](@cite).