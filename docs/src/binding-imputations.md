# Binding imputations together
If you have a number of `Mids` objects that were produced in the same way (e.g. through [multithreading](@ref Multithreading)), you can bind them together into a single `Mids` object using the function `bindImputations`. Note that the log of events might not make sense in the resulting object: it is better to inspect the logs of the individual objects before binding them together.

```@docs
bindImputations
```

```@raw html
<br> <div align="right"> Funded by Wellcome &nbsp;&nbsp;&nbsp; <img src="../wellcome-logo-white.png" style="vertical-align:middle" alt="Wellcome logo" width="50" height="50"> </div>
```