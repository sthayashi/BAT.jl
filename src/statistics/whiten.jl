# This file is a part of BAT.jl, licensed under the MIT License (MIT).


"""
    abstract type WhiteningAlgorithm

Abstract type for sample whitening algorithms.
"""
abstract type WhiteningAlgorithm end
export WhiteningAlgorithm


"""
    struct NoWhitening <: WhiteningAlgorithm

No-op whitening transformation, leaves samples unchanged.

Constructors:

* ```$(FUNCTIONNAME)()```
"""
struct NoWhitening <: WhiteningAlgorithm end
export NoWhitening



"""
    struct CholeskyWhitening <: WhiteningAlgorithm

Whitening transformation based on a Cholesky transformation of the empirical
sample covariance matrix.

Constructors:

* ```$(FUNCTIONNAME)()```
"""
struct CholeskyWhitening <: WhiteningAlgorithm end
export CholeskyWhitening



"""
    struct CholeskyPartialWhitening <: WhiteningAlgorithm

Whitening transformation based on a Cholesky transformation of the empirical
sample covariance matrix.

Only transforms dimensions (degrees of freedom) for which the marginalized
distribution asymptotically approaches zero in the positive and negative
direction.

Constructors:

* ```$(FUNCTIONNAME)()```
"""
struct CholeskyPartialWhitening <: WhiteningAlgorithm end
export CholeskyPartialWhitening


"""
    struct StatisticalWhitening <: WhiteningAlgorithm

Whitening transformation based statistical whitening.
CholeskyPartialWhitening

Constructors:

* ```$(FUNCTIONNAME)()```
"""
struct StatisticalWhitening <: WhiteningAlgorithm end
export StatisticalWhitening



# ToDo: Move whitening code from AHMI here, add new public API
# function `bat_whiten`. Remove `_amhi_whitening_func` and
# change AHMI to use `bat_whiten` directly.

_amhi_whitening_func(algorithm::NoWhitening) = no_whitening
_amhi_whitening_func(algorithm::CholeskyWhitening) = cholesky_whitening!
_amhi_whitening_func(algorithm::CholeskyPartialWhitening) = cholesky_partial_whitening!
_amhi_whitening_func(algorithm::StatisticalWhitening) = statistical_whitening!
