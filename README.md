# DPMNeal3.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://igutierrezm.github.io/DPMNeal3.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://igutierrezm.github.io/DPMNeal3.jl/dev)
[![Build Status](https://github.com/igutierrezm/DPMNeal3.jl/workflows/CI/badge.svg)](https://github.com/igutierrezm/DPMNeal3.jl/actions)
[![Coverage](https://codecov.io/gh/igutierrezm/DPMNeal3.jl/branch/master/graph/badge.svg?token=8LyURKsV9M)](https://codecov.io/gh/igutierrezm/DPMNeal3.jl)

DPMNeal3.jl provides an interface for Dirichlet process mixture (DPM) models in Julia using Neal's algorithm 3 [1]. Particularly, this package provides:

- A method for performing one iteration of Neal's algorithm 3.

- Several methods for accessing the current state of the chain (e.g. the vector of cluster labels).

## References

[1] Neal, R. (2000). Markov chain sampling methods for Dirichlet process mixture models. *Journal of Computational and Graphical Statistics*, 9(2). [https://doi.org/10.1080/10618600.2000.10474879](https://doi.org/10.1080/10618600.2000.10474879)
