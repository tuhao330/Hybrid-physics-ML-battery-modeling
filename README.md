# Hybrid-physics-ML-battery-modeling
Repository for the implementation of the paper:

["Integrating Physics-Based Modeling with Machine Learning for Lithium-Ion Batteries"](https://www.sciencedirect.com/science/article/pii/S030626192201546X)

by Hao Tu, Scott Moura, Yebin Wang, Huazhen Fang

Applied Energy

DOI: [10.1016/j.apenergy.2022.120289](https://doi.org/10.1016/j.apenergy.2022.120289)

## Repository overview

Various battery models are built for battery voltage prediction. They include:

* PureML: Pure feedforward neural network (FNN) models
* DFN: Doyle-Fuller-Neumann model
* SPMT: Single particle model with thermal dynamics (SPMT)
* SPMTNet: SPMT + FNN
* NDC: Nonlinear double capacitor (NDC) model
* NDCNet: NDC + FNN
* Aging NDCNet: Aging-aware NDCNet model
