# Hybrid-physics-ML-battery-modeling
Repository for the implementation of the paper:

["Integrating Physics-Based Modeling with Machine Learning for Lithium-Ion Batteries"](https://www.sciencedirect.com/science/article/pii/S030626192201546X)

by Hao Tu, Scott Moura, Yebin Wang, Huazhen Fang

Applied Energy

DOI: [10.1016/j.apenergy.2022.120289](https://doi.org/10.1016/j.apenergy.2022.120289)

## Repository overview

Different models are built for battery voltage prediction. They include:

* PureML: Pure feedforward neural network (FNN) models
* DFN: Doyle-Fuller-Neumann model
* SPMT: Single particle model with thermal dynamics (SPMT)
* SPMTNet: SPMT + FNN
* NDC: Nonlinear double capacitor (NDC) model
* NDCNet: NDC + FNN
* AgingNDCNet: Aging-aware NDCNet model

Data generated for simulation and experimental validation are in the "Dataset" folder.

## Acknowledgement
The authors would like to acknowledge the use of the following repositories for the implementation of the DFN model and SPMT model:

https://github.com/scott-moura/fastDFN

https://github.com/scott-moura/SPMeT

The authors would also like to thank Ning Tian for providing the code of the NDC model.
