# TurboLearner
A toolset for the analysis &amp; improvement of CFD simulations

# Description
* Designed to enhance the accuracy of turbulence models
* Samples are generated programatically using different model parameters
* Results of the samples are later analyzed using an algorithm to decide which samples display the best performance against experimental data (or any data deemed correct)
* The algorithm can also be used to optimize simulation times

# Contents
* The `Core` folder contains the code of the OpenFoam solver used in this project, which is a modified third-party implementation of the Myong-Kasagi model, but this can also be easily adapted for other solvers.
* `Simulations/Sampling` contains the sampling code which takes a list of parameters and calls the solver with them, organizing their results when done
* `Simulations/Evaluation` determines which samples show the best matches, for a given set of parameters
* `Simulations/Optimization` applies the same algorithm as the evaluation step, but takes all sample parameters into consideration (useful optimizing simulation times)
