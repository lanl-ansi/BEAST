# Bilinear Ensemble Actuation Synthesis Toolkit (BEAST)

The Bilinear Ensemble Actuation Synthesis Toolkit (BEAST) is a computational platform that optimizes time-varying control signals to achieve a specified transfer of states governed by bilinear ensemble systems in which model parameters are subject to uncertainty.  The computational method is based on Legendre polynomial approximation of the ensemble state in parameter space and discretization of the evolution equations in the time domain using a product of matrix exponentials corresponding to zero-order hold controls over the discretized time intervals.  The dynamics are sequentially linearized about control and trajectory iterates to formulate a sequence of quadratic programs for computing perturbations to the control that successively minimize the error to the target state and then minimize control energy until the iteration converges.  The input quantities for the code are state transition and control-to-state matrices, initial and target states, control design time, and compact intervals that express the ranges of uncertainties in the model parameters.  The output quantities are time-varying signals for each control variable and the time evolution of the collection of states corresponding to all possible values of parameters in the uncertainty intervals.

# Publications

L. S. Baker, A. L. P. de Lima, A. Zlotnik, and J.-S. Li, “Convergence of Iterative Quadratic Programming for Robust Fixed-Endpoint Transfer of Bilinear Systems,” in 63nd IEEE Conference on Decision and Control (CDC), pp. 8740–8747, IEEE, 2024.

# License

This code is provided under a BSD license for Bilinear Ensemble Actuation Synthesis Toolkit (BEAST), which was developed at the Center for Nonlinear Studies at Los Alamos National Laboratory under U.S. Department of Energy Contract No. 89233218CNA000001.
