# Bayesian Functional Spatial Partitioning for Single Lesion Detection in 3D

Example is included as `sim_dat.RDS`. The first three elements of the list store spatial coordinates (x,y,z). The fourth element contains the voxel-wise probabilities (mcoord). The last element provides the index for the true region.

The file `BFSP_3D.R` contains the functions and required libraries for implementing the BFSP-3D, as well as functions for computing summary statistics such as lesion volume and dimensions.

The file `run_BFSP_3D.R` generates lesion masks for using the example data.

## bfsp3d.partition()

Estimates region boundaries and Gaussian spatial processes in each region via MCMC.

**Arguments**

*data*: list. The first three elements are vectors of length n giving voxel spatial coordinate x,y,z coordinates. The last element is a vector of voxel-wise values (mcoord).

*iterations*: integer. Total number of MCMC iterations.

*M*: integer. Number of basis functions.

*burn*: integer. Number of MCMC iterations to discard (must be less than iterations).

*vary_centriod*: boolean. If TRUE, the algorithm updates the center at each iteration; if FALSE, assumes the center is fixed.

*center*: vector. Initial center for the algorithm.

*spatial*: boolean. If TRUE, models autocorrelation with exponential gaussian process; if FALSE, assumes data follow independent normal distributions

*initial:* list. Customized initial value. The three elements are vectors of initial values for the boundary parameters (intercept, R, R_minus) and the last element is a vector of initial values for the Gaussian Process parameters.

**Value**

A list containing MCMC samples of estimated parameters

## run_summary()

Computes cluster results and summary statistics (e.g. lesion volume, dimensions) from the model, along with their associated upper and lower quantiles.

**Arguments**

*out*: list. Output from bfsp3d.partition().

*burn*: integer. Number of MCMC iterations to discard (must be less than iterations).

*voxel_size*: vector. Dimensions of the image (x,y,z).

*data*: dataset containing voxel-wise parameter values across the 3D volume.

*name*: vector. Character string specifying the paramter in data from which summary statistics should be computed. 

**Value**

A list with two components:

cluster: clustered results from the model, including the estimated clusters and their associated upper and lower quantiles.

measurement: a list containing elements of summary statistics computed from the clustered results.
