Bayesian Functional Spatial Partitioning for Single Lesion Detection in 3D

Example is included as "sim_dat.RDS". The coordinates are saved as the first three element in the list. The fourth element contains the mcoord voxel-wise probabilities. The last element contains the index for the true region.

The file BFSP_3D.R contains the functions and required libraries for implementing the BFSP-3D algorithm, as well as functions for computing summary statistics such as volume and dimensions.

The file "run_BFSP_3D.R" generates lesion masks for example data.

bfsp3d.partition()

data,iterations=10000,M=5,burn=1,vary_centroid=TRUE,center,spatial=FALSE,initial=NULL
Estimates boundaries and gaussian spatial processes in each region via MCMC.

Arguments

data: list, the first three elements are vector of length n giving voxel spatial coordinate x,y,z coordinates and the last element is vector of voxel-wise values (mcoord).

iterations: total number of MCMC iterations

M: number of basis functions 

burn: number of MCMC iterations to discard (must be less than iterations)

vary_centriod: boolean, if TRUE: allow the algorithm to update center within eath iterations, if FALSE: assumes fixed center 

center: the initial center

spatial: boolean, if TRUE: models autocorrelation with exponential gaussian process, if FALSE: assumes data follow independent normal distributions

initial: list, customized initial value where the three elements are vectors of initial values for the boundary parameters (intercept, R, R_minus) and the last element is a vector of initial values for the Gaussian Process parameters    

Value

list containing MCMC samples of estimated parameters

run_summary()

Arguments

out: output from bfsp3d

burn: number of MCMC iterations to discard (must be less than iterations)

voxel_size: dimensions of the image (x,y,z)

data: Data used

Value

list 
