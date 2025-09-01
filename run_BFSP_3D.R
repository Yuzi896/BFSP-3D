library(readr)

source("BFSP-3D/bfsp_3D.R")

dat <- readRDS("BFSP-3D/sim_dat.RDS")

center <- c(0.5,0.5,0.5)

result <- bfsp3d.partition(dat,iterations=50000,M=M,burn=20000,vary_centroid=TRUE,center,spatial=FALSE)

run_summary(result, burn = 20000,voxel_size = c(0.067, 0.067, 0.067))