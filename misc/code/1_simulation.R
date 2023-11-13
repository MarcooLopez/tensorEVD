#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=simulation
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=84G
#SBATCH --constraint=intel18
#SBATCH --array=1-810

rm(list=ls())

library(SFSI)
library(tensorEVD)
library(data.table)

# Set working directory
setwd("/mnt/scratch/quantgen/TENSOR_EVD/pipeline")

# Reading the job id from the array
job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

eps <- .Machine$double.eps  # a numeric zero

# Read the jobs
load("parms/JOBS1.RData")

alpha <- as.vector(JOBS[job,"alpha"])
nG <- as.vector(JOBS[job,"nG"])
nE <- as.vector(JOBS[job,"nE"])
n <- as.vector(JOBS[job,"n"])
replicate <- as.vector(JOBS[job,"replicate"])

# Read data from G2F
load("data/GRM.RData")
load("data/ERM.RData")

# Select subsets of the data at random
seeds <- as.integer(seq(1E3, 1E8, length=2000))
set.seed(seeds[job])

# 1. Select hybrids
index <- sample(1:nrow(KG), size=nG, replace=FALSE)
K1 <- KG[index,index]

# 2. Select environments
index <- sample(1:nrow(KE), size=nE, replace=FALSE)
K2 <- KE[index,index]

# 3. Expand genotypes and environments according to sample size
ID1 <- sample(rownames(K1), size=n, replace=TRUE)
ID2 <- sample(rownames(K2), size=n, replace=TRUE)

nComb <- length(unique(paste0(ID1,":",ID2))) # number of unique combinations

# Hadamard (tensor) matrix (using IDs)
K <- Hadamard(K1, K2, ID1, ID2)

# ============================================
# 1) EVD of full Hadamard matrix
# ============================================
time_eigen <- as.numeric(system.time({
  EVD <- eigen(K, symmetric=TRUE)
  cumvar <- cumsum(EVD$values/sum(EVD$values))
  tmp <- abs(cumvar-alpha)
  index <- 1:min(which(abs(tmp - min(tmp)) <= eps))
  PC <- sweep(EVD$vectors[,index],2,sqrt(EVD$values[index]),FUN="*")
})[3])

# Number of eigenvalues>0
tmp <- abs(cumvar-1)
nPosEigen <- min(which(abs(tmp - min(tmp)) <= eps))

nPC_eigen <- ncol(PC)
pPC_eigen <- nPC_eigen/nPosEigen
K_eigen <- tcrossprod(PC)
rm(PC, EVD); gc()

# ============================================
# 2) EVD tensor
# ============================================
time_tensorEVD <- as.numeric(system.time({
  EVD <- tensorEVD(K1, K2, ID1, ID2, alpha=alpha)
  PC <- sweep(EVD$vectors,2,sqrt(EVD$values),FUN="*")
})[3])

nPC_tensorEVD <- ncol(PC)
pPC_tensorEVD <- nPC_tensorEVD/nPosEigen
K_tensorEVD <- tcrossprod(PC)
rm(PC, EVD); gc()

# Convert to correlation matrix
cov2cor2(K, inplace=TRUE)
cov2cor2(K_eigen, inplace=TRUE)
cov2cor2(K_tensorEVD, inplace=TRUE)

# Frobenius norm
Frobenius_eigen <- sqrt(sum((K-K_eigen)^2))
Frobenius_tensorEVD <- sqrt(sum((K-K_tensorEVD)^2))

# Correlation Matrix Distance
F_K <- sqrt(sum((K^2)))
F_K_eigen <- sqrt(sum((K_eigen^2)))
F_K_tensorEVD <- sqrt(sum((K_tensorEVD^2)))

# a) Full rank vs eigen
trace1 <- 0
for(j in 1:n){
  trace1 <- trace1 + crossprod(K[,j],K_eigen[,j])[1,1]
}
CMD_eigen <- 1 - trace1/(F_K*F_K_eigen)

# b) Full rank vs tensorEVD
trace2 <- 0
for(j in 1:n){
  trace2 <- trace2 +  crossprod(K[,j],K_tensorEVD[,j])[1,1]
}
CMD_tensorEVD <- 1 - trace2/(F_K*F_K_tensorEVD)

outfolder <- "output/simulation"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

tmp <- data.frame(job, alpha, nG, nE, n, replicate, nComb, nPosEigen,
                  time_eigen, Frobenius_eigen, CMD_eigen, nPC_eigen, pPC_eigen,
                  time_tensorEVD, Frobenius_tensorEVD, CMD_tensorEVD, nPC_tensorEVD, pPC_tensorEVD)

outfile <- paste0(outfolder,"/simulation_results.txt")
fwrite(tmp, file=outfile, col.names=!file.exists(outfile), append=TRUE)
