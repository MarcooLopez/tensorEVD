#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=comps
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=06:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=160G
#SBATCH --constraint=intel18
#SBATCH --array=1-20

rm(list=ls())

library(tensorEVD)
library(data.table)

setwd("/mnt/scratch/quantgen/TENSOR_EVD/pipeline")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# Read the jobs
load("parms/JOBS2.RData")

alpha <- as.vector(JOBS[job,"alpha"])
replicate <- as.vector(JOBS[job,"replicate"])

# Load data
load("data/pheno.RData")
load("data/GRM.RData")
load("data/ERM.RData")

outfolder <- "output/genomic_prediction/model_comps"
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

IDG <- as.character(PHENO$genotype)
IDE <- as.character(PHENO$year_loc)

eps <- .Machine$double.eps  # a numeric zero

if(job == 1){ # This code is performed only once for main effects G and E
  message("Genotype data ...")
  EVD <- eigen(KG, symmetric=TRUE)
  rownames(EVD$vectors) <- rownames(KG)
  index <- which(EVD$values>1E-8)
  XG <- sweep(EVD$vectors[,index], 2, sqrt(EVD$values[index]), FUN='*')
  XG <- XG[IDG,]

  save(XG, file=paste0(outfolder,"/XG.RData"))

  message("Environmental data ...")
  EVD <- eigen(KE, symmetric=TRUE)
  rownames(EVD$vectors) <- rownames(KE)
  index <- which(EVD$values>1E-8)
  XE <- sweep(EVD$vectors[,index], 2, sqrt(EVD$values[index]), FUN='*')
  XE <- XE[IDE,]

  save(XE, file=paste0(outfolder,"/XE.RData"))
  rm(XG, XE, EVD); gc()
}

message("Standard EVD for alpha ",alpha,"...")
KGE <- Hadamard(KG, KE, IDG, IDE)
time_eigen <- as.numeric(system.time({
  EVD <- eigen(KGE, symmetric=TRUE)
  cumvar <- cumsum(EVD$values/sum(EVD$values))
  tmp <- abs(cumvar-alpha)
  index <- 1:min(which(abs(tmp - min(tmp)) <= eps))
  XGE <- sweep(EVD$vectors[,index], 2, sqrt(EVD$values[index]), FUN='*')
})[3])

if(replicate == 1){  # Save it just once
  save(XGE, file=paste0(outfolder,"/XGE_",100*alpha,"_eigen.RData"))
}
rm(XGE, EVD); gc()

time_tensorEVD <- NA
if(alpha < 1.0){  # do not run for alpha=1
  message("Tensor EVD for alpha ",alpha,"...")
  time_tensorEVD <- as.numeric(system.time({
    EVD <- tensorEVD(KG, KE, IDG, IDE, alpha=alpha, verbose=TRUE)
    XGE <- sweep(EVD$vectors, 2, sqrt(EVD$values), FUN='*')
  })[3])

  if(replicate == 1){  # Save it just once
    save(XGE, file=paste0(outfolder,"/XGE_",100*alpha,"_tensorEVD.RData"))
  }
}

# Save timing
tmp <- data.frame(alpha, replicate, time_eigen, time_tensorEVD)
outfile <- paste0(outfolder,"/timing_EVD.txt")
fwrite(tmp, file=outfile, col.names=!file.exists(outfile), append=TRUE)
