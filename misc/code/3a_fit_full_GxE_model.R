#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=ANOVA
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=20:59:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=96G
#SBATCH --constraint=intel16
#SBATCH --array=1-140

rm(list=ls())

library(BGLR)

setwd("/mnt/scratch/quantgen/TENSOR_EVD/pipeline")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# Read the jobs
load("parms/JOBS3.RData")

trait <- as.vector(JOBS[job,"trait"])
method <- as.vector(JOBS[job,"method"])
alpha <- as.vector(JOBS[job,"alpha"])
replicate <- as.vector(JOBS[job,"replicate"])

print(JOBS[job,])

# Load data
load("data/pheno.RData")

infolder <- "output/genomic_prediction/model_comps"
load(paste0(infolder,"/XG.RData"), verbose=T)
load(paste0(infolder,"/XE.RData"), verbose=T)
load(paste0(infolder,"/XGE_",100*alpha,"_",method,".RData"), verbose=T)

ETA <- list(
         G=list(X=XG,model='BRR',saveEffects=TRUE),
         E=list(X=XE,model='BRR',saveEffects=TRUE),
         GE=list(X=XGE,model='BRR',saveEffects=TRUE)
       )
rm(XG, XE, XGE); gc()

y <- scale(PHENO[,trait])

# Set seeds for reproducibility
seeds <- as.integer(seq(1E3, 1E8, length=2000))
set.seed(seeds[replicate])

# Output folder
outfolder <- paste0("output/genomic_prediction/ANOVA/",trait,"/",
                    method,"/alpha_",100*alpha,"/rep_",replicate)
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

# Run the model
nIter <- 50000
burnIn <- 5000

fm <- BLRXy(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn, saveAt=paste0(outfolder,"/"))
names(fm$ETA) <- names(ETA)

save(fm, file=paste0(outfolder,"/fm.RData"))
