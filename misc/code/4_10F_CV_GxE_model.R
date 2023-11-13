#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=10F_CV
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=20:59:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=96G
#SBATCH --constraint=intel16
#SBATCH --array=1-240

rm(list=ls())

library(BGLR)

setwd("/mnt/scratch/quantgen/TENSOR_EVD/pipeline")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# Read the jobs
load("parms/JOBS4.RData")

trait <- as.vector(JOBS[job,"trait"])
method <- as.vector(JOBS[job,"method"])
fold <- as.vector(JOBS[job,"fold"])
alpha <- as.vector(JOBS[job,"alpha"])

# Load data
load("data/pheno.RData")

infolder <- "output/genomic_prediction/model_comps"
load(paste0(infolder,"/XG.RData"), verbose=T)
load(paste0(infolder,"/XE.RData"), verbose=T)
load(paste0(infolder,"/XGE_",100*alpha,"_",method,".RData"), verbose=T)

tst <- which(PHENO$CV_10fold == fold)
trn <- which(PHENO$CV_10fold != fold)

ETA <- list(
         G=list(X=XG[trn,],model='BRR',saveEffects=FALSE),
         E=list(X=XE[trn,],model='BRR',saveEffects=FALSE),
         GE=list(X=XGE[trn,],model='BRR',saveEffects=FALSE)
       )

y <- PHENO[,trait]

# Output folder
outfolder <- paste0("output/genomic_prediction/10F_CV/",trait,"/",
                    method,"/alpha_",100*alpha)
if(!file.exists(outfolder)) dir.create(outfolder,recursive=TRUE)

# Run the model
nIter <- 18000
burnIn <- 5000

fm <- BLRXy(y=y[trn], ETA=ETA, nIter=nIter, burnIn=burnIn,
            saveAt=paste0(outfolder,"/fold_",fold,"_"))
names(fm$ETA) <- names(ETA)

# Compute Effects in testing
out <- data.frame(PHENO[tst,c("year_loc","genotype","CV_10fold")], y=y[tst])
out$INTERCEPT <- rep(fm$mu,length(tst))
out$G <- as.vector(XG[tst,] %*% fm$ETA$G$b)
out$E <- as.vector(XE[tst,] %*% fm$ETA$E$b)
out$GE <- as.vector(XGE[tst,] %*% fm$ETA$GE$b)
out$yHat <- apply(out[,c("INTERCEPT",names(ETA))],1,sum)

save(out, file=paste0(outfolder,"/results_fold_",fold,".RData"))

unlink(paste0(outfolder,"/fold_",fold,"_*.dat"))
