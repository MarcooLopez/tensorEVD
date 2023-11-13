#!/usr/bin/env Rscript
#!/mnt/research/quantgen/tools/scripts/Rscript_login_shell_wrapper
#SBATCH --job-name=get_var
#SBATCH --output=../log/%x_%A_%a
#SBATCH --time=03:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=68G
#SBATCH --constraint=intel16
#SBATCH --array=1-140

rm(list=ls())

source("https://raw.githubusercontent.com/MarcooLopez/tensorEVD/main/misc/functions.R")

setwd("/mnt/scratch/quantgen/TENSOR_EVD/pipeline")

job <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

# Read the jobs
load("parms/JOBS3.RData")

trait <- as.vector(JOBS[job,"trait"])
method <- as.vector(JOBS[job,"method"])
alpha <- as.vector(JOBS[job,"alpha"])
replicate <- as.vector(JOBS[job,"replicate"])

# Load data
infolder <- "output/genomic_prediction/model_comps"
load(paste0(infolder,"/XG.RData"), verbose=T)
load(paste0(infolder,"/XE.RData"), verbose=T)
load(paste0(infolder,"/XGE_",100*alpha,"_",method,".RData"), verbose=T)

# Read results
outfolder <- paste0("output/genomic_prediction/ANOVA/",trait,"/",
                     method,"/alpha_",100*alpha,"/rep_",replicate)

# Get variance explained
VC <- c()
message("G component variance ...")
B <- as.matrix(readBinMat(paste0(outfolder,"/ETA_G_b.bin")))
VC <- rbind(VC, data.frame(source="G", t(get_variance(XG, B))))

message("E component variance ...")
B <- as.matrix(readBinMat(paste0(outfolder,"/ETA_E_b.bin")))
VC <- rbind(VC, data.frame(source="E", t(get_variance(XE, B))))
rm(B, XG, XE); gc()

message("GE component variance ...")
B <- as.matrix(readBinMat(paste0(outfolder,"/ETA_GE_b.bin")))
VC <- rbind(VC, data.frame(source="GE", t(get_variance(XGE, B))))

# Error
load(paste0(outfolder,"/fm.RData"))
VC <- rbind(VC, data.frame(source="Error",mean=as.vector(fm$varE),SD=as.vector(fm$SD.varE)))

message("Saving results ...")
VC <- data.frame(trait, method, alpha, replicate, VC)
save(VC, file=paste0(outfolder,"/VC.RData"))
message("Done ...")
