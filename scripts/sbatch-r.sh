#!/bin/bash
######### Sbatch configuration.
#
# NOTE: we do not specify account, partition, or QOS in this file,
# in order to allow easier customization. Instead those settings
# should be specified in the command line via the calling file.
#
###SBATCH --job-name=ck37_vim
#SBATCH --mail-user=ck37@berkeley.edu
###SBATCH --workdir=/global/home/users/alhubbard/Trauma
#
# Job output
#SBATCH --output=slurm.out
#SBATCH --error=slurm.out
#
# Wall clock limit:
#SBATCH --time=48:00:00
#
#### Done configuring sbatch.

# Output to current directory by default. Overriden by --dir option.
dir_output=.

# Extract command line arguments
for i in "$@"
do
case $i in
-f=*|--file=*)
file="${i#*=}"
;;
-d=*|--dir=*)
dir_output="${i#*=}"
;;
esac
done

# Load R if we are using the built-in R module:
# Here we are using a custom compiled version of R, so we don't load the r module.
# module load r

# Load a newer version of gcc than the default.
module load gcc/4.8.5

# Load Java if any R packages need RJava (bartMachine, h2o, etc.)
module load java

# Load a better linear algebra system.
module load lapack

# GPU computation modules. CUDA is 7.5, cudnn is 4.0.
module load cuda cudnn

R CMD BATCH --no-save --no-restore ${file} ${dir_output}/${file}-${SLURM_JOB_ID}.out
