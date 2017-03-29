#!/bin/bash
######### Sbatch configuration.
#
# NOTE: we do not specify account, partition, or QOS in this file,
# in order to allow easier customization. Instead those settings
# should be specified in the command line via the calling file.
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
    file_raw="${i#*=}"
    shopt -s extglob    # Turn on extended pattern support
    # Remove .Rmd if it's included in the filename.
    file=${file_raw%.Rmd}
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

# knitr does not support subdirectories - need to use cd.
cd $dir_output
# This assumes we are in a subdirectory; remove "../" if not.
Rscript -e "knitr::knit('../$file.Rmd', '$file.md')" 2>&1

# Check if the markdown file was generated.
if [ -f "$file.md" ]; then
  # Convert markdown to html
  Rscript -e "markdown::markdownToHTML('$file.md', '$file.html')"
else
  echo "Error: Markdown file $file.md does not exist. Can't create html file."
fi