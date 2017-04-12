
########################################
# General setup

# Directory where sbatch-r.sh, sbatch-rmd.sh, etc. can be found.
SCRIPT_DIR=scripts

# Directory to store command results.
OUTPUT_DIR=output

# How do we want to run tasks? Can be slurm or bash currently.
# Use SLURM by default, but support running directly in R
# e.g. we run in BASH: "export USE_JOBS=shell" and code will work on bluevelvet.
ifndef USE_JOBS
  # TODO: detect automatically based on sbatch being found in path.
  # Other possible values: shell
	USE_JOBS=slurm
endif

########################################
# Savio configuration.

# This allows us to use environmental variables to override this default.
ifndef ACCOUNT
	ACCOUNT=co_biostat
endif

# This allows us to use environmental variables to override this default.
ifndef PARTITION
	PARTITION=savio2
endif

# This allows us to override the default QOS by setting an environmental variable.
# e.g. we run in BASH: "export QOS=biostat_normal"
ifndef QOS
	# Choose one QOS and comment out the other, or use environmental variables.
	QOS=biostat_savio2_normal
	#QOS=savio_lowprio
endif

########################################
# Execution engines.

# Sbatch runs a SLURM job, e.g. on Savio or XSEDE.
SBATCH=sbatch -A ${ACCOUNT} -p ${PARTITION} --qos ${QOS}

# Setup R to run commands in the background and keep running after logout.
R=nohup nice -n 19 R CMD BATCH --no-restore --no-save

# TODO: support Sun Grid Engine (SGE) for grizzlybear2.
# Or just convert to batchtools?

########################################
# Tasks that can be run.

# Example job:
#data-prep: 1-data-prep.Rmd
#	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-rmd.sh --file=$< --dir=${OUTPUT_DIR}


# Install necessary packages; only needs to be run once per machine.
setup: setup.R
ifeq (${USE_JOBS},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

# Import 2016 data.
import-2016: import-2016.R
ifeq (${USE_JOBS},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

# Analyze 2016 data using targeted_learning.R
# Depends on import-2016.R results.
analyze-2016: analyze-2016.R
ifeq (${USE_JOBS},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

# Test estimate_att() on single 2016 file.
# Depends on import-2016.R results.
test-2016: test-2016.R
ifeq (${USE_JOBS},slurm)
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}
else
	${R} $< ${OUTPUT_DIR}/$<.out &
endif

# Start a bash session with 2 nodes, for up to 12 hours.
bash:
	srun -A ${ACCOUNT} -p ${PARTITION} -N 2 -t 12:00:00 --pty bash

# Next line ensures that this rule works even if there's a file named "clean".
.PHONY : clean
clean:
	rm -f *.Rout
	rm -f slurm*.out
	rm -f install*.out
	rm -f cache/*
