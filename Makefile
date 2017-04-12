#######################
# Setup

######
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

SBATCH=sbatch -A ${ACCOUNT} -p ${PARTITION} --qos ${QOS}

######
# Makefile configuration.
SCRIPT_DIR=scripts
OUTPUT_DIR=output

# Example job:
#data-prep: 1-data-prep.Rmd
#	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-rmd.sh --file=$< --dir=${OUTPUT_DIR}

# Example job:
setup: setup.R
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}

# Import 2016 data.
import-2016: import-2016.R
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}

# Test command-line execution on 2016 data.
# Depends on import-2016.R results.
test-2016: test-2016.R
	${SBATCH} --nodes 1 --job-name=$< ${SCRIPT_DIR}/sbatch-r.sh --file=$< --dir=${OUTPUT_DIR}

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
