#! /bin/bash
# mothurAlpha.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SHARED=${1:?ERROR: Need to define SHARED.}

# Setting threshold for minimum number of reads to subsample
export SUBTHRESH=1000



###############################
# Calculate Diversity Metrics #
###############################

# Setting variables to determine number of reads for subsampling
echo PROGRESS: Setting subsampling parameters.

# Getting shared file prefix so it can be used to pull the right count file for use in setting subsample size
PREFIX=$(echo "${SHARED}" | rev | cut -d "." -f 2- | rev) # Removes .shared suffix and stores remaining string



# Setting number of reads to subsample to
if [ -e "${PREFIX}".count.summary ]; then

	# Pulling smallest number of reads greater than or equal to $SUBTHRESH for use in subsampling 
	READCOUNT=$(awk -v SUBTHRESH="${SUBTHRESH}" '$2 >= SUBTHRESH { print $2}' "${PREFIX}".count.summary | sort -n | head -n 1)

else

	# If a count file doesn't already exist, subsample to $SUBTHRESH reads
	READCOUNT="${SUBTHRESH}"

fi



# Run diversity analysis on new aligned data set
echo PROGRESS: Calculating "${SHARED}" alpha diversity and subsampling to "${READCOUNT}" reads.

# Calculating alpha and beta diversity
# If a sample doesn't have enough reads, it'll be eliminated from the analysis
mothur "#summary.single(shared="${SHARED}", calc=nseqs-coverage-invsimpson-shannon-sobs, subsample="${READCOUNT}")"


