#! /bin/bash
# mothurSubsampleShared.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SHAREDLIST=${@} # Takes in all of the shared files supplied via command line

# Setting threshold for minimum number of reads to subsample
export SUBTHRESH=1000



############################
# Subsampling Shared Files #
############################

# # Creating list of all group-specific shared files (ex: sample.final.shared, etc.)
# SHAREDLIST=$(ls "${WORKDIR}"/ | grep "\.final\.shared$")

# Subsampling reads based on count tables
echo PROGRESS: Subsampling shared files.

# Using $SUBTHRESH as threshold for choosing number of reads to subsample
for SHARED in $(echo "${SHAREDLIST[@]}"); do

	# Getting shared file prefix so it can be used to pull the right count file for use in setting subsample size
	PREFIX=$(echo "${SHARED}" | rev | cut -d "." -f 2- | rev) # Removes .shared suffix and stores remaining string

	# Pulling smallest number of reads greater than or equal to $SUBTHRESH for use in subsampling
	READCOUNT=$(awk -v SUBTHRESH="${SUBTHRESH}" '$2 >= SUBTHRESH { print $2}' "${PREFIX}".count.summary | sort -n | head -n 1)

	# If $READCOUNT is empty (no samples with reads above $SUBTHRESH), use the lowest number of reads from count file
	# If set as $SUBTHRESH instead (comment/uncomment as needed), will not generate subsampled file for groups with less than $SUBTHRESH reads
	if [ -z "${READCOUNT}" ]; then
		# READCOUNT="${SUBTHRESH}"
		READCOUNT=$(awk '{ print $2}' "${PREFIX}".count.summary | sort -n | head -n 1)
	fi

	# Debugging message
	echo PROGRESS: Subsampling "${SHARED}" to "${READCOUNT}" reads.

	# Subsampling reads based on $READCOUNT
	mothur "#sub.sample(shared="${SHARED}", size="${READCOUNT}")"

done
