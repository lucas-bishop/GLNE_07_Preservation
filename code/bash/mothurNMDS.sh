#! /bin/bash
# mothurNMDS.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export DIST=${@:?ERROR: Need to define DIST.}

# export ACCNOS=$2

# Other variables
export SEED=20170415



###################
# NMDS Ordination #
###################

# Calculating NMDS ordination
for FILE in $(echo "${DIST[@]}"); do

	# Pulling distance type from file name
	DISTTYPE=$(echo "${FILE}" | sed 's/.*\.\([a-z]*\)\.0.*/\1/')

	echo PROGRESS: Calculating NMDS ordination and metrics using "${DISTTYPE}" distance.

	# Calculate axes for the whole distance file
	mothur "#nmds(phylip="${FILE}", seed="${SEED}")"

done


# # If $ACCNOS isn't provided
# if [ -z "${ACCNOS}" ]; then

# 	# Calculate axes for the whole distance file
# 	mothur "#nmds(phylip="${DIST}", seed="${SEED}")"

# else

# 	# Use $ACCNOS to remove sequences from the distance file then calculate axes
# 	mothur "#remove.dists(phylip="${DIST}", accnos="${ACCNOS}");
# 		nmds(phylip=current, seed="${SEED}")"

# 	FILEPATH=$(echo "${DIST}" | sed 's;\(.*\/\).*;\1;')
# 	REMGROUP=$(echo "${ACCNOS}" | sed 's/.*\.final\.\([a-z1-9]*\).*/\1/')

# 	mv "${FILEPATH}"/*.pick.dist $(echo "${FILEPATH}"/*.pick.dist | sed 's/\.pick\./&'"${REMGROUP}"'\./')

# 	for FILE in $(ls "${FILEPATH}" | grep "\.pick\.nmds\.[a-z]*$"); do
# 		NEWNAME=$(echo "${FILE}" | sed 's/\.pick\.nmds\./&'"${REMGROUP}"'\./')
# 		mv "${FILEPATH}"/"${FILE}" "${FILEPATH}"/"${NEWNAME}"
# 	done

# fi


