#! /bin/bash
# mothurPCoA.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export DIST=${@:?ERROR: Need to define DIST.}
# export ACCNOS=${2:-} # Can provide a file of sequences to remove but isn't necessary



###################
# PCoA Ordination #
###################

# Calculating PCoA ordination
for FILE in $(echo "${DIST[@]}"); do

	# Pulling distance type from file name
	DISTTYPE=$(echo "${FILE}" | sed 's/.*\.\([a-z]*\)\.0.*/\1/')

	echo PROGRESS: Calculating PCoA ordination and metrics using "${DISTTYPE}" distance.

	# Run diversity analysis on new aligned data set
	mothur "#pcoa(phylip="${FILE}")"

done


# # If $ACCNOS isn't provided
# if [ -z "${ACCNOS}" ]; then

# 	# Calculate axes for the whole distance file
# 	mothur "#pcoa(phylip="${DIST}")"

# else

# 	# Use $ACCNOS to remove sequences from the distance file then calculate axes
# 	mothur "#remove.dists(phylip="${DIST}", accnos="${ACCNOS}");
# 		pcoa(phylip=current)"

# 	FILEPATH=$(echo "${DIST}" | sed 's;\(.*\/\).*;\1;')
# 	REMGROUP=$(echo "${ACCNOS}" | sed 's/.*\.final\.\([a-z1-9]*\).*/\1/')

# 	mv "${FILEPATH}"/*.pick.dist $(echo "${FILEPATH}"/*.pick.dist | sed 's/\.pick\./&'"${REMGROUP}"'\./')

# 	for FILE in $(ls "${FILEPATH}" | grep "\.pick\.pcoa\.[a-z]*$"); do

# 		NEWNAME=$(echo "${FILE}" | sed 's/\.pick\.pcoa\./&'"${REMGROUP}"'\./')
# 		mv "${FILEPATH}"/"${FILE}" "${FILEPATH}"/"${NEWNAME}"

# 	done

# fi


