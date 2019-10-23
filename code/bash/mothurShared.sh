#! /bin/bash
# mothurShared.sh
# William L. Close
# Schloss Lab
# University of Michigan

##################
# Set Script Env #
##################

# Set the variables to be used in this script
export SAMPLEDIR=${1:?ERROR: Need to define SAMPLEDIR.}
export MOCKGROUPS=${2:?ERROR: Need to define MOCKGROUPS.} # List of mock groups in raw data dir separated by '-'
export CONTROLGROUPS=${3:?ERROR: Need to define CONTROLGROUPS.} # List of control groups in raw data dir separated by '-'
export SILVAV4=${4:?ERROR: Need to define SILVAV4.}
export RDPFASTA=${5:?ERROR: Need to define RDPFASTA.}
export RDPTAX=${6:?ERROR: Need to define RDPTAX.}


# Other variables
export OUTDIR=data/mothur/process
export COMBINEDGROUPS=$(echo "${MOCKGROUPS}"-"${CONTROLGROUPS}") # Combines the list of mock and control groups into a single string separated by '-'



###################
# Run QC Analysis #
###################

echo PROGRESS: Assembling, quality controlling, clustering, and classifying sequences.

# Making output dir
mkdir -p "${OUTDIR}"

# Making contigs from fastq.gz files, aligning reads to references, removing any non-bacterial sequences, calculating distance matrix, making shared file, and classifying OTUs
mothur "#make.file(type=gz, inputdir="${SAMPLEDIR}", outputdir="${OUTDIR}");
	make.contigs(file=current);
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275, maxhomop=8);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference="${SILVAV4}");
	screen.seqs(fasta=current, count=current, start=1968, end=11550);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=T);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference="${RDPFASTA}", taxonomy="${RDPTAX}", cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	dist.seqs(fasta=current, cutoff=0.03);
	cluster(column=current, count=current);
	make.shared(list=current, count=current, label=0.03);
	classify.otu(list=current, count=current, taxonomy=current, label=0.03);
	get.oturep(fasta=current, count=current, list=current, label=0.03, method=abundance)"



# Renaming output files for use later
mv "${OUTDIR}"/*.precluster.pick.pick.fasta "${OUTDIR}"/errorinput.fasta
mv "${OUTDIR}"/*.vsearch.pick.pick.count_table "${OUTDIR}"/errorinput.count_table
mv "${OUTDIR}"/*.opti_mcc.shared "${OUTDIR}"/final.shared
mv "${OUTDIR}"/*.cons.taxonomy "${OUTDIR}"/final.taxonomy
mv "${OUTDIR}"/*.0.03.rep.fasta "${OUTDIR}"/final.rep.seqs
mv "${OUTDIR}"/*.0.03.rep.count_table "${OUTDIR}"/final.rep.count_table



####################################
# Make Group-Specific Shared Files #
####################################

# Sample shared file
echo PROGRESS: Creating sample shared file.

# Removing all mock and control groups from shared file leaving only samples
mothur "#remove.groups(shared="${OUTDIR}"/final.shared, groups="${COMBINEDGROUPS}")"

# Renaming output file
mv "${OUTDIR}"/final.0.03.pick.shared "${OUTDIR}"/sample.final.shared



# Mock shared file
echo PROGRESS: Creating mock shared file.

# Removing non-mock groups from shared file
mothur "#get.groups(shared="${OUTDIR}"/final.shared, groups="${MOCKGROUPS}")"

# Renaming output file
mv "${OUTDIR}"/final.0.03.pick.shared "${OUTDIR}"/mock.final.shared



# Control shared file
echo PROGRESS: Creating control shared file.

# Removing any non-control groups from shared file
mothur "#get.groups(shared="${OUTDIR}"/final.shared, groups="${CONTROLGROUPS}")"

# Renaming output file
mv "${OUTDIR}"/final.0.03.pick.shared "${OUTDIR}"/control.final.shared



###############
# Cleaning Up #
###############

echo PROGRESS: Cleaning up working directory.

# Making dir for storing intermediate files (can be deleted later)
mkdir -p "${OUTDIR}"/intermediate/

# Deleting unneccessary files
#rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.fasta")
#rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.map")
#rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.count_table")

# Moving all remaining intermediate files to the intermediate dir
mv $(find "${OUTDIR}"/ -regex ".*\/stability.*") "${OUTDIR}"/intermediate
