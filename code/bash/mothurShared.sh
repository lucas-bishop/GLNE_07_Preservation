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
export SILVAV4=${2:?ERROR: Need to define SILVAV4.}
export RDPFASTA=${3:?ERROR: Need to define RDPFASTA.}
export RDPTAX=${4:?ERROR: Need to define RDPTAX.}
# export MOCKGROUPS=${5:?ERROR: Need to define MOCKGROUPS.} # List of mock groups in raw data dir separated by '-'
# export CONTROLGROUPS=${6:?ERROR: Need to define CONTROLGROUPS.} # List of control groups in raw data dir separated by '-'

# Other variables
export OUTDIR=data/mothur/process/
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
	classify.otu(list=current, count=current, taxonomy=current, label=0.03)"



# Renaming output files for use later
mv "${OUTDIR}"/*.precluster.pick.pick.fasta "${OUTDIR}"/errorinput.fasta
mv "${OUTDIR}"/*.vsearch.pick.pick.count_table "${OUTDIR}"/errorinput.count_table
mv "${OUTDIR}"/*.opti_mcc.shared "${OUTDIR}"/final.shared
mv "${OUTDIR}"/*.cons.taxonomy "${OUTDIR}"/final.taxonomy


###############
# Cleaning Up #
###############

echo PROGRESS: Cleaning up working directory.

# Making dir for storing intermediate files (can be deleted later)
mkdir -p "${OUTDIR}"/intermediate/

# Deleting unneccessary files
rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.fasta")
rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.map")
rm $(find "${OUTDIR}"/ -regex ".*filter.unique.precluster..*.count_table")

# Moving all remaining intermediate files to the intermediate dir
mv $(find "${OUTDIR}"/ -regex ".*\/stability.*") "${OUTDIR}"/intermediate
