# Begum Topcuoglu edited from Snakefile by William L. Close
# Pat Schloss Lab
# University of Michigan

# Snakemake file for running mothur pipeline in GLNE07_Preservation_Bias_XXXX_20XX

# NOTE: Change these settings before running workflow
mothurMock = ['Mock1', 'Mock2', 'Mock3', 'Mock4', 'Mock5', 'Mock6', 'Mock7', 'Mock8']

mothurControl = ['Neg1', 'Neg2', 'Neg3', 'Neg4', 'Neg5', 'Neg6', 'Neg7', 'Neg8']

mothurBuffer = ['Water5a', 'Water5b', 'Water6a', 'Water6b', 'Water7a', 'Water7b', 'bufferP1.1', 'bufferP1.2', 'bufferP2.1', 'bufferP2.1b', 'bufferP2.2', 'bufferP2.2b', 'bufferP2.3', 'bufferP2.4', 'bufferP2.5']

mothurAlpha = ['nseqs','coverage','invsimpson','shannon','sobs']

mothurBeta = ['sharedsobs','thetayc','braycurtis']

# Leave these settings as is
readNum = ['R1', 'R2']
mothurSamples = list(set(glob_wildcards(os.path.join('/nfs/turbo/schloss-lab/begumtop/glne_sequences/', '{sample}_{readNum, R[12]}_001.fastq.gz')).sample))

mothurGroups = ['sample','mock','control', 'buffer']



rule all:
	input:
		expand("data/mothur/process/{group}.final.count.summary",
			group = mothurGroups),
		expand("data/mothur/process/{group}.final.0.03.subsample.shared",
			group = ['sample','mock']),
		"data/mothur/process/sample.final.groups.rarefaction",
		"data/mothur/process/sample.final.groups.ave-std.summary",
		expand("data/mothur/process/sample.final.{beta}.0.03.lt.ave.dist",
			beta = mothurBeta),
		expand("data/mothur/process/sample.final.{beta}.0.03.lt.ave.nmds.axes",
			beta = mothurBeta),
		expand("data/mothur/process/sample.final.{beta}.0.03.lt.ave.pcoa.axes",
			beta = mothurBeta),
		"data/mothur/process/error_analysis/errorinput.pick.error.summary"
	shell:
		"""
		rm data/mothur/process/*rabund
		mkdir -p logs/mothur/
		mv mothur*logfile logs/mothur/
		"""





##################################################################
#
# Part 1: Generate Reference and Mock Control Files
#
##################################################################

rule get16SReferences:
	input:
		script="code/bash/mothurReferences.sh"
	output:
		silvaV4="data/mothur/references/silva.v4.align",
		rdpFasta="data/mothur/references/trainset16_022016.pds.fasta",
		rdpTax="data/mothur/references/trainset16_022016.pds.tax"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script}"


rule get16SMock:
	input:
		script="code/bash/mothurMock.sh",
		silvaV4=rules.get16SReferences.output.silvaV4
	output:
		mockV4="data/mothur/references/zymo.mock.16S.v4.fasta"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script}"





##################################################################
#
# Part 2: Generate Shared Files
#
##################################################################

rule make16SShared:
	input:
		script="code/bash/mothurShared.sh",
		raw=expand('data/mothur/raw/{mothurSamples}_{readNum}_001.fastq.gz',
			mothurSamples = mothurSamples, readNum = readNum),
		refs=rules.get16SReferences.output
	output:
		shared="data/mothur/process/final.shared",
		taxonomy="data/mothur/process/final.taxonomy",
		errorFasta="data/mothur/process/errorinput.fasta",
		errorCount="data/mothur/process/errorinput.count_table"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} data/mothur/raw/ {input.refs}"


rule split16SShared:
	input:
		script="code/bash/mothurSplitShared.sh",
		shared=rules.make16SShared.output.shared
	output:
		shared=expand("data/mothur/process/{group}.final.shared",
			group = mothurGroups)
	params:
		mockGroups='-'.join(mothurMock), # Concatenates all mock group names with hyphens
		controlGroups='-'.join(mothurControl) # Concatenates all control group names with hyphens
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {params.mockGroups} {params.controlGroups} {params.bufferGroups}"


rule count16SShared:
	input:
		script="code/bash/mothurCountShared.sh",
		shared="data/mothur/process/{group}.final.shared"
	output:
		count="data/mothur/process/{group}.final.count.summary"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared}"


rule subsample16SShared:
	input:
		script="code/bash/mothurSubsampleShared.sh",
		shared="data/mothur/process/{group}.final.shared",
		count="data/mothur/process/{group}.final.count.summary"
	output:
		subsampleShared="data/mothur/process/{group}.final.0.03.subsample.shared"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared} {input.count}"





##################################################################
#
# Part 3: Diversity Metrics
#
##################################################################

rule rarefy16SReads:
	input:
		script="code/bash/mothurRarefaction.sh",
		shared="data/mothur/process/sample.final.shared"
	output:
		rarefaction="data/mothur/process/sample.final.groups.rarefaction"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared}"


rule calc16SAlphaDiversity:
	input:
		script="code/bash/mothurAlpha.sh",
		shared="data/mothur/process/sample.final.shared",
		count="data/mothur/process/sample.final.count.summary"
	output:
		alpha="data/mothur/process/sample.final.groups.ave-std.summary"
	params:
		alpha='-'.join(mothurAlpha) # Concatenates all alpha metric names with hyphens
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared} {input.count} {params.alpha}"


rule calc16SBetaDiversity:
	input:
		script="code/bash/mothurBeta.sh",
		shared="data/mothur/process/sample.final.shared",
		count="data/mothur/process/sample.final.count.summary"
	output:
		dist=expand("data/mothur/process/sample.final.{beta}.0.03.lt.ave.dist",
			beta = mothurBeta)
	params:
		beta='-'.join(mothurBeta) # Concatenates all beta metric names with hyphens
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared} {input.count} {params.beta}"





##################################################################
#
# Part 4: Ordination
#
##################################################################

rule calc16SPCoA:
	input:
		script="code/bash/mothurPCoA.sh",
		dist="data/mothur/process/sample.final.{beta}.0.03.lt.ave.dist"
	output:
		loadings="data/mothur/process/sample.final.{beta}.0.03.lt.ave.pcoa.loadings",
		axes="data/mothur/process/sample.final.{beta}.0.03.lt.ave.pcoa.axes"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.dist}"


rule calc16SNMDS:
	input:
		script="code/bash/mothurNMDS.sh",
		dist="data/mothur/process/sample.final.{beta}.0.03.lt.ave.dist"
	output:
		stress="data/mothur/process/sample.final.{beta}.0.03.lt.ave.nmds.stress",
		axes="data/mothur/process/sample.final.{beta}.0.03.lt.ave.nmds.axes"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.dist}"





##################################################################
#
# Part 5: Quality Control
#
##################################################################

rule calc16SError:
	input:
		script="code/bash/mothurError.sh",
		errorFasta=rules.make16SShared.output.errorFasta,
		errorCount=rules.make16SShared.output.errorCount,
		mockV4=rules.get16SMock.output.mockV4
	output:
		summary="data/mothur/process/error_analysis/errorinput.pick.error.summary"
	params:
		mockGroups='-'.join(mothurMock) # Concatenates all mock group names with hyphens
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.errorFasta} {input.errorCount} {input.mockV4} {params.mockGroups}"
