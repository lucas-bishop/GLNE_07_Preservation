# Begum Topcuoglu edited from Snakefile by William L. Close
# Pat Schloss Lab
# University of Michigan

# Snakemake file for running mothur pipeline in GLNE07_Preservation_Bias_XXXX_20XX

readNum = ['R1', 'R2']
mothurSamples = list(set(glob_wildcards(os.path.join('data/mothur/raw/', '{sample}_{readNum, R[12]}_001.fastq.gz')).sample))
mothurMock = ['Mock1-Mock2-Mock3-Mock4']
mothurWater= ['Water1-Water2-Water3-Water4']
mothurBuffer = ['buffer1-buffer2']
mothurControl = ['Neg1-Neg2-Neg3-Neg4']
mothurMockCells = ['MockCells1-MockCells2']

rule all:
	input:
		"data/mothur/process/sample.final.0.03.subsample.shared",
		"data/mothur/process/sample.final.groups.rarefaction",
		"data/mothur/process/sample.final.groups.ave-std.summary",
		"data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.dist",
		"data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.nmds.axes",
		"data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.pcoa.axes",
		"data/mothur/process/error_analysis/errorinput.pick.error.summary"






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
		silvaV4=rules.get16SReferences.output.silvaV4,
		rdpFasta=rules.get16SReferences.output.rdpFasta,
		rdpTax=rules.get16SReferences.output.rdpTax
	output:
		sampleShared="data/mothur/process/sample.final.shared",
		mockShared="data/mothur/process/mock.final.shared",
		controlShared="data/mothur/process/control.final.shared",
		errorFasta="data/mothur/process/errorinput.fasta",
		errorCount="data/mothur/process/errorinput.count_table"
	params:
		mockGroups='-'.join(mothurMock), # Concatenates all mock group names with hyphens
		controlGroups='-'.join(mothurControl) # Concatenates all control group names with hyphens
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} data/mothur/raw/ {params.mockGroups} {params.controlGroups} {input.silvaV4} {input.rdpFasta} {input.rdpTax}"


rule count16SShared:
	input:
		script="code/bash/mothurCountShared.sh",
		sampleShared=rules.make16SShared.output.sampleShared,
		mockShared=rules.make16SShared.output.mockShared,
		controlShared=rules.make16SShared.output.controlShared
	output:
		sampleCount="data/mothur/process/sample.final.count.summary",
		mockCount="data/mothur/process/mock.final.count.summary",
		controlCount="data/mothur/process/control.final.count.summary"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.sampleShared} {input.mockShared} {input.controlShared}"


rule subsample16SShared:
	input:
		script="code/bash/mothurSubsampleShared.sh",
		sampleShared=rules.make16SShared.output.sampleShared,
		mockShared=rules.make16SShared.output.mockShared,
		controlShared=rules.make16SShared.output.controlShared,
		count=rules.count16SShared.output
	output:
		sampleSubsampleShared="data/mothur/process/sample.final.0.03.subsample.shared",
		mockSubsampleShared="data/mothur/process/mock.final.0.03.subsample.shared",
		controlSubsampleShared="data/mothur/process/control.final.0.03.subsample.shared"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.sampleShared} {input.mockShared} {input.controlShared}"





##################################################################
#
# Part 3: Diversity Metrics
#
##################################################################

rule rarefy16SReads:
	input:
		script="code/bash/mothurRarefaction.sh",
		shared=rules.make16SShared.output.sampleShared
	output:
		rarefaction="data/mothur/process/sample.final.groups.rarefaction"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared}"


rule calc16SAlphaDiversity:
	input:
		script="code/bash/mothurAlpha.sh",
		shared=rules.make16SShared.output.sampleShared
	output:
		alpha="data/mothur/process/sample.final.groups.ave-std.summary"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared}"


rule calc16SBetaDiversity:
	input:
		script="code/bash/mothurBeta.sh",
		shared=rules.make16SShared.output.sampleShared
	output:
		sobsDist="data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.dist",
		thetaycDist="data/mothur/process/sample.final.thetayc.0.03.lt.ave.dist",
		braycurtisDist="data/mothur/process/sample.final.braycurtis.0.03.lt.ave.dist"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.shared}"





##################################################################
#
# Part 4: Ordination
#
##################################################################

rule calc16SPCoA:
	input:
		script="code/bash/mothurPCoA.sh",
		dist=rules.calc16SBetaDiversity.output
	output:
		sobsLoadings="data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.pcoa.loadings",
		sobsAxes="data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.pcoa.axes",
		thetaycLoadings="data/mothur/process/sample.final.thetayc.0.03.lt.ave.pcoa.loadings",
		thetaycAxes="data/mothur/process/sample.final.thetayc.0.03.lt.ave.pcoa.axes",
		braycurtisLoadings="data/mothur/process/sample.final.braycurtis.0.03.lt.ave.pcoa.loadings",
		braycurtisAxes="data/mothur/process/sample.final.braycurtis.0.03.lt.ave.pcoa.axes"
	conda:
		"envs/mothur.yaml"
	shell:
		"bash {input.script} {input.dist}"


rule calc16SNMDS:
	input:
		script="code/bash/mothurNMDS.sh",
		dist=rules.calc16SBetaDiversity.output
	output:
		sobsStress="data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.nmds.stress",
		sobsAxes="data/mothur/process/sample.final.sharedsobs.0.03.lt.ave.nmds.axes",
		thetaycStress="data/mothur/process/sample.final.thetayc.0.03.lt.ave.nmds.stress",
		thetaycAxes="data/mothur/process/sample.final.thetayc.0.03.lt.ave.nmds.axes",
		braycurtisStress="data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.stress",
		braycurtisAxes="data/mothur/process/sample.final.braycurtis.0.03.lt.ave.nmds.axes"
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
