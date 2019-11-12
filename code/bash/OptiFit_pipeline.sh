mothur "#make.contigs(file=$WORKDIR/glne007.files);
	summary.seqs(fasta=current, processors=4);
	screen.seqs(fasta=current, group=current, maxambig=0, maxlength=275);
	summary.seqs(fasta=current);
	unique.seqs(fasta=current);
	count.seqs(name=current, group=current);
	align.seqs(fasta=current, reference=$REF/silva.seed.align);
	summary.seqs(fasta=current, count=current);
	screen.seqs(fasta=current, count=current, start=13862, end=23444, maxhomop=8);
	filter.seqs(fasta=current, vertical=T, trump=.);
	unique.seqs(fasta=current, count=current);
	summary.seqs(fasta=current, count=current);
	pre.cluster(fasta=current, count=current, diffs=2);
	chimera.vsearch(fasta=current, count=current, dereplicate=t);
	remove.seqs(fasta=current, accnos=current);
	classify.seqs(fasta=current, count=current, reference=$REF/trainset14_032015.pds.fasta, taxonomy=$REF/trainset14_032015.pds.tax, cutoff=80);
	remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota);
	remove.groups(fasta=current, count=current, taxonomy=current, groups=mock1-mock2-mock5-mock6-mock7)"

mothur "#unique.seqs(fasta=data/process/baxter/glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, output=count);
dist.seqs(fasta=current, cutoff=0.03);
cluster(column=current, count=current)"

mothur "#cluster.fit(fasta=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, column=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=data/mothur/process/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, reffasta=glne007.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.fasta, refcolumn=old_samples.dist, reflist=old_samples.opti_mcc.list, method=closed)"
