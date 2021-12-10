#Project HACKATHON - reproducible pipeline
#Call command: snakemake -s pipeline.smk --cores 2 --use-singularity

### Beware : this pipeline needs to have the sample and chromosome names as empty file to work.
### > Choose appropriate launcher before running.
### Beware2 : be ready to press 'x' after rule configure ran (only interaction required in blue panel of SRA configuration) 

sample = ['SRR628585', 'SRR628586', 'SRR628587', 'SRR628588', 'SRR628582', 'SRR628584', 'SRR628583', 'SRR628589'] #please, do not change the order, it is essential for Deseq2
knames = ['1', '2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','MT']

rule all:
	""" When does the pipeline stop ? Here are the final output files we want."""
	input:
		"WT_vs_mutated.csv",
		"WT_vs_mutated_annotated.csv",
		"Rplots.pdf",
		"Rplots1.pdf",
		"Rplots2.pdf",
		"Rplots3.pdf"

rule SRA:
	""" To download runs from SRA - parallelized """
	input: "{sample}"
	output: "{sample}.sra"
	shell: "wget -O {input}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/{input}/{input}.1"

rule configure:
	"""Configuration of the docker"""
	output: "configDone.txt" #indicate when configuration has been done
	singularity: "docker://evolbioinfo/sratoolkit:v2.10.8"
	shell:
		"""
		vdb-config --interactive
		touch {output} 
		"""
#BEWARE : user will need to press the "x" key once because of the forced interaction to configure
#(or ctrl+x if the sole "x" key does not suffice)

rule Fastq:
	""" To convert the runs in fastq file" - parallelized """
	input: 
		config="configDone.txt", #make sure configuration of the docker has been done
		ech="{sample}.sra"
	output: "{sample}_1.fastq.gz","{sample}_2.fastq.gz"
	singularity: "docker://evolbioinfo/sratoolkit:v2.10.8"
	shell: "fastq-dump --gzip --split-files ./{input.ech}"

rule Kseq:
	""" To download all human chromosomes + mitochondrial DNA, needed to build the map (quiet mode) - parallelised"""
	input: "{knames}"
	output: "{knames}.fa.gz"
	shell: "wget -q -O {output} ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.{output}"

rule regroupKseq:
	"""To group up the chromosomes sequences in a reference genome"""
	input: expand("{k}.fa.gz",k=knames) #need ALL files before starting to regroup
	output: "ref.fa"
	shell: "gunzip -c *.fa.gz > {output}"

rule GTF:
	""" Downloading genome annotations for indexing (quiet mode) """
	output: "chr.gtf.gz"
	shell: "wget -q -O {output} ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz"

rule index:
	""" To create the index from human genome in a directory named ref """
	input:
		ref = "ref.fa",
		gtf = "chr.gtf.gz" 
	output: marker="indexDone.txt" #a file indicating when index is done. Because a directory type as input (for rule mapping) is not valid
	singularity: "docker://evolbioinfo/star:v2.7.6a"
	shell:
		"""
		gunzip -c {input.gtf} > chr.gtf
		STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ref --genomeFastaFiles {input.ref} --sjdbGTFfile chr.gtf 
		touch {output.marker}
		"""

#fastqBig needed because of issues with mapping otherwise (although it is not optimized in space...)
rule fastqBig:
	""" Unzipping fastq file. while keeping original file - parallelised """
	input:
		f1="{sample}_1.fastq.gz",
		f2="{sample}_2.fastq.gz"
	output:
		"{sample}_1.fastq",
		"{sample}_2.fastq"
	shell: 
		"""
		gunzip -k {input.f1}
		gunzip -k {input.f2}
		"""
rule mapping:
	""" To map the reads onto the index - parallelised """
	input:
		indexation="indexDone.txt", #Make sure index has been FULLY done before starting to map
		read1="{sample}_1.fastq",
		read2="{sample}_2.fastq",
		nameFile="{sample}" #Use empty file to create varying tmp directory name
	output: "{sample}.bam"
	singularity: "docker://evolbioinfo/star:v2.7.6a"
	shell:
		"""STAR --outSAMstrandField intronMotif \
		--outFilterMismatchNmax 4 \
		--outFilterMultimapNmax 10 \
		--genomeDir ref \
		--readFilesIn {input.read1} {input.read2} \
		--runThreadN 10 \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outStd BAM_SortedByCoordinate \
		--genomeLoad NoSharedMemory \
		--limitBAMsortRAM 0 \
		--outTmpDir ./{input.nameFile}_tmp > {output} """

rule indexBam:
	""" To index the bam file - parallelised """
	input: "{sample}.bam"
	output: "{sample}.bam.bai" #needed if one wants to visualise content of BAM file on IGV
	singularity: "docker://evolbioinfo/samtools:v1.11"
	shell: "samtools index {input}"

rule featureCounts:
	""" To count reads """
	input:
		expand("{sample}.bam", sample=sample), #Need ALL files before starting
		expand("{sample}.bam.bai",sample=sample), #indicated .bai files here. Otherwise, would have been useless to generate
		"chr.gtf.gz"
	output: "counts.txt"
	singularity: "docker://evolbioinfo/subread:v2.0.1"
	shell:
		"""
		gunzip -c chr.gtf.gz > chr.gtf
		featureCounts -T 10 -t gene -g gene_id -s 0 -a chr.gtf -o ./counts.txt ./*.bam
		"""
		# strandness should be set to 0 (-s 0), tested on 0, 1 and 2 by Anna

rule deseq2:
	"""Differential expression analysis"""
	input:
		"counts.txt",
		"deseq2_annot_visualisation.R"
	output:
		"WT_vs_mutated.csv",
		"WT_vs_mutated_annotated.csv",
		"Rplots.pdf",
		"Rplots1.pdf",
		"Rplots2.pdf",
		"Rplots3.pdf"
	singularity: "docker://evolbioinfo/deseq2:v1.28.1"
	shell: "Rscript deseq2_annot_visualisation.R"

#End