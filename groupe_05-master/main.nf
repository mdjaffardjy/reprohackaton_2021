#!/usr/bin/env nextflow
// ====================================================== Main workflow for RNA-seq data analysis ======================================================= // 

channel_SRP = Channel.from('SRR628589', 'SRR628588', 'SRR628587', 'SRR628586', 'SRR628585', 'SRR628584', 'SRR628583', 'SRR628582') 

process downloadSRR {
    /*
    Input : 8 samples from patients with uveal melanoma from a SRProject on SRArchives (NCBI, API key provided in .config file)
    Output : 8 samples in SRA format and SRR code of samples (i.e. SRR6285xx for xx in 82 .. 89), sra identifiers and files
    Function : Download files
    */
    publishDir "results/SRR/"

    input:
    val srr_id from channel_SRP

    output:
    tuple val(srr_id), file("${srr_id}.sra") into channel_srr 
    
    script:
    """
    wget -O ${srr_id}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${srr_id}/${srr_id}.1
    """
}


process fastqConversion {
    /*
    Input  : Data from downloadSRR output channel
    Output : 8 samples in FASTQ format. Since sequencing is paired-end, two FASTQ files are generated from one SRA sample. 
    	     We also keep in output channel SRA Ids. Every emission in the channel gives 3 values.
    Function : Converts SRA format to FASTQ format, without dezipping
    */
    publishDir "results/FastQ/"
    
    input:
    tuple val(sraid), file(srr) from channel_srr
    
    output:
    tuple val(sraid), file("*1.fastq.gz"), file("*2.fastq.gz") into fastqFiles
    
    script:
    """    
    fastq-dump --gzip --split-files ${srr}
    """
}


//chromosomes channel
channel_chr = Channel.from('1','2' ,'3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','MT','X','Y')

process downloadHumanDNA {
    /*
    Input : Values corresponding to each chromosome from the chromosomes list "channel_chr"
    Output : A set of fasta files with each one corresponding to a chromosome in gz format transiting in the output channel channel_chromosome
    Function : Downloads the fasta files of the human DNA for 25 chromosomes through the link below
    */
    publishDir "results/genome/"

    input:
    val chr from channel_chr

    output: 
    file "${chr}.fa.gz" into channel_chromosome
    
    script:
    """
    wget -O ${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    """
}


process concatenationGenome {
    /*
    Input : Collection of all the downloaded chromosomes in the "downloadHumanDNA" process coming from the chanell and sending them all at once to the 
   	    same output channel
    Output : A fasta file grouping all the chromosomes forming the entire Human genome
    Function : Decompression and concatenation of all the fa.gz files associated with the different chromosomes merged in a single fasta file
              corresponding to the Human Genome
    */
    publishDir "results/genome/"

    input:
    file all_chr from channel_chromosome.collect()
		
    output: 
    file "HumanGenome.fa" into channel_human
    
    script:
    """
    gunzip -c ${all_chr} > HumanGenome.fa
    """
}


process genomeIndexation {
    /*
    Input : Fasta file containing the Human genome
    Output : A genome index that will be used during the mapping 
    Function : Generating genome index files from the fasta file containing the Human genome.
    */
    publishDir "results/GenomeIndex/"
    
    input: 
    file (genome) from channel_human.collect()

    output:
    path "ref" into channel_genome_index

    script:
    """
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
    """
}
  
  
process downloadGTF {
    /*
    Output : Human genome annotation in .gtf format
    Function : Downloading and decompressing the Human genome annotation available through the link below
    */
    publishDir "results/"
    
    output:
    file "human.gtf" into channel_gtf

    script:
    """
    wget -O human.gtf.gz ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz
    gunzip -c human.gtf.gz > human.gtf 
    """
}


process mappingFastQFiles {
   /* 
   Input : SRA ids and fastq files
   Output : .bam files 
   Function :  RNA-Seq reads mapping on human genome
   Details : Arguments used
   --outSAMstrandField intronMotif : "For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute"
   --outFilterMismatchNmax : "alignment will be output only if it has no more mismatches than the selected value"
   --outFilterMultimapNmax : "max number of multiple alignments allowed for a read"
   --genomeDir : path to genome index
   --readFilesIn : fastq files
   --runThreadN : number of CPU needed
   --outSAMunmapped : "output of unmapped reads in the SAM format"
   --outSAMtype : "type of signal output" 
   --outStd : "which output will be directed to stdout"
   --genomeLoad : "mode of shared memory usage for the genome files"
   --limitBAMsortRAM : "maximum available RAM for BAM sorting"
   --outFileNamePrefix : output files name prefix  
   */
    publishDir "results/Mapping/"
    
    input:
    tuple val(sraid), file(r1), file(r2) from fastqFiles
    path index from channel_genome_index
    
    output:
    file "${sraid}.bam" into bam_channel, bam_channel2
    
    script:
    """
    STAR --outSAMstrandField intronMotif \
    --outFilterMismatchNmax 4 \
    --outFilterMultimapNmax 10 \
    --genomeDir ${index} \
    --readFilesIn <(gunzip -c ${r1}) <(gunzip -c ${r2}) \
    --runThreadN ${task.cpus} \
    --outSAMunmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --outStd BAM_SortedByCoordinate \
    --genomeLoad NoSharedMemory \
    --limitBAMsortRAM ${task.memory.toBytes()} \
    > ${sraid}.bam
    """
}


process bamIndex {
    /*
    Input : Bam files from the previous process
    Output : Indexed bam files 
    Function : Indexing bam files
    */
    publishDir "results/bamIndex/"
    
    input: 
    file bam from bam_channel
    
    output:
    file "${bam}.bai" into channel_bam_index
    
    script:
    """
    samtools index $bam
    """
}
 
 
process countingReads {
    /*
    Input : Collecting bam files and their indexes from the previous channels and the Human genome annotation
    Output : a large matrix in .txt format where the columns represent the samples and the rows correspond to the genes
    Function : Counting the aligned reads on each gene during the mapping step
    Details : Function parameters :
    	-T : threads or cpus number
	-t : feature type to count 
	-g : attribute type in GTF annotation, the gene_id available in the annotation file
	-s : unstranded read counting
	-a : Human genome annotation file 
	-o : name of the output file including read counts
    */
    publishDir "results/featureCounts/"
    
    input:
    file bam_index from channel_bam_index.collect()
    file bam from bam_channel2.collect()
    val gtf from channel_gtf 
    
    output:
    file "matrice_featureCounts.txt" into channel_featureCounts
    
    script: 
    """
    featureCounts -T ${task.cpus} -t gene -g gene_id -s 0 -a $gtf -o matrice_featureCounts.txt ${bam}
    """
}


process exon_count{
        input:
        	file bam_index from channel_bam_index.collect()
    file bam from bam_channel2.collect()
    val gtf from channel_gtf
        output:
                file "matrice_exonsCounts.txt" into channel_exonsCounts
        shell:
                "featureCounts -p -T ${task.cpus} -t exon -g exon_id -s 0 -a ${gtf} -o matrice_exonsCounts.txt ${bam}"

}

metaData_channel = Channel.fromPath('./metaData.txt')

process statAnalysis {
    /*
    Input : Matrix from countingReads channel output and associated meta-sata from files and statAnalysis.R script
    Output : PDF files containing outputs graphics summarizing statistic analyses and TXT containing DESeq2 results
    Function : DESeq2 analysis and PCA : See statAnalysis.R for more details
    */
    publishDir "results/statResults/"
    
    input:
    file countData from channel_featureCounts
    file countExon from channel_exonsCounts
    file metaData from metaData_channel
    
    output:
    tuple file("DESeq2_results.txt"), file("analysis_result*") into channel_end
    
    script:
    """
    Rscript ${workflow.projectDir}/bin/statAnalysis.R ${metaData} ${countData} "analysis_result" 
    """
}


workflow.onComplete = {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
} 

 
workflow.onError = {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
 
