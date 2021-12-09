#!/usr/bin/env nextflow

// Ce workflow a pour objectif de réaliser une analyse de séquences "RNA-seq", afin d'étudier l'expression différentielle de l'ensemble des gènes de nos échantillons. Il est écrit pour effectuer cette analyse sur des données brutes "SRA" spécifiques, disponibles sur un dépôt du ncbi (lien :https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP017413&o=acc_s%3Aa).


/////////////////////////////////////////////////////////////////////
// DEFINITION DES CHANNELS                                         //
/////////////////////////////////////////////////////////////////////


// Identifiants SRA (SRAID) des echantillons
SRAID = Channel.from("SRR628582","SRR628583", "SRR628584", "SRR628585", "SRR628586", "SRR628587", "SRR628588", "SRR628589")

// Chromosomes d'interet
liste_chromosomes = Channel.from("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16" , "17", "18", "19", "20", "21", "22", "MT", "X", "Y")

gtf_URL="ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.chr.gtf.gz"

/////////////////////////////////////////////////////////////////////
// DEFINITION DES PROCESS                                          //
/////////////////////////////////////////////////////////////////////

// Téléchargement des données de séquencage brutes .sra
process download_SRA{
	publishDir "data/SRA/"

	input:
	val sraid from SRAID

	output:
	tuple val(sraid), file("${sraid}.sra") into SRAFiles

	
	"""
	wget -O ${sraid}.sra https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos1/sra-pub-run-5/${sraid}/${sraid}.1
	"""
}

// Conversion des SRA en fastq
process conversion_FASTQ {
    publishDir "data/FastQ/"
    
    input:
    tuple val(sraid), file("${sraid}.sra") from SRAFiles
    
    output:
    tuple val(sraid), file("*1.fastq.gz"), file("*2.fastq.gz") into FASTQ_files, FASTQ_files_qc
    
    
    """    
    fastq-dump --gzip --split-files ${sraid}.sra
    """
}


// Téléchargement du génome de référence -- séquences de chaque chromosome pour procéder à l'alignement - version GRCh38
process get_chr_seq {
	publishDir "data/Genome/"

    input:
    val chr from liste_chromosomes
    
    output:
    file "*.fa.gz" into chr_fa
    
    
    """
    wget -O chr${chr}.fa.gz ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.${chr}.fa.gz
    """
}

// Contrôle qualité des reads
process fastqc {
    publishDir "results/fastqc_results/", mode: "copy"

    input:
    tuple val(sraid), file(read_fw), file(read_rc) from FASTQ_files_qc

    output:
    tuple file("${SRAID}_1_fastqc.html"), file("${SRAID}_2_fastqc.html")

    script:
    """
    fastqc ${read_fw} 
    fastqc ${read_rc}
    """

}


process concatenation_genome {
    publishDir "data/genome/"

    input:
    file '*.fa.gz' from chr_fa.collect()
		
    output: 
    file "HumanGenome.fa" into genome_hum
    
    script:
    """
    gunzip -c *.fa.gz > HumanGenome.fa
    """
}

// Alignement avec STAR. Il faut tout d'abord créer un index :
process make_STAR_index {
    publishDir "results/GenomeIndex/"
    
    input:
    file (genome) from genome_hum.collect()

    output:
    path "ref" into genome_idx

    script:
    """
    STAR --runThreadN ${task.cpus} --runMode genomeGenerate --genomeDir ref/ --genomeFastaFiles ${genome}
    """
}

// Obtention de fichier d'annotation du génome sous le format gff (.gtf) pour compter les reads
process download_genome_annotations{
    publishDir "data/GenomeAnnotation/"
    
    input:
    val gtfURL from gtf_URL
    
    output:
    file "annot.gtf" into gtf
    
    """
    wget ${gtfURL}
    gunzip -c *.gtf.gz > annot.gtf
    """
}

// Mapping des fichiers FASTQ grâce à l'index précédemment créé avec STAR
process map_FASTQ {
    publishDir "data/bam/"

    input:
    tuple val(sar), file(fq_fw), file(fq_rv) from FASTQ_files
    path ref from genome_idx
 
    output:
    file "${sar}.bam" into mapped_fq_1, mapped_fq_2
 
    script:
    """
    STAR --outSAMstrandField intronMotif \
        --outFilterMismatchNmax 4 \
        --outFilterMultimapNmax 10 \
        --genomeDir ${ref} \
        --readFilesIn <(gunzip -c ${fq_fw}) <(gunzip -c ${fq_rv}) \
        --runThreadN ${task.cpus} \
        --outSAMunmapped None \
        --outSAMtype BAM SortedByCoordinate \
        --outStd BAM_SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --limitBAMsortRAM ${task.memory.toBytes()} \
        --outFileNamePrefix ${sar} 
        >${sar}.bam    
    """
}

// Indexation des fichiers bam
process index_BAM {
    publishDir "data/bam_index/"

    input:
    file bam from mapped_fq_1
 
    output:
    file "*.bai" into bam_index
 
    script:
    """
    samtools index ${bam}
    """
}

// Comptage des reads avec featureCounts
process count_reads{
    publishDir "results/counts/"
    
    input:
    file gtf_file from gtf
    file bam from mapped_fq_2.collect()
    
    output:
    file "counts.txt" into file_count
    file "output.counts.summary" into logsFileCount
    
    script:
    """ 
     featureCounts ${bam} -p -T ${task.cpus} -t gene -g gene_id -s 0 -a ${gtf} -o counts.txt 
    """
} 


// Analyse de l'expression différentielle avec un script R et génération des graphes et tableaux

process counts_analysis{
	
	publishDir "results/analysis/"
	
	input:
	file count from file_count	
	
	output:
	tuple file("results.txt"), file("*.csv"), file("analyse*") into chann_end
	
	script:
    """
    Rscript ${workflow.projectDir}/bin/diff_expr_analysis.R ${count}
    """
}



workflow.onComplete = {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
} 


workflow.onError = {
    println "The pipeline had trouble executing : ${workflow.errorMessage}"
}











































