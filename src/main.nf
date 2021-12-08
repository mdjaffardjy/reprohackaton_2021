#!/usr/bin/env nextflow
// ====================================================== Main workflow for RNA-seq data analysis ======================================================= // 

channel_SRP = Channel.from('SRR628589') 

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
