#!/usr/bin/env nextflow

params.normal_cram = "/path/to/default/normal.cram"
params.tumor_cram = "/path/to/default/tumor.cram"
params.runDir = "${launchDir}/NF_Strelkatest_singularity"
params.simg = "/path/to/strelka2.sif"
params.referenceFasta = "/path/to/default/reference.fa"
params.genome_fasta = params.referenceFasta // Assign referenceFasta to genome_fasta if that's the intended use

process RunStrelka {
    // Other directives...
    
    input:
    path normal_cram from params.normal_cram
    path tumor_cram from params.tumor_cram
    path genome_fasta from params.genome_fasta // Now genome_fasta should be correctly defined

    // Output and script directives...
}
