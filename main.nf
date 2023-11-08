#!/usr/bin/env nextflow
user="$USER"
date=new Date().format( 'yyMMdd' )



params.normal_cram        = "/lnx01_data2/shared/testdata/AV1_CRAM/107578340086_AV1_CV6.hg38.V3.BWA.MD.cram"
params.tumor_cram         = "/lnx01_data2/shared/testdata/AV1_CRAM/107578340086_AV1_CV6.hg38.V3.BWA.MD.cram"
params.runDir             = "${launchDir.baseName} /NF_Strelkatest_singularity"
params.simg               = "/data/shared/programmer/simg/strelka2_2.9.10.sif"
params.referenceFasta     = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
params.genome_fasta       = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"


process RunStrelka {
    bind_paths = "/data/shared/programmer/simg/"

    input:
    path normal_cram from params.normal_cram
    path tumor_cram from params.tumor_cram
    path genome_fasta from params.referenceFasta

    output:
    path "${params.runDir}/*"

    script:
    """
    singularity run -B ${bind_paths} ${params.simg} /tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_cram} \
    --tumorBam ${tumor_cram} \
    --referenceFasta \
    --exome \
    --runDir ${params.runDir}

    singularity run -B ${bind_paths} ${params.simg} python2 ${params.runDir}/runWorkflow.py -j 10 -m local
    """
}
