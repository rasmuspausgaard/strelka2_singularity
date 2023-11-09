#!/usr/bin/env nextflow


params.normalBam= "/lnx01_data2/shared/testdata/AV1_CRAM/107578340086_AV1_CV6.hg38.V3.BWA.MD.cram"
params.tumorBam = "/lnx01_data2/shared/testdata/AV1_CRAM/107578340086_AV1_CV6.hg38.V3.BWA.MD.cram"
params.runDir = "${launchDir}/NF_Strelkatest_singularity" 
params.referenceFasta = "/data/shared/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
params.genome_fasta = hg19.fa


process RunStrelka {
    input:
    path normalBam from params.normalBam
    path tumorBam from params.tumorBam
    path referenceFasta from params.referenceFasta

    output:
    path "${params.runDir}/*"

    script:
    """
    singularity run -B s_bind /data/shared/programmer/simg/strelka2_2.9.10.sif /data/shared/programmer/simg/tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam          \
    --tumorBam           \
    --referenceFasta     \
    --exome              \
    --runDir 

    singularity run -B ${params.bind_paths} ${params.simg} python2 ${params.runDir}/runWorkflow.py -j 10 -m local
    """
}
