#!/usr/bin/env nextflow

date=new Date().format( 'yyMMdd' )
user="$USER"
runID="${date}.${user}"
params.help = null
params.fastq = null


///////////////////////////////USE///////////////////////////////
// singularity run  data/shared/programmer/simg/strelka2_2.9.10.sif 
// cd .. indtil man er i en mappe hvor man kan køre ls og se tools
// tools/strelka2/bin/configureStrelkaSomaticWorkflow.py --normalBam ../"cramfilen med sti" --tumorBam ../"cramfilen med sti" 
// "/lnx01_data2/shared/testdata/AV1_CRAM/107578340086_AV1_CV6.hg38.V3.BWA.MD.cram"


// unset parameters
params.normalCram           =null 
params.normalCrai           =null
params.tumorCram            =null 
params.tumorCrai            =null

// preset parameters
params.hg38v1               =null  // primary hg38 full, original hg38 assembly, no decoys, etc.
params.hg38v2               =null  // UCSC NGS set
params.hg38v3               =null  // DEFAULT: NGC version, masked regions. 
params.cnvpytorHis1=200
params.cnvpytorHis2=1000
params.cnvpytorHis3=10000

params.genome               ="hg38"
params.outdir               ="${launchDir.baseName}.NF_Strelkatest_singularity"
params.runDir               ="${launchDir.baseName}"

params.server               ="lnx01"


switch (params.server) {
    case 'lnx01':
        syspath="/data/shared";
        s_bind="/data/:/data/,/lnx01_data2/:/lnx01_data2/";
        simgpath="/data/shared/programmer/simg";
        tmpDIR="/data/TMP/TMP.${user}/";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        multiqc_config="/data/shared/programmer/configfiles/multiqc_config.yaml";
        //tank_storage="/home/mmaj/tank.kga2/data/data.storage.archive";
        tank_storage="/lnx01_data2/shared/dataArchive";
        genomes_dir="/lnx01_data2/shared/dataArchive/tank_kga_external_archive/genomes";
    break;
    case 'kga01':
        simgpath="/data/shared/programmer/simg";
        syspath="/data/shared";
        s_bind="/data/:/data/";
        tmpDIR="/data/TMP/TMP.${user}/";
        params.intervals_list="/data/shared/genomes/hg38/interval.files/WGS_splitIntervals/wgs_splitinterval_BWI_subdivision3/*.interval_list";
        gatk_exec="singularity run -B ${s_bind} ${simgpath}/gatk4261.sif gatk";
        //tank_storage="/home/mmaj/tank.kga/data/data.storage.archive";
        tank_storage="/data/shared/dataArchive";
       // genomes_dir="/home/mmaj/tank.kga/data/data.storage.archive/genomes";
        genomes_dir="/data/shared/dataArchive/tank_kga_external_archive/genomes";
    break;
}


switch (params.genome) {
    case 'hg19':
        assembly="hg19"
        // Genome assembly files:
        genome_fasta = "${syspath}/genomes/hg19/human_g1k_v37.fasta"
        genome_fasta_fai = "${syspath}/genomes/hg19/human_g1k_v37.fasta.fai"
        genome_fasta_dict = "${syspath}/genomes/hg19/human_g1k_v37.dict"
        genome_version="V1"
        sve_genome = "${syspath}/genomes/hg19/SVE/SVE.human_g1k_v37.fa"
        // Gene and transcript annotation files:

        gencode_gtf = "${syspath}/genomes/GRCh37/gene.annotations/gencode.v19.annotation.gtf"
        
         //Program  files:
        msisensor_list="${syspath}/genomes/hg19/human_g1k_v37.microsatellites.list"

        cnvradar_anno="${syspath}/genomes/hg19/databases/vcfs/All_20180423.vcf.gz"
        
        cnvradar_anno_idx="${syspath}/genomes/hg19/databases/vcfs/All_20180423.vcf.gz.tbi"
        
        cnvradar_ROI="${syspath}/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.sorted.bed"
        
        cnvradar_roisum_dir="${syspath}/genomes/hg19/databases/cnvradar_roi_summaries/"
        expansionhunter_db="${syspath}/programmer/ExpansionHunter-v4.0.2-linux_x86_64/variant_catalog/grch37/variant_catalog.json"
        // Somatic calling files GATK Mutect2 pipeline
        gatk_wgs_pon="${syspath}/genomes/hg19/databases/vcfs/somatic-b37_Mutect2-WGS-panel-b37.vcf"
        mutect_gnomad="${syspath}/genomes/hg19/databases/vcfs/somatic-b37_af-only-gnomad.raw.sites.vcf"
        gatk_contamination_ref="${syspath}/genomes/hg19/databases/vcfs/somatic-b37_small_exac_common_3.vcf"
       
        // Program indexes
        pcgr_assembly="grch37"
        sequenza_cg50_wig="${syspath}/genomes/hg19/human_g1k_v37.cg50base.wig.gz"

        // Regions / variants:
        dbsnp= "${syspath}/genomes/hg19/databases/vcfs/dbsnp147_All_20160601.vcf"
        qualimap_ROI="${syspath}/genomes/hg19/interval.files/200108.NCBIrefseq.codingexons.nocontig.20bp.merged.sorted.6col.bed"
        
        ROI="${syspath}d/genomes/hg19/interval.files/WES/IDT.exomes.EV7/EV7.ROI.bed"

        excluderegions="${syspath}/genomes/hg19/interval.files/WGS/gaplist.b37.cleaned.bed"

        KGindels="${syspath}/genomes/hg19/databases/vcfs/1000G_phase1.indels.b37.vcf"
        KGindels_idx="${syspath}/genomes/hg19/databases/vcfs/1000G_phase1.indels.b37.vcf.idx"
        KGmills="${syspath}/genomes/hg19/databases/vcfs/Mills_and_1000G_gold_standard.indels.b37.vcf"
        KGmills_idx="${syspath}/genomes/hg19/databases/vcfs/Mills_and_1000G_gold_standard.indels.b37.vcf.idx"

        KG_p1_High_snps="${syspath}/genomes/hg19/databases/vcfs/1000G_phase1.snps.high_confidence.b37.vcf"

        hapmap="${syspath}/genomes/hg19/databases/vcfs/hapmap_3.3.b37.vcf"
        omni="${syspath}/genomes/hg19/databases/vcfs/1000G_omni2.5.b37.vcf"

        break;


    case 'hg38':
        assembly="hg38"
        smncaller_assembly="38"
        // Genome assembly files:
        if (params.hg38v1) {
        genome_fasta = "${syspath}/genomes/hg38/GRCh38.primary.fa"
        genome_fasta_fai = "${syspath}/genomes/hg38/GRCh38.primary.fa.fai"
        genome_fasta_dict = "${syspath}/genomes/hg38/GRCh38.primary.dict"
        genome_version="V1"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_germline_PON/jgmr_45samples.reference.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v1_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="${syspath}/genomes/hg38/inhouse_DBs/hg38v1_primary/"
        }
        
        if (params.hg38v2){
        genome_fasta = "${syspath}/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa"
        genome_fasta_fai = "${syspath}/genomes/hg38/ucsc.hg38.NGS.analysisSet.fa.fai"
        genome_fasta_dict = "${syspath}/genomes/hg38/ucsc.hg38.NGS.analysisSet.dict"
        genome_version="V2"
        }

        // Current hg38 version (v3): NGC with masks and decoys.
        if (!params.hg38v2 && !params.hg38v1){
        genome_fasta = "${syspath}/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa"
        genome_fasta_fai = "${syspath}/genomes/hg38/GRCh38_masked_v2_decoy_exclude.fa.fai"
        genome_fasta_dict = "${syspath}/genomes/hg38/GRCh38_masked_v2_decoy_exclude.dict"
        genome_version="V3"
        //cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/hg38v3_109samples.cnvkit.reference.cnn"
        cnvkit_germline_reference_PON="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/WGS_109samples_hg38v3_cnvkitRef.cnn"
        cnvkit_inhouse_cnn_dir="/data/shared/genomes/hg38/inhouse_DBs/hg38v3_primary/cnvkit/wgs_persample_cnn/"
        inhouse_SV="${syspath}/genomes/hg38/inhouse_DBs/hg38v3_primary/"
        }

        // Gene and transcript annotation files:

        gencode_gtf = "${syspath}/genomes/hg38/gene.annotations/gencode.v36.annotation.gtf"
        gencode_gff3 = "${syspath}/genomes/hg38/gene.annotations/gencode.v36.annotation.gff3"
     
        //Program  files:
        msisensor_list="${syspath}/genomes/hg38/program_DBs/msisensor/hg38_msisensor_scan.txt"
        
        accucopy_config="${syspath}/genomes/hg38/accucopy/accucopy.docker.nextflow.conf"
        cnvradar_anno="${syspath}/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz"
        cnvradar_anno_idx="${syspath}/genomes/hg38/program_DBs/cnvradar/All_20180418.vcf.gz.tbi"
        cnvradar_ROI="${syspath}/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed" 

        cnvradar_roisum_dir="${syspath}/genomes/hg38/program_DBs/cnvradar/inhouse_roi_summaries/"
        
        //Structural variants
        delly_exclude="${syspath}/genomes/hg38/program_DBs/delly/human.hg38.excl.tsv"
        
        smoove_exclude="${syspath}/genomes/hg38/interval.files/smoove/smoove.hg38.excluderegions.bed"
        smoove_gff="${syspath}/genomes/hg38/gene.annotations/GRCh38_latest_genomic.gff.gz"


        //Repeat Expansions:
        expansionhunter_catalog="/data/shared/genomes/hg38/program_DBs/expansionHunter/expansionHunter_hg38_stripy.variant_catalog.json"
        hipSTR_bed="${syspath}/genomes/hg38/interval.files/STRs/GRCh38.hipstr_reference.bed"

        // Somatic calling files (GATK Mutect2 pipeline):
        gatk_wgs_pon="${syspath}/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz"
        mutect_gnomad="${syspath}/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_af-only-gnomad.hg38.vcf.gz"
        gatk_contamination_ref="${syspath}/genomes/hg38/program_DBs/GATK/somatic/somatic-hg38_small_exac_common_3.hg38.vcf.gz"

        // Program indexes:
        pcgr_assembly="grch38"
        sequenza_cg50_wig="${syspath}/genomes/hg38/program_DBs/sequenza/GRCh38.primary.cg50.sequenza.wig.gz"
        spliceai_assembly="grch38"

        // Regions & variants:
        qualimap_ROI="${syspath}/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.6col.bed"
        gencode_exons_ROI="${syspath}/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.SM.bed"

        ROI="${syspath}/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        //ROI="${syspath}/genomes/hg38/interval.files/210129.hg38.gencode36.codingexons.20bp.SM.bed"

        callable_regions="${syspath}/genomes/hg38/interval.files/GATK.hg38.callable.regions.bed"
        manta_callable_regions="${syspath}/genomes/hg38/interval.files/manta/GATK.hg38.callable.regions.bed.gz"

        dbsnp="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
        KGindels="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz"
        KGindels_idx="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"

        KGmills="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        KGmills_idx="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        KG_p1_High_snps="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf.gz"

        hapmap="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz"
        omni="${syspath}/genomes/hg38/program_DBs/GATK/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
        
        AV1_ROI="${syspath}/genomes/${params.genome}/interval.files/panels/av1.hg38.ROI.bed"
        WES_ROI="${syspath}/genomes/hg38/interval.files/exome.ROIs/211130.hg38.refseq.gencode.fullexons.50bp.SM.bed"
        
        break;
}




///////////////// Input CRAM and CRAI Files Channels //////////////////////

// Create a channel for normal CRAM files and split the file name to get the base name.
Channel
.fromPath(params.normal_cram)
.map { tuple(it.baseName.tokenize('.').get(0),it) }
.into {normalCram1; normalCram2; normalCram3}

// Create a channel for normal CRAI files and split the file name to get the base name.
Channel
.fromPath(params.normal_Crai)
.map { tuple(it.baseName.tokenize('.').get(0),it) }
.into {normalCrai1; normalCrai2; normalCrai3}

// Create a channel for tumor CRAM files and split the file name to get the base name.
Channel
.fromPath(params.tumor_Cram)
.map { tuple(it.baseName.tokenize('.').get(0),it) }
.into {tumorCram1; tumorCram2; tumorCram3}

// Create a channel for tumor CRAI files and split the file name to get the base name.
Channel
.fromPath(params.tumor_Crai)
.map { tuple(it.baseName.tokenize('.').get(0),it) }
.into {tumorCrai1; tumorCrai2; tumorCrai3}



// Join CRAI files into CRAM files
// Join normal CRAM and CRAI files based on their base names.
normalCram1 
    .join(normalCrai1)
    .into {normal_cram_crai1; normal_cram_crai2; normal_cram_crai3}

// Join tumor CRAM and CRAI files based on their base names.
tumorCram1 
    .join(tumorCrai1)
    .into {tumor_cram_crai1; tumor_cram_crai2; tumor_cram_crai3}




// Script execution for Strelka workflow
// Define the Strelka workflow process with inputs and output. The script includes configuration and execution of the Strelka somatic variant caller.

process strelka2_singularity {

    input:
    // Input specification for the Strelka workflow
    // Define inputs for the Strelka workflow: metadata, normal CRAM and index, and tumor CRAM and index.
    tuple val(meta), path(normal_cram), path(normal_index) from normal_cram_crai1
    tuple val(meta), path(tumor_cram), path(tumor_index) from tumor_cram_crai1

    // path normalBam from params.normalCram
    // path tumorBam from params.tumorCram
    // path referenceFasta from params.referenceFasta
    // path genome_fasta from params.genome_fasta

    output:
    

    script:
    """
    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif /tools/strelka2/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam ${normal_cram} \
    --tumorBam ${tumor_cram} \
    --referenceFasta ${genome_fasta} \
    --exome \
    --runDir NF_Strelkatest_singularity

    singularity run -B /data/shared/programmer/simg/strelka2_2.9.10.sif python2 ${params.runDir}/runWorkflow.py -j 10 -m local

    singularity run -B ${s_bind} ${simgpath}/strelka2_2.9.10.sif python2 NF_Strelkatest_singularity/runWorkflow.py \
    -j 10 \
    -m local
    """
}
