#!/bin/bash
: '
    Google Deep Variant swarm creator:
        Run deep variant on a collection of bam files in the common parallel fashion on NIH biowulf.

        This is setup to run DeepVariant in the cpu version and not the gpu version because of the bottleneck that is
        getting enough GPUs on biowulf.

        This is also setup to run on RNAseq data originally, which is not the most common use case for deep variant. You
        can revert back to one of the basic internal models of deepvariant by removing `--customized_model` argument, and
        editting `model_type` to your liking.
'
### SET THESE CONSTANTS BEFORE RUNNING ####
BAM_DIR="/dir/to/bam/files"
BAM_GLOB="*.bam" # glob subsetting bam directory
OUTPUT_DIRECTORY="/deepvariant/output/dir"
GENOME_FASTA="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_30/ref.fa"  # hg38 default genome, change if want different
MODEL_DIR="/data/OpenOmics/references/deepvariant/RNAseek_DV_model" # rna seq deep variant model, change if want different
JOB_MEMORY="128" # gigabytes of memory per swarm
JOB_CPU="64" # number of threads per swarm
JOB_WALLTIME="1000" # swarm job walltime minutes

### DO NOT CHANGE THESE UNLESS YOU NEED TO ###
SNG_CMD="singularity exec -C"
SNG_BINDS="-B /data/OpenOmics,/data/NIAMS_IDSS,/data/CCBR_Pipeliner,/lscratch/\$SLURM_JOB_ID/:/tmp,$OUTPUT_DIRECTORY:/output"
DV_IMG="/data/OpenOmics/SIFs/deepvariant_1.6.0.sif"
REGIONS="/data/NIAMS_IDSS/dev/NIAMS-40/gencode_annotations_cds.bed"
echo "#SWARM -t $JOB_CPU -g $JOB_MEMORY --time $JOB_WALLTIME --gres=lscratch:500 --module singularity --partition norm" > dv_varcall.swarm
for bam in $BAM_DIR/$BAM_GLOB
do
    bname=$(basename -- "$bam")
    DV_CMD="run_deepvariant --model_type=WES --customized_model=$MODEL_DIR/model.ckpt --ref=$GENOME_FASTA --regions=$REGIONS --reads=$bam --output_vcf=/output/vcfs/$bname.output.vcf.gz --intermediate_results_dir /output/int_res/$bname"
    FULL_SNG_CMD="$SNG_CMD $SNG_BINDS $DV_IMG $DV_CMD"
    echo "$FULL_SNG_CMD" >> dv_varcall.swarm
done

swarm dv_varcall.swarm
