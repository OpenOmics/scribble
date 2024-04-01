#!/bin/bash
BAM_DIR="/data/NHLBI_IDSS/dev/NIAMS-40/RNA_hg38_42/bams"
BAM_GLOB="*.star_rg_added.sorted.dmark.bam"
GENOME_FA="/data/CCBR_Pipeliner/db/PipeDB/Indices/hg38_30/ref.fa"
MOD_DIR="/data/OpenOmics/references/deepvariant/RNAseek_DV_model"
REGIONS="/data/OpenOmics/references/genome-seek/wes_bedfiles/gencode_v44_protein-coding_exons.bed"
echo "#SWARM -t 24 -g 64 --time 10-00:00:00 --gres=lscratch:500 --module singularity,deepvariant --partition norm" > dv_varcall_se.swarm
for bam in $BAM_DIR/$BAM_GLOB
do
    bname=$(basename -- "$bam")
    arrIN=(${bname//./ })
    sid=${arrIN[0]}
    outdir="/data/NHLBI_IDSS/dev/NIAMS-40/PE_var_calling/$sid"
    DV_CMD="run_deepvariant --model_type=WES --num_shards=12 --make_examples_extra_args=\"split_skip_reads=true,channels=''\" --customized_model=$MOD_DIR/model.ckpt --ref=$GENOME_FA --regions=$REGIONS --reads=$bam --output_gvcf=$outdir/$sid.output.g.vcf.gz --output_vcf=$outdir/$sid.output.vcf.gz --intermediate_results_dir $outdir/int_res/$sid"
    echo "mkdir -p $outdir/int_res/$sid; $DV_CMD" >> dv_varcall_se.swarm
done

