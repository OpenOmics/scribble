#! /bin/bash
#SBATCH --job-name=Bam2Cram
#SBATCH --mail-user=markowitzte@nih.gov
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=1-00:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=18

#### 
## Variables
# FOLDER: the full path to the folder you want converted
# DIRECTION: options are "Bam2Cram" or "Cram2Bam"
# GENOME: the path to the reference fasta file
# example usage: sbatch --export FOLDER="/data/NCBR/projects/X",DIRECTION="Bam2Cram",GENOME="Y.fa" Bam2Cram.sbatch
#####

module load samtools

cd $FOLDER

if [[ $DIRECTION == "Bam2Cram" ]]; then
   RE="(.+)\.bam"
   
   for BAM in `ls *bam`; do
       if [[ $BAM =~ $RE ]]; then
       	     ROOT=${BASH_REMATCH[1]}
		fi

			CRAM="$ROOT.cram"
				samtools view -@ 16 -C -o $CRAM -T $GENOME $BAM
                                rm $BAM
				done
else
	RE="(.+)\.cram"
	
	for CRAM in `ls *cram`; do
	    if [[ $CRAM =~ $RE ]]; then
   	       	  ROOT=${BASH_REMATCH[1]}
			fi

				BAM="$ROOT.bam"
					samtools view -@ 16 -o $BAM -T $GENOME $CRAM
                                        rm $CRAM
					done
fi

