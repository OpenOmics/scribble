#!/bin/bash

# Script to pair FASTQ files and validate R1/R2 matches
# Usage: ./<SCRIPT>.sh file1_R1.fastq file1_R2.fastq file2_R1.fastq.gz file2_R2.fastq.gz ...

set -euo pipefail
shopt -s extglob

function err() { cat <<< "$@" 1>&2; }
function fatal() { cat <<< "$@" 1>&2; exit 1; }
function usage() { err "Usage: $0 <gencode.gtf> <gencode_exons.bed>"; }

function require(){
  # Requires an executable is in $PATH, 
  # as a last resort it will attempt to load
  # the executable or dependency as a module
  # @INPUT $@ = List of executables to check

  for exe in "${@}"; do
    # Check if executable is in $PATH
    command -V ${exe} &> /dev/null && continue;
    # Try to load exe as lua module
    module load ${exe} &> /dev/null || \
      fatal "Failed to find or load '${exe}', not installed on target system."
  done
}

function main() {
    # Check if any arguments provided
    if [ $# -eq 0 ]; then
        echo "Error: No FASTQ files provided"
        echo "Usage: $0 <fastq_files...>"
        echo "Example: $0 sample1_R1.fastq sample1_R2.fastq sample2_R1.fastq.gz sample2_R2.fastq.gz"
        exit 1
    fi

    # Associative arrays to store R1 and R2 files
    declare -A r1_files
    declare -A r2_files
    declare -A all_samples

    echo "Processing ${#} FASTQ files..."

    # Process each file argument
    for file in "$@"; do
        # Check if file exists
        if [[ ! -f "$file" ]]; then
            echo "Error: File '$file' does not exist"
            exit 1
        fi
        
        # Check if file has valid FASTQ extension
        if [[ ! "$file" =~ \.(fastq|fq)(\.gz)?$ ]]; then
            echo "Error: File '$file' does not have a valid FASTQ extension (.fastq, .fq, .fastq.gz, .fq.gz)"
            exit 1
        fi
        
        # Extract basename and determine if R1 or R2
        basename=$(basename "$file")
        
        # Try different common R1/R2 naming patterns
        sample_name=""
        read_type=""
        read_regex="(.+)[_\.]R(1|2)[_\.].+\.fastq\.?(gz)?"

        [[ $basename =~ $read_regex ]]
        sample_name="${BASH_REMATCH[1]}"
        rtype="${BASH_REMATCH[2]}"
        read_type="R${rtype}"
        gz="${BASH_REMATCH[3]}"
        echo $read_type
        echo $gz
        echo $sample_name

        if [[ $read_type =~ "R[1-2]" ]]; then
            echo "Error: Cannot determine R1/R2 from filename '$basename'"
            echo "Expected patterns: *_R1_*.fastq[.gz], *_R1_*.fastq[.gz], *.R1.fastq[.gz], *.R2.fastq[.gz], *_1.fastq[.gz], *_2.fastq[.gz]"
            exit 1
        fi 
        
        echo "  Found: $sample_name -> $read_type ($file)"
        
        # Store the file path
        if [[ "$read_type" == "R1" ]]; then
            if [[ -n "${r1_files[$sample_name]:-}" ]]; then
                echo "Error: Duplicate R1 file found for sample '$sample_name'"
                echo "  First:  ${r1_files[$sample_name]}"
                echo "  Second: $file"
                exit 1
            fi
            r1_files["$sample_name"]="$file"
        else
            if [[ -n "${r2_files[$sample_name]:-}" ]]; then
                echo "Error: Duplicate R2 file found for sample '$sample_name'"
                echo "  First:  ${r2_files[$sample_name]}"
                echo "  Second: $file"
                exit 1
            fi
            r2_files["$sample_name"]="$file"
        fi
        
        # Track all sample names
        all_samples["$sample_name"]=1
    done

    echo ""
    echo "Validating paired samples..."

    # Check that every sample has both R1 and R2
    missing_pairs=0
    orphaned_files=()

    for sample in "${!all_samples[@]}"; do
        has_r1=false
        has_r2=false
        
        if [[ -n "${r1_files[$sample]:-}" ]]; then
            has_r1=true
        fi
        
        if [[ -n "${r2_files[$sample]:-}" ]]; then
            has_r2=true
        fi
        
        if [[ "$has_r1" == true && "$has_r2" == true ]]; then
            echo "  ✓ $sample: PAIRED"
            echo "    R1: ${r1_files[$sample]}"
            echo "    R2: ${r2_files[$sample]}"
        elif [[ "$has_r1" == true && "$has_r2" == false ]]; then
            echo "  ✗ $sample: MISSING R2"
            echo "    R1: ${r1_files[$sample]}"
            echo "    R2: (missing)"
            orphaned_files+=("${r1_files[$sample]}")
            ((missing_pairs++))
        elif [[ "$has_r1" == false && "$has_r2" == true ]]; then
            echo "  ✗ $sample: MISSING R1"
            echo "    R1: (missing)"
            echo "    R2: ${r2_files[$sample]}"
            orphaned_files+=("${r2_files[$sample]}")
            ((missing_pairs++))
        fi
    done

    echo ""

    if ! [[ $missing_pairs -eq 0 ]]; then
        echo "ERROR: Found $missing_pairs sample(s) with missing pairs!"
        echo ""
        echo "Orphaned files:"
        for file in "${orphaned_files[@]}"; do
            echo "  $file"
        done
        echo ""
        echo "Please ensure all samples have both R1 and R2 files."
        exit 1
    fi

    common_parent_dir=$(./common_parent.sh "$@")
    if [[ $common_parent_dir == "/" ]]; then
        echo "Common parent path between files too short! Please relocate files to same directory!"
    fi

    # swarm setup
    swarm_file="centrifudge.swarm"
    [ -f $swarm_file ] && rm $swarm_file

    swarm_head="#SWARM -t 12 -g 16 --time 4:00:00 --partition norm"
    echo $swarm_head >> $swarm_file

    swarm_head="#SWARM --sbatch '--mail-type=FAIL --chdir=$PWD'"
    echo $swarm_head >> $swarm_file

    swarm_head="#SWARM --logdir $PWD/swarm_logs"
    echo $swarm_head >> $swarm_file

    swarm_head="#SWARM --module singularity"
    echo $swarm_head >> $swarm_file

    index_path="/data/OpenOmics/references/centrifuger/refseq_custom"
    img="/data/OpenOmics/SIFs/centrifudger_sylph_0.0.4.sif"
    binds="$common_parent_dir:$common_parent_dir,$PWD:/work2,/data/OpenOmics/references/centrifuger/refseq_custom:/data/OpenOmics/references/centrifuger/refseq_custom"

    for sample in "${!all_samples[@]}"; do
        # echo $sample
        this_r1="${r1_files[$sample]}"
        this_r2="${r1_files[$sample]}"
        singularity_base="singularity run -c --cwd /work2 -B $binds $img "
        singularity_cmd="centrifuger -t 10 -1 $this_r1 -2 $this_r2 -x /data/OpenOmics/references/centrifuger/refseq_custom/refseq_abv > ${sample}_classification.tsv; "
        singularity_cmd+="centrifuger-quant -x /data/OpenOmics/references/centrifuger/refseq_custom/refseq_abv -c ${sample}_classification.tsv > ${sample}_centrifudge_report.tsv; "
        singularity_cmd+="centrifuger-quant -x /data/OpenOmics/references/centrifuger/refseq_custom/refseq_abv --output-format 1 -c ${sample}_classification.tsv > ${sample}_metaphlan_report.tsv"
        singularity_full="$singularity_base \"$singularity_cmd\""
        echo $singularity_full >> $swarm_file
    done
}

main "$@"