#!/bin/bash

# Script to pair FASTQ files and validate R1/R2 matches
# Usage: ./pair_fastq.sh file1_R1.fastq file1_R2.fastq file2_R1.fastq.gz file2_R2.fastq.gz ...

set -euo pipefail

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

function process_files() {

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
        
        if [[ "$basename" =~ ^(.+)_R1\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R1"
        elif [[ "$basename" =~ ^(.+)_R2\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R2"
        elif [[ "$basename" =~ ^(.+)\.R1\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R1"
        elif [[ "$basename" =~ ^(.+)\.R2\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R2"
        elif [[ "$basename" =~ ^(.+)_1\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R1"
        elif [[ "$basename" =~ ^(.+)_2\.(fastq|fq)(\.gz)?$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            read_type="R2"
        else
            echo "Error: Cannot determine R1/R2 from filename '$basename'"
            echo "Expected patterns: *_R1.fastq[.gz], *_R2.fastq[.gz], *.R1.fastq[.gz], *.R2.fastq[.gz], *_1.fastq[.gz], *_2.fastq[.gz]"
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


    # do swarming
}

main "$@"