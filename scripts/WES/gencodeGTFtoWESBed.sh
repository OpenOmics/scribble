#!/usr/bin/env bash

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


function main() {
  # Main entry of script

  # Check if asked for usage
  { [ "${1:-}" == "-h" ] || [ "${1:-}" == "--help" ]; } && \
    usage && exit 0

  # Check if user provided a
  # input gtf file
  { [ -z "${1:-}" ] || [ -z "${2:-}" ]; } && usage && \
   err && \
   fatal "Error: Please provide a GENCODE GTF and output BED filename as input!"

  # Check for required tools
  require bedtools gawk tee

  # Create output directory
  # as needed
  outdir=$(dirname "$2")
  if [ ! -d "$outdir" ]; then
    mkdir -p "$outdir" || \
    fatal "Failed to create output directory: ${outdir}"
  fi

  # BED files are 0-based while
  # GTF files are 1-based coordinates,
  # NR gets the line number of the
  # entry, useful for cross-lookups,
  # if needed later
  awk -F '\t' -v OFS='\t' '$3 == "exon"  && $4 < $5 {print $1,int($4)-1,$5,NR,$6,$7,$9}' "$1" | \
    # Find protein coding entries
    grep 'transcript_type "protein_coding"' | \
    # Remove deplicates, keep first
    awk '!dup[$1,$2,$3,$6]++' | \
    # Parse exon_id, gene_id, and transcript_id
    # with gawk regex captured groups
    awk -F '\t' -v OFS='\t' '{match($0, /exon_id "(\w+.\w+)";/, e); match($0, /gene_id "(\w+.\w+)";/, g); match($0, /transcript_id "(\w+.\w+)";/, t); print $1,$2,$3,$4"_"e[1],$5,$6,g[1],t[1],$4}' | \
    # Sort the resulting bed-like file
    bedtools sort | \
    # Merge overlapping regions on the 
    # same strand, collapse and output 
    # metadata fields:
    #  4. set(GTFEntryLine#_ExonIDs)
    #  5. set(Score)
    #  6. set(Strand)
    #  7. set(GeneIDs)
    #  8. set(TranscriptIDs)
    #  9. min(set(GTFEntryLine))
    #  10. max(set(GTFEntryLine))
    #  11. count(set(GTFEntryLine))
    bedtools merge -s  -c 4,5,6,7,8,9,9,9 -o distinct,distinct,distinct,distinct,distinct,min,max,count -delim '_' | \
    tee "${2%.bed}.gtf2bed_build.log" | \
    # Create name (4th) column in BED
    # file, formatted as follows:
    # Count_MergedExons_set(GeneID)_Line#FirstMergedExonEntry_Line#LastMergedExonEntry
    awk -F '\t' -v OFS='\t' '{print $1,$2,$3,$11"_"$7"_"$9"_"$10,$5,$6}' > "$2"
}


# Call main method
main "$@"
