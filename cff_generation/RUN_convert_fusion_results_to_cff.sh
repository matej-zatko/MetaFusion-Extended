#!/bin/bash

# Usage instructions
usage() {
  echo "Usage: $0 -c <caller_file_dir> -s <sampleinfo> -d <dataset> -o <outdir> -t \"<tool1 tool2 ...>\""
  exit 1
}

# Parse command-line arguments
while getopts "c:s:d:o:t:" opt; do
  case "$opt" in
    c) caller_file_dir="$OPTARG" ;;
    s) sampleinfo="$OPTARG" ;;
    d) dataset="$OPTARG" ;;
    o) outdir="$OPTARG" ;;
    t) tools=($OPTARG) ;;  # Convert space-separated tools into an array
    *) usage ;;
  esac
done

# Ensure all required arguments are provided
if [[ -z "$caller_file_dir" || -z "$sampleinfo" || -z "$dataset" || -z "$outdir" || ${#tools[@]} -eq 0 ]]; then
  echo "Error: Missing required arguments."
  usage
fi

# make output dir if doesnt exist
mkdir -p "$outdir"

# Iterate over each specified fusion tool
for tool in "${tools[@]}"; do 
  raw_file_dir="$caller_file_dir/$dataset/$tool"
  
  # Ensure the directory exists
  if [[ ! -d "$raw_file_dir" ]]; then
    echo "Warning: Directory $raw_file_dir does not exist. Skipping $tool."
    continue
  fi

  result_files=$(ls "$raw_file_dir"/*/* 2>/dev/null)

  echo "Generating CFF for $tool"

  for result_file in ${result_files[@]}; do
    sample=$(basename "$(dirname "$result_file")")
    echo "$sample, $tool, $result_file, $outdir"
    
    # Run the Python script
    python /MetaFusion/cff_generation/convert_fusion_results_to_cff_extended.py \
      --sample "$sample" \
      --sample_info_file "$sampleinfo" \
      --tool "$tool" \
      --fusion_result_file "$result_file" \
      --outdir "$outdir"
  done
done

# Merge all .cff files into one "merged.cff"
merged_file="$outdir/merged.cff"
files=$(ls "$outdir"/*cff 2>/dev/null | grep -v "merged.cff")

if [[ -n "$files" ]]; then
  cat $files > "$merged_file"
  echo "Merged CFF file created: $merged_file"
else
  echo "No CFF files found to merge."
fi
