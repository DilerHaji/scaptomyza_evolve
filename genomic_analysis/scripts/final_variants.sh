#!/bin/bash

output_file=$1
shift
input_files=$@

sort -u -m <(
    for file in $input_files; do
        awk -F, 'NR>1 {print $1 "\t" $2}' "$file"
    done
) > "$output_file"