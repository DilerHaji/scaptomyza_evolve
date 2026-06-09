#!/bin/bash

if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_variants> <output_prefix> <num_chunks> <output_dir>"
    exit 1
fi

input_variants="$1"
output_prefix="$2"
num_chunks="$3"
output_dir="$4"

mkdir -p "$output_dir"

total_lines=$(wc -l < "$input_variants")

lines_per_chunk=$(( (total_lines + num_chunks - 1) / num_chunks ))

split -l "$lines_per_chunk" --numeric-suffixes=1 --additional-suffix=.txt "$input_variants" "${output_prefix}"

for i in $(seq 1 "$num_chunks"); do
    old_name=$(printf "${output_prefix}%02d.txt" $i)
    new_name="${output_prefix}${i}.txt"
    if [ -f "$old_name" ] && [ "$old_name" != "$new_name" ]; then
        mv "$old_name" "$new_name"
    elif [ ! -f "$new_name" ]; then
        touch "$new_name"
    fi
done

for i in $(seq $((num_chunks + 1)) 99); do
    if [ -f "${output_prefix}${i}.txt" ]; then
        rm "${output_prefix}${i}.txt"
    fi
done
