#!/bin/bash
output_file=$1
shift
input_files=("$@")

header=$(head -n 1 "${input_files[0]}")
echo "$header" > "$output_file"
awk -F, 'NR>1 {print $1}' "${input_files[0]}" | sort -u > tmp_intersect.txt

for f in "${input_files[@]:1}"; do
    awk -F, 'NR>1 {print $1}' "$f" | sort -u > tmp_next.txt
    comm -12 tmp_intersect.txt tmp_next.txt > tmp_new.txt
    mv tmp_new.txt tmp_intersect.txt
done

grep -Fwf tmp_intersect.txt <(cat "${input_files[@]}" | awk 'NR==1 || FNR>1') | sort -u >> "$output_file"

rm tmp_intersect.txt tmp_next.txt 2>/dev/null
