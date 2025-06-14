#!/bin/bash

input_file="co-ordinates_1.txt"

mkdir -p plots

while IFS= read -r line; do
    gene=$(echo "$line" | cut -d':' -f1)
    chrom=$(echo "$line" | cut -d':' -f2)
    bounds=$(echo "$line" | cut -d':' -f3)
    left=$(echo "$bounds" | cut -d'-' -f1)
    right=$(echo "$bounds" | cut -d'-' -f2)

    output_name="${gene}_${chrom}_${left}_${right}.svg"

    # Run the plot (assume default output is alignment_plot.svg)
    klumpy alignment_plot \
        --alignment_map Nanopore_merged.sorted.bam \
        --reference "$chrom" \
        --candidates Nanopore_merged.sorted_Candidate_Regions.tsv \
        --window_size 10000 \
        --window_step 5000 \
        --color red \
        --vertical_line_gaps \
        --vertical_line_klumps \
        --format svg \
        --leftbound "$left" \
        --rightbound "$right" \
        --annotation GCF_047335585.1_Pvi2.1_genomic.gtf \
        --gap_file GCF_047335585.1_Pvi2.1_genomic_gaps.tsv \
        --height 1500 \
        --width 4000 \
        --min_len 10000 \
        --vertical_line_exons \
        --number

    # Rename output file
    mv alignment_plot.svg "plots/$output_name"

    echo "Saved: plots/$output_name"

done < "$input_file"
