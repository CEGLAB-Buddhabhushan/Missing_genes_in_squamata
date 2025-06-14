#!/bin/bash

# Input files
GAP_FILE="GCF_047335585.1_Pvi2.1_genomic_gaps.tsv"
GENE_FILE="co-ordinates.txt"
OUTPUT="gene_with_gap_counts.tsv"

# Output header
echo -e "GeneName\tChromosome\tCoordinates\tGapCount" > "$OUTPUT"

# Read each gene line
while IFS= read -r line; do
    gene=$(echo "$line" | cut -d':' -f1)
    chrom=$(echo "$line" | cut -d':' -f2)
    coords=$(echo "$line" | cut -d':' -f3)
    start=$(echo "$coords" | cut -d'-' -f1)
    end=$(echo "$coords" | cut -d'-' -f2)

    # Count overlaps with gaps in the same chromosome
    count=$(awk -v chr="$chrom" -v gs="$start" -v ge="$end" -F'\t' '
        $1 == chr && $3 >= gs && $2 <= ge { count++ }
        END { print count+0 }
    ' "$GAP_FILE")

    # Append result to output
    echo -e "${gene}\t${chrom}\t${start}-${end}\t${count}" >> "$OUTPUT"

done < "$GENE_FILE"
