
<h1 style="color:skyblue;">Missing genes in squamata</h1>

This GitHub repository contains the data for the paper **"Gene loss through chromosomal rearrangements and segmental deletions in squamates"**

Buddhabhushan Girish Salve, Sonal H, Nagarjun Vijay

Computational Evolutionary Genomics Lab, Department of Biological Sciences, IISER Bhopal, Bhauri, Madhya Pradesh, India

*Author for Correspondence: Nagarjun Vijay, nagarjun@iiserb.ac.in

____________________________________________________________________________________________________________________________________________________

____________________________________________________________________________________________________________________________________________________
**Folder Structure**

├── CSF1R-Selection_analysis/                          # Selection analysis of CSF1R **Fig. 7**

├── GC_content_of_IL34-LAPTM5-STAP1-TNIP2/             # GC content analysis **Fig. 3**

├── IL34_assembly_verification/                        # Assembly verification for _IL34_ **Fig. 6**

├── Macrophage_related/                                # Genes and pathways involved in macrophage regulation **Fig. 8**

├── Missing_genes_in_squamates-BLASTn-Genome/          # BLASTn of missing genes against 365 squadmates genome assemblies, **Fig. 1** and **Fig. 3**

├── Missing_genes_in_squamates-BLASTn-RNA-seq_database/# BLASTn of missing genes against RNA-seq database

├── Missing_genes_in_squamates-BioMart-Paralogs_of_missing_genes/ # Paralogs of missing genes from Ensembl BioMart **Fig. 4**

├── Missing_genes_in_squamates-GC_content_seq_divergence/         # GC content and sequence divergence analysis  **Fig. 3**

├── Missing_genes_in_squamates-Gene_status_in_tuatara/            # Gene presence/absence check in tuatara

├── Missing_genes_in_squamates-NCBI_GDV-chicken_green_anole/      # Visualization in NCBI GDV for chicken and anole

├── Missing_genes_in_squamates-Pogona_vitticeps_synteny_of_13_genes/ # Synteny mapping of 13 genes in Pogona vitticeps

├── Missing_genes_in_squamates-Screening/                         # Screening of candidate gene losses

├── Missing_genes_in_squamates-pogona_vitticeps_assembly_verfication/ # Assembly assessment of Pogona vitticeps

├── Missing_genes_in_squamates-summary/                           # Summary files and results

├── Segmental_deletion-IL34-STAP1-LAPTM5/             # Evidence of segmental deletions involving _IL34_, _STAP1_, _LAPTM5_ **Fig. 5**

├── Synteny_of_IL34-LAPTM5-STAP1-TNIP2-CSF1R/         # Synteny maps of the focal gene cluster **Fig. 5** and **Fig. 8**

├── TNIP2-Interchromosomal_rearrangement/             # Rearrangement events of _TNIP2_ **Fig. 8**

├── Cross-species_RNA-seq_expression.7z               # RNA-seq expression analysis of _IL34_, _CSF1R_, _CSF1_, etc.

├── Missing_genes_in_squamates-HybPiper.7z            # HybPiper output for missing gene analysis

├── Missing_genes_in_squamates-Pogona_vitticeps_RNA-seq_mapping.7z # Cross-species RNA-seq mapping 

├── Missing_genes_in_squamates-Segmental_deleltion.7z # Archive with segmental deletion evidence

├── Missing_genes_in_squamates-based_on_chain_files_of_chicken.7z  # Chain file–based alignments (chicken) **Fig. 1** and **Fig. 2**

├── Missing_genes_in_squamates-based_on_chain_files_of_human.7z    # Chain file–based alignments (human)

├── Tail_regeneration.7z                            # Tail regeneration-related gene expression or analysis

____________________________________________________________________________________________________________________________________________________
**Prerequisites:**
- NCBI datasets
- GFA v2
- MafFilter (1.3.1)
- minimap2 (2.29-r1283)
- FLAIR (2.2.0)
- seqkit
- EMBOSS-needle
- PAML (4.9f)
- BLASTN (2.13.0)
- phastBias(1.6)
- mapnh(1.3.0)
- HYPHY 2.5.48(MP) and 2.5.62(MP)
- Samtools (1.3)
- BEDTools (2.27.1)
- seqtk (1.2-r94)
- Guidance (v2.01)
- PRANK (v.170427)
- IGV (2.8.13)
- igv-reports (1.12.0)
- STAR (2.7.0d)
- R (4.1.0)
- HybPiper (2.3.1)
- IQ-TREE (2.3.6)
- klumpy (1.0.11)
- kallisto (0.51.0)

**R packages:**
- ape (5.5)
- pheatmap
- phylolm
- GENESPACE (1.3.1)
- tidyheatmap (1.11.6)
- phytools (0.7-90)
- ggplot2 (3.3.5)
- ggrepel (0.9.1)
- cowplot (1.1.1)
- dplyr (1.0.7)
- ggplotify (0.1.0)
- grid (4.1.1)
- gridExtra (2.3)
- reshape2 (1.4.4)
