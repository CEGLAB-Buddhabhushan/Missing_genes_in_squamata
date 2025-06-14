
# Missing genes in squamata
This GitHub repository contains the data for the paper **"Gene loss through chromosomal rearrangements and segmental deletions in squamates"**

Buddhabhushan Girish Salve, Sonal H, Nagarjun Vijay

Computational Evolutionary Genomics Lab, Department of Biological Sciences, IISER Bhopal, Bhauri, Madhya Pradesh, India

*Author for Correspondence: Nagarjun Vijay, nagarjun@iiserb.ac.in

____________________________________________________________________________________________________________________________________________________

____________________________________________________________________________________________________________________________________________________
**Folder Structure**
- **Chromosome_wise_Chains_and_TOGA_results:** Contains the output files of make_lastz_chains and TOGA.
- **Figure1:** Contains the gene loss confirmation genes losses in thylacine and other marsupial species, IGV-reports, and gene loss history. **(Figure 1, figure S2, figure S5, figure S7, figure S8 and figure S9)**
- **Figure2:** Contains the figure 2 figures and the inputs required to get figures like phylogenetic trees from TimeTree and iTOLs. **(Figure 2)**
- **Figure3:**
- (a) Contains the input and output files of GC-content variation in thylacine. **(Figure 3a)**
- (b) Selection pressure analysis of "clearly lost" genes (found using TOGA), tested using RELAX. **(Figure 3b)**
- (c) Expression of lost genes in Devil and Dunnart transcriptome. **(Figure 3c)**
- **Gene_Tree_SAMD9-9L:** Contains the alignment and tree files used to obtain gene trees. **(figure S3)**
- **Hybpiper_patchwork:** Contains the input files for Hybpiper and patchwork. The subset of fastq reads gives hits (BLASTn) to lost genes (as s query sequences).
- **Hybpiper_SAMD9-9L:**  In-frame stop codon in _SAMD9L_ reassessed using short-read data of thylacine, using _SAMD9_ and _SAMD9L_ codding sequence of other 20 marsupial species.
- **IGV_reports-BWA+BLASTn:**  Contains the IGV-reports for clearly lost genes in chromosome-wise folders **(Figure 1)**
- **Loss_confirmation:** Contains the output files of short-read BLASTn, bam-read count, and gene-loss events found by TOGA to generate a summary table and confirm gene loss. **(Figure 1, table S1)**
- **Main_Figures:** contains the main figures in PNG and JPG formats.
- **Myrmecobius_fasciatus_RNA_seq:** Contains an IGV screenshot for the _VWA7_ gene, subseted bam files, and BED files. **(figure S10)**
- **PAML_omega_SAMD9-9L:** Evolutionary rate dN/dS ratio (omega) across the _SAMD9_ and _SAMD9L_ genes. Contains the input and R script. **(figure S4)**
- **Selection_analysis:**  Selection pressure analysis of lost gene across 21 marsupial species, tested using RELAX, codeML, MEME, FEL, aBSREL, BUSTED, and gBGC.  **(Figure 3 and table S5)**
- **Thylacine_transcriptome-miRNA:** Contains BLASTn search output for lost genes and control genes in the miRNA database of thylacine.
- **TOGA_output_sorted:** Contains the sorted output files of TOGA, such as inact_mut_data.txt and loss_summ_data.tsv, codon.fasta, nucleotide.fasta.tgz, prot.fasta, query_annotation.bed, and query_gene_spans.bed. Also, a mutation plot was generated using TOGA. **(Figure 1)**
- **Synteny:** Contains the NCBI synteny images for _SAMD9L_, _HSD17B13_, _CUZD1_ and _VWA7_. **(figure S5 and tables S2)**
- **SAMD9_gene_loss_in_carnivore:** Contain the script for gene loss events identification by TOGA. Assembly verification output files of klumpy and igv reports of events by short read for dingo and tiger. **(figure S6)**
____________________________________________________________________________________________________________________________________________________
**Prerequisites:**
- TOGA (1.1.7)
- make_lastz_chains (https://github.com/hillerlab/make_lastz_chains.git)
- PAML (4.9f)
- BLASTN (2.13.0)
- phastBias(1.6)
- mapnh(1.3.0)
- HYPHY 2.5.48(MP) and 2.5.62(MP)
- BWA(0.7.17-r1188)
- Samtools (1.3)
- bedtools (v2.27.1)
- ea-utils (https://github.com/ExpressionAnalysis/ea-utils.git)
- seqtk (1.2-r94)
- Guidance (v2.01)
- PRANK (v.170427)
- IGV (2.8.13)
- igv-reports (1.12.0)
- STAR (2.7.0d)
- MegaX
- R (4.1.0)
- bam-readcount (1.0.1)
- HybPiper (2.3.1)
- Patchwork (0.5.2)
- IQ-TREE (2.3.6)
- klumpy (1.0.11)
- kallisto (0.51.0)
- Bowtie2 (2.5.4)

**R packages:**
- ape (5.5)
- phytools (0.7-90)
- ggplot2 (3.3.5)
- ggrepel (0.9.1)
- cowplot (1.1.1)
- dplyr (1.0.7)
- ggplotify (0.1.0)
- grid (4.1.1)
- gridExtra (2.3)
- reshape2 (1.4.4)
