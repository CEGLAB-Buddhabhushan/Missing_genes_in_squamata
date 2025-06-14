
library(GENESPACE)

genomeRepo <- "/home/buddha/work/BMCBIO-EoI/TNIP2_EBR/rawGenomes/"
wd <- "/home/buddha/work/BMCBIO-EoI/TNIP2_EBR/"
path2mcscanx <- "~/software/GENESPACE/MCScanX/"

urls <- c(
Gallus_gallus = "016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_",
Alligator_mississippiensis = "030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_",
Chelonia_mydas = "015/237/465/GCF_015237465.2_rCheMyd1.pri.v2/GCF_015237465.2_rCheMyd1.pri.v2_",
Pogona_vitticeps = "900/067/755/GCF_900067755.1_pvi1.1/GCF_900067755.1_pvi1.1_",
Anolis_carolinensis = "035/594/765/GCF_035594765.1_rAnoCar3.1.pri/GCF_035594765.1_rAnoCar3.1.pri_",
Crotalus_tigris = "016/545/835/GCF_016545835.1_ASM1654583v1/GCF_016545835.1_ASM1654583v1_",
Zootoca_vivipara = "963/506/605/GCF_963506605.1_rZooViv1.1/GCF_963506605.1_rZooViv1.1_",
Podarcis_muralis = "004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_")

Gallus_gallus https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/699/485/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_translated_cds.faa.gz
Alligator_mississippiensis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/030/867/095/GCF_030867095.1_rAllMis1/GCF_030867095.1_rAllMis1_translated_cds.faa.gz
Chelonia_mydas https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/237/465/GCF_015237465.2_rCheMyd1.pri.v2/GCF_015237465.2_rCheMyd1.pri.v2_translated_cds.faa.gz
Pogona_vitticeps https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/067/755/GCF_900067755.1_pvi1.1/GCF_900067755.1_pvi1.1_translated_cds.faa.gz
Anolis_carolinensis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/594/765/GCF_035594765.1_rAnoCar3.1.pri/GCF_035594765.1_rAnoCar3.1.pri_translated_cds.faa.gz
Crotalus_tigris https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/545/835/GCF_016545835.1_ASM1654583v1/GCF_016545835.1_ASM1654583v1_translated_cds.faa.gz
Zootoca_vivipara https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/963/506/605/GCF_963506605.1_rZooViv1.1/GCF_963506605.1_rZooViv1.1_translated_cds.faa.gz
Podarcis_muralis https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/329/235/GCF_004329235.1_PodMur_1.0/GCF_004329235.1_PodMur_1.0_translated_cds.faa.gz

genomes2run <- names(urls)
urls <- file.path("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/", urls)
translatedCDS <- sprintf("%stranslated_cds.faa.gz", urls)
geneGff <- sprintf("%sgenomic.gff.gz", urls)

names(translatedCDS) <- genomes2run
names(geneGff) <- genomes2run
writeDirs <- file.path(genomeRepo, genomes2run)
names(writeDirs) <- genomes2run
for(i in genomes2run){
  print(i)
  if(!dir.exists(writeDirs[i]))
    dir.create(writeDirs[i])
  download.file(
    url = geneGff[i], 
    destfile = file.path(writeDirs[i], basename(geneGff[i])))
  download.file(
    url = translatedCDS[i], 
    destfile = file.path(writeDirs[i], basename(translatedCDS[i])))
}

genomes2run <- c("Gallus_gallus", "Alligator_mississippiensis", "Chelonia_mydas", "Pogona_vitticeps", "Anolis_carolinensis", "Crotalus_tigris", "Zootoca_vivipara", "Podarcis_muralis")
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = genomes2run,
  genomeIDs = genomes2run,
  presets = "ncbi",
  genespaceWd = wd)

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx)
  
out <- run_genespace(gpar, overwrite = T)
ripd <- plot_riparian(
  gsParam = out,
  refGenome = "Gallus_gallus", 
  useRegions = FALSE)
  

  
  ripDat <- plot_riparian(
  gsParam = out, 
  genomeIDs = c("Sarcophilus_harrisii", "Thylacinus_cynocephalus"), 
  refGenome = "Sarcophilus_harrisii", 
  inversionColor = "green")
  
  ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "black", "darkblue", "purple", "darkred", "salmon"))
ripDat <- plot_riparian(
  gsParam = out, 
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  genomeIDs = c("Sarcophilus_harrisii", "Thylacinus_cynocephalus"), 
  refGenome = "Sarcophilus_harrisii")
  
  
  
roi <- data.frame(
  genome = c("Gallus_gallus"),
  chr = c("4"),
  start = c(81631120),
  end = c(82274877))  

qreturn <- query_hits(gsParam = out, bed = roi, synOnly = TRUE)
test <- query_pangenes(
  gsParam = out, bed = roi)
roibed <- roi[,c("genome", "chr")]
roibed$color <- c("pink")
ripd <- plot_riparian(
  gsParam = out, 
  useRegions = FALSE, 
  highlightBed = roibed)
  
ripd <- plot_riparian(
  gsParam = out,  
  useRegions = FALSE,
  highlightBed = roibed, 
  backgroundColor = NULL)
  
  
ripd <- plot_riparian(
  gsParam = out,
  refGenome = "Gallus_gallus", 
  useRegions = FALSE) 


  ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))
customPal <- colorRampPalette(
  c("darkorange", "skyblue", "darkblue", "purple", "darkred", "salmon"))
ripDat <- plot_riparian(
  gsParam = out, 
  palette = customPal,
  braidAlpha = .75,
  chrFill = "lightgrey",
  addThemes = ggthemes,
  refGenome = "Gallus_gallus")
  
jpg(file="Malassezia.jpg",width=16, height=9, units="in", res=600)
print (ripDat)
dev.off()

