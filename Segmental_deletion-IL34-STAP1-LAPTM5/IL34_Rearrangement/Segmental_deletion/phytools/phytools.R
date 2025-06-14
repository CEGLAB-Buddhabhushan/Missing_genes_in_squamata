setwd("/home/ceglab358/BUDDHA/BMCBIO-EoI/IL34_Rearrangement/Segmental_deletion/phytools/Final_phytools/")
library(phytools)
a <- read.table("gene_diatance.removed_outliers.csv", header=F, row.names=1, sep=',')
a
# Calculate row-wise sum of V3, V4, and V5
row_sums <- rowSums(a[, c("V3", "V4", "V5")])

# Compute percentages
a$V3 <- (a$V3 / row_sums) * 100
a$V4 <- (a$V4 / row_sums) * 100
a$V5 <- (a$V5 / row_sums) * 100
# Print result
print(a)
#write.table(a, "gene_distance_percentages.csv", sep=",", row.names=TRUE, col.names=TRUE, quote=FALSE)

t=read.tree("Species.no_internodes.nwk")
a$V3 = a$V3/100
a$V4 = a$V4/100
a$V5 = a$V5/100
column_names=c("COG4.SF3B3", "SF3B3.MTSS2",	"MTSS2.VAC14")
head(a)
i=1
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c=as.data.frame(fit)
colnames(c)=column_names[i]
i=2
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c[ , ncol(c) + 1] <- fit
i=3
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c[ , ncol(c) + 1] <- fit
colnames(c)=column_names
b=a[,-1]
colnames(b)=column_names
d=rbind(b,c)
ma=as.matrix(d[,])
geo = factor(a$V2)
mycol <- c("red","blue")[geo]
png("1Gene_diatance_phylogeny_plot.png",height=30,width=50,units = "in",res=600)
par(mfrow= c(1, 2),mar=c(3.5,2.5,4,1))
plot(t, tip.color = mycol,edge.width=3,adj=1,cex=1.5)
#num_tips <- length(t$tip.label)
#num_internal_nodes <- Nnode(t)
#total_nodes <- num_tips + num_internal_nodes
#total_nodes
nodelabels(node=1:209,pie=ma[,],piecol=c("red","yellow","blue"),cex=0.2,adj=c(-1.75,0.5))
library(phylolm)
a <- read.table("gene_diatance.removed_outliers.csv", header=F, sep=',')
# Calculate row-wise sum of V3, V4, and V5
row_sums <- rowSums(a[, c("V3", "V4", "V5")])

# Compute percentages
a$V3 <- (a$V3 / row_sums) * 100
a$V4 <- (a$V4 / row_sums) * 100
a$V5 <- (a$V5 / row_sums) * 100
head(a)
colnames(a)=c("species","gene_status","COG4.SF3B3", "SF3B3.MTSS2",	"MTSS2.VAC14")
rownames(a)=a[,1]
###############MPLE###############
#fitwmple=phyloglm(gene_status~SF3B3.MTSS2,a,t, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc1=coef(fitwmple)
summary(fitwmple)
mpleinter=round(cc1[1],2)
mplew=round(cc1[2],2)
#####################IG10#########
#fitwig10=phyloglm(gene_status~SF3B3.MTSS2,a,t, method = c("logistic_IG10"), btol = 20, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc2=coef(fitwig10)
summary(fitwig10)
ig10inter=round(cc2[1],2)
ig10w=round(cc2[2],2)
####################################
t(table(a$gene_status,a$SF3B3.MTSS2))->M
plot(rep(c(0,13.6463361123310,20.030165912518910,49.260756310069201,52.478839177750910,52.538291108756410,59.810161920714710,66.042247907532910,67.829642496846101,70.26601501476501,71.168834245374601,71.690968550984801,74.875782112624201,75.232714463182601,75.666485310119710,77.421236872812110,77.557884203796901,78.331589176259501,78.453966437405801,79.293424926398401,79.357136352089101,79.39782268107501,79.432073861432201,79.623768653490501,79.768649669499501,82.315294117647101,82.357891515568801,82.932713018006401,83.308144416456810,83.678226876591501,83.821815154038301,85.219079312257301,85.756150753083901,85.838860230901501,86.034015025041701,86.109833117236201,86.345785003196601,86.436372872745701,86.6565388023801,86.716562293162301,87.090546367205101,87.294689077778501,87.59770713913501,87.728376568847701,88.026088410725201,88.389798248953201,88.634149509644701,89.909498040679201,89.94992586110501,89.953066027881801,93.718974358974401,93.979697285650501,94.640080105271101,95.026306180707701,96.065588770153401,96.842376205987601,96.988849993517401,97.046152113310601,97.110428231562301,97.184904560687701,97.377312327667701,97.386941016820801,97.389084324115701,97.452238761624601,97.455025392362201,97.459199625966101,97.464125395540501,97.507555434819501,97.571327182398901,97.593232621662301,97.601137238642401,97.681386130738501,97.725588222470901,97.737475078783201,97.790368271954701,97.812809226464701,97.844585678098501,97.851877954952601,97.867187270645601,97.886215161911901,97.895825345370501,97.938353175839901,97.941733480001401,98.009914157901101,98.013162656052101,98.016165569749201,98.047085311536401,98.054077013950601,98.063084742977201,98.086526188503501,98.112578977599101,98.11433633612901,98.120775545524701,98.128114043181401,98.159634421526501,98.166911496980401,98.217218580232201,98.328954851989301,98.344516358886101,98.355512753166401,98.478359036111601,98.80478087649401,99.86595536208301,99.897010407652601,99.925065673176501),2),c(rep(1,105),rep(0,105)),cex=c(as.vector(M[,2]),as.vector(M[,1]))/3,pch=19, xlab="",ylab='',xlim=c(0,100),axes=F)
axis(1)
axis(2, at = seq(0,1,by=1), las = 1)
box()
mtext("Loss",side=2,line=1.5,at=0)
mtext("Diastance between SF3B3 and MTSS2 (in %)",side=1,line=2,at=50,cex=2)
mtext(substitute(paste(italic('IL34'), " gene")),side=2,line=2,at=0.5,cex=2)
mtext("Retention",side=2,line=1.5,at=1,cex=2)
curve(plogis(cc1[1]+cc1[2]*x),col="blue",add=TRUE,lwd=2)
textleg=paste("logistic_MPLE\nIntercept = ",mpleinter,"\nSlope = ",mplew,"\np = 0.01387\nn = 105")
legend(65,0.9,textleg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n',col="blue")
curve(plogis(cc2[1]+cc2[2]*x),col="red",add=TRUE,lwd=2)
text2leg=paste("logistic_IG10\nIntercept = ",ig10inter,"\nSlope = ",ig10w,"\np = 0.01399\nn = 105")
legend(55,0.3,text2leg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n',col="red")
dev.off()







#################
library(phytools)
a=read.table("wir_data_ordered.txt",header=F,row.names=1)
t=read.tree("species24forfigure.nwk")
a$V3 = a$V3/100
a$V4 = a$V4/100
a$V5 = a$V5/100
column_names=c("W","I","R")
i=1
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c=as.data.frame(fit)
colnames(c)=column_names[i]
i=2
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c[ , ncol(c) + 1] <- fit
i=3
b=as.data.frame(a[,i+1])
rownames(b)=rownames(a)
colnames(b)="bv"
b<-setNames(b$bv,rownames(b))
fit<-fastAnc(t,b,vars=TRUE,CI=FALSE)[1]
c[ , ncol(c) + 1] <- fit
colnames(c)=column_names
b=a[,-1]
colnames(b)=column_names
d=rbind(b,c)
ma=as.matrix(d[,])
geo = factor(a$V2)
mycol <- c("red","blue")[geo]
png("wir_phylogeny_plot.png",height=12,width=22,units = "in",res=200)
par(mfrow= c(1, 2),mar=c(3.5,2.5,4,1))
plot(t, tip.color = mycol,edge.width=3,adj=1,cex=1.5)
nodelabels(node=1:47,pie=ma[,],piecol=c("white","pink","red"),cex=0.75,adj=c(-1.75,0.5))
library(phylolm)
a=read.table("wir_data_ordered.txt",header=F)
head(a)
colnames(a)=c("species","gene_status","W","I","R")
rownames(a)=a[,1]
###############MPLE###############
fitwmple=phyloglm(gene_status~W,a,t, method = c("logistic_MPLE"), btol = 10, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc1=coef(fitwmple)
summary(fitwmple)
mpleinter=round(cc1[1],2)
mplew=round(cc1[2],2)
#####################IG10#########
fitwig10=phyloglm(gene_status~W,a,t, method = c("logistic_IG10"), btol = 20, log.alpha.bound = 4,start.beta=NULL, start.alpha=NULL,boot = 2000, full.matrix = TRUE)
cc2=coef(fitwig10)
summary(fitwig10)
ig10inter=round(cc2[1],2)
ig10w=round(cc2[2],2)
####################################
t(table(a$gene_status,a$W))->M
plot(rep(c(0,6.2,11,12,15.6,17,18,24,25,27,37.5,51,60,67,72.3,76,95),2),c(rep(1,17),rep(0,17)),cex=c(as.vector(M[,2]),as.vector(M[,1]))/3,pch=19, xlab="",ylab='',xlim=c(0,100),axes=F)
axis(1)
axis(2, at = seq(0,1,by=1), las = 1)
box()
mtext("Loss",side=2,line=1.5,at=0)
mtext("Percentage of white muscle (FG) fiber",side=1,line=2,at=50,cex=1.4)
mtext(substitute(paste(italic('COA1'), " gene")),side=2,line=2,at=0.5,cex=1.5)
mtext("Retention",side=2,line=1.5,at=1)
curve(plogis(cc1[1]+cc1[2]*x),col="red",add=TRUE,lwd=2)
textleg=paste("logistic_MPLE\nIntercept = ",mpleinter,"\nSlope = ",mplew,"\np = 0.039\nn = 24")
legend(65,0.8,textleg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n')
curve(plogis(cc2[1]+cc2[2]*x),col="red",add=TRUE,lwd=2)
text2leg=paste("logistic_IG10\nIntercept = ",ig10inter,"\nSlope = ",ig10w,"\np = 0.024\nn = 24")
legend(55,0.4,text2leg,xjust = 0.5,yjust = 0.5,x.intersp = 0.2,y.intersp = 0.8,adj = c(0, 0.5),bty='n')
dev.off()
####################

library(phytools)
library(phylolm)
library(ape)
library(ggplot2)
library(dplyr)

# Load tree
tree <- read.tree("gene_distances.sorted.species.no_internodes.nwk")

# Load data
gene_data <- read.table("gene_distances.sorted.tsv", header = TRUE, sep = "\t")

# Ensure species names in data match those in tree
rownames(gene_data) <- gene_data$species  # Assuming the first column is 'Species'

tree <- drop.tip(tree, setdiff(tree$tip.label, gene_data$species))  # Remove species not in data

gene_trait <- gene_data$SF3B3.MTSS2  # Replace 'Trait' with the relevant column for analysis
names(gene_trait) <- gene_data$species
setdiff(names(gene_trait), tree$tip.label)
setdiff(tree$tip.label, names(gene_trait))
a[, c("V4", "V12", "V16")] <- lapply(a[, c("V4", "V12", "V16")], as.numeric)

# Now apply the percentage transformation
a_percent <- a  # Copy original data
a[, c("V4", "V12", "V16")] <- sweep(a[, c("V4", "V12", "V16")], 2, apply(a[, c("V4", "V12", "V16")], 2, max), FUN = "/") * 100
# Check for missing or non-numeric values
summary(gene_trait)
sum(is.na(gene_trait))
gene_trait <- log1p(gene_trait)  # log(1 + x) to handle zeros  
anc_states <- fastAnc(tree, gene_trait)  
summary(tree$edge.length)  
tree$edge.length[tree$edge.length == 0] <- 1e-6  

# Ancestral state reconstruction
anc_states <- fastAnc(tree, gene_trait)
table(gene_data$gene_status, gene_data$SF3B3.MTSS2)
gene_data$SF3B3.MTSS2 <- gene_data$SF3B3.MTSS2 + rnorm(nrow(gene_data), 0, 1e-6)
cor(gene_data[, c("gene_status", "SF3B3.MTSS2")], use = "complete.obs")

# Phylogenetic logistic regression
phylo_logit_model <- phyloglm(gene_status ~ SF3B3.MTSS2, data = gene_data, phy = tree, method = "logistic_MPLE")  # Replace 'Predictor' with actual column
summary(phylo_logit_model)
phylo_logit_model <- phyloglm(
  gene_status ~ SF3B3.MTSS2, 
  data = gene_data, 
  phy = tree, 
  method = "logistic_MPLE", 
  btol = 100000  # Increase boundary tolerance
)
# Plot ancestral state reconstruction results
df <- data.frame(Species = names(anc_states), Ancestral_State = anc_states)
ggplot(df, aes(x = Species, y = Ancestral_State)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
