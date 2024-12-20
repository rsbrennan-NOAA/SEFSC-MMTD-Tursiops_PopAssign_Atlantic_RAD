# first plot the plink result:
library(ggplot2)
library(dplyr)

dat <- read.table("analysis/structure/all_indivs_PCA.eigenvec", header=F)
eigenval <- read.table("analysis/structure/all_indivs_PCA.eigenval", header=F)

# first convert to percentage variance explained
pve <- data.frame(PC=1:20, pve=round(eigenval$V1/sum(eigenval$V1)*100,1))

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

# plot the PC's
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + 
  ylab("Percentage variance explained") + 
  theme_classic()
a

ggsave("figures/pca_allindivs_var.png",
       a, w=4, h=4)


####################
# plot the PCA
####################

# rename our columns, just for ease
colnames(dat) <- c("ID", "ID2", "PC1", "PC2", "PC3", "PC4", colnames(dat)[7:ncol(dat)])

# plot the PCA
dat$population <- substr(dat$ID, 1,2)
d <- ggplot(dat, aes(PC1, PC2, label=ID, color=population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC2: ",pve$pve[2],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC2")
d

ggsave("figures/PCA_1_2_plink.png",
       d, w=6, h=4)


d <- ggplot(dat, aes(PC1, PC3, label=ID, color=population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",pve$pve[1],"% variance")) +
  ylab(paste0("PC3: ",pve$pve[3],"% variance")) +
  theme_bw() +
  ggtitle("plink: PC1, PC3")
d

ggsave("figures/PCA_1_3_plink.png",
       d, w=6, h=4)


####---------------------------------------------------------------------------
# snpRelate pca

library(SNPRelate)
library(SeqArray)
setwd("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/analysis")

pops <- read.csv("../SW_Metadata.csv")

# try with missing data allowed. then with no missing data. 
filename = "filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)
print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))

print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))

summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))
samples <- SeqArray::seqGetData(gdsin, "sample.id")
cbind(samples,m2)[order(-m2),]
hist(m2,breaks=50)

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05) 

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- gsub("b", "",pca.out$sample.id)

result <- dat %>%
  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location, Sex), by = c("IDs" = "Lab.ID.."))

# with ids
d <- ggplot(result, aes(PC1, PC2, label=IDs)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('snpRelate: PC1, PC2')
d

# with points
d <- ggplot(result, aes(PC1, PC2, fill=Pop.Structure.Location)) +
  geom_point(size =2, color="black", shape=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('snpRelate: PC1, PC2')
d

ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_2_snpRelate.png",
       d, w=4, h=4)


#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=IDs)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('snpRelate: PC1, PC3')
d
ggsave("C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/figures/PCA_1_3_snpRelate.png",
       d, w=4, h=4)