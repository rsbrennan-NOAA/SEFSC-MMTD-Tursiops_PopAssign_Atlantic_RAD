# first plot the plink result:
library(ggplot2)
library(dplyr)
library(SNPRelate)
library(SeqArray)
library(stringr)

####---------------------------------------------------------------------------
# snpRelate pca

#pops <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

filename = "analysis/filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)

print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))
# [1] "The number of SAMPLES in data: 380"
print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))
# [1] "The number of SNPs in data: 8759"
summary(m1 <- SeqArray::seqMissing(gdsin, per.variant=TRUE))
summary(m2 <- SeqArray::seqMissing(gdsin, per.variant=FALSE))

samples <- SeqArray::seqGetData(gdsin, "sample.id")

pop.ids <- str_extract(samples, "^[^Tt]+")
pop.ids <- sub("^(FB).*$", "\\1", pop.ids)


#--------------------------
### all samples
snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, 
                                     missing.rate=NaN,
                                     remove.monosnp = T)  
# 4,691  markers are selected in total.
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id) # filtering for pruned SNPs


eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- pop.ids

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC2')
d

ggsave("figures/PCA_allIndivs_1_2.png",
       d, w=7, h=5)


# with ids
d <- ggplot(dat, aes(PC1, PC3, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('All individuals: PC1, PC3')
d

ggsave("figures/PCA_allIndivs_1_3.png",
       d, w=7, h=5)


#--------------------------
### drop brazil samples, 41Tt
keep <- samples[!pop.ids == "41"]
# 356 samples
snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate=NaN
                                  
                                     )  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Drop Brazil: PC1, PC2')
d

ggsave("figures/PCA_dropBrazil_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=IDs, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Drop Brazil: PC1, PC3')
d

ggsave("figures/PCA_dropBrazil_1_3.png",
       d, w=7, h=5)

#------------------------------------------------------------------------------
# which are definitely atlantic samples? color those:
# I think from that table that Nikki sent. 

known_samps <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

atl_samps <- known_samps$Lab_ID

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)
dat$Population[dat$IDs %in% atl_samps] <- "Atlantic"

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_point(data=subset(dat, Population == "Atlantic"),color="grey5") +
  geom_text(data=subset(dat, Population != "Atlantic"), size=3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Drop Brazil: PC1, PC2')
d

ggsave("figures/PCA_dropBrazil_withAtlanticPts_1_2.png",
       d, w=7, h=5)


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# drop the outliers from PC3, they're from south atlantic
# all start with 157
# also drop 42193 and 78068 # why??? I forget

keep2 <- keep[grep("157", keep, invert=T)] # these are some other species
length(keep2)
keep3 <- keep2[grep("42193|78068", keep2, invert=T)] # from somewhere else. San diego?
length(keep3)# 350

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep3,
                                     missing.rate=NaN
                                     )  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep3) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Drop S. Atlantic: PC1, PC2')
d

ggsave("figures/PCA_dropSAtlantic_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=IDs, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Drop S. Atlantic: PC1, PC3')
d

ggsave("figures/PCA_dropSAtlantic_1_3.png",
       d, w=7, h=5)

#------------------------------------------------------------------------------
# which are definitely atlantic samples? color those:
# I think from that table that Nikki sent. 

known_samps <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

atl_samps <- known_samps$Lab_ID

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)
dat$Population[dat$IDs %in% atl_samps] <- "Atlantic"

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_point(data=subset(dat, Population == "Atlantic"),color="grey5") +
  geom_text(data=subset(dat, Population != "Atlantic"), size=3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Black = Atlantic Samples: PC1, PC2')
d

ggsave("figures/PCA_dropSAtlantic_withAtlanticPts_1_2.png",
       d, w=7, h=5)


# with ids
d <- ggplot(dat, aes(PC1, PC3, label=Population, color=Population)) +
  geom_point(data=subset(dat, Population == "Atlantic"),color="grey5") +
  geom_text(data=subset(dat, Population != "Atlantic"), size=3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Black = Atlantic Samples: PC1, PC3')
d

ggsave("figures/PCA_dropSAtlantic_withAtlanticPts_1_3.png",
       d, w=7, h=5)


#--------------------------------------------------------------------------------------
### subset to just the atlantic nearshore?
samples <- SeqArray::seqGetData(gdsin, "sample.id")
pop.ids <- str_extract(samples, "^[^Tt]+")
pop.ids <- sub("^(FB).*$", "\\1", pop.ids)

sum(dat$PC1 > 0 & dat$PC2 < 0)
#1711
keep <- dat$IDs[dat$PC1 > 0 & dat$PC2 < 0]

length(keep)

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, 
                                     remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate=NaN)  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Nearshore Atlantic: PC1, PC2')
d

ggsave("figures/PCA_NearshoreAtl_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=IDs, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Nearshore Atlantic: PC1, PC3')
d

ggsave("figures/PCA_NearshoreAtl_1_3.png",
       d, w=7, h=5)


##--------------------------------------------------------------------------------------
### drop the 44 samples
samples <- SeqArray::seqGetData(gdsin, "sample.id")
pop.ids <- str_extract(samples, "^[^Tt]+")
pop.ids <- sub("^(FB).*$", "\\1", pop.ids)

#
keep <- dat$IDs[dat$PC1 > 0 & dat$PC2 < 0]
keep2 <- keep[grep("44",keep, invert=T)]
length(keep2)

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep2,
                                     missing.rate=NaN
                                     )  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep2) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Nearshore Atlantic, no 44: PC1, PC2')
d

ggsave("figures/PCA_NearshoreAtlNo44_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=IDs, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Nearshore Atlantic, no 44: PC1, PC3')
d

ggsave("figures/PCA_NearshoreAtlNo44_1_3.png",
       d, w=7, h=5)


#-------------------------------------------------------------------------------
# just Nikki's atlantic samples
#--------------------------------------------------------------------------------------
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))

filename = "analysis/filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
#SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

#gdsin = SeqArray::seqOpen(filename.gds)

# only the atlantic samples from Nikki
keep <- pops$indiv

length(keep)

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate=NaN)  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
out <- merge(x=dat, y=pops, by.x="IDs", by.y="indiv")
dat<- out
#dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, fill=pop)) +
  geom_point(size=4, pch=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw(base_size = 14) +
  ggtitle('Nikki NC Pops: PC1, PC2') +
  scale_fill_manual(values=c("darkgreen","lawngreen","orange3", "red3"))
d

ggsave("figures/PCA_NikkiLabels_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, fill=pop)) +
  geom_point(size=4, pch=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Nikki NC Pops: PC1, PC3') +
  scale_fill_manual(values=c("darkgreen","lawngreen","orange3", "red3"))

d

ggsave("figures/PCA_NikkiLabels_1_3.png",
       d, w=7, h=5)




#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# just  atlantic samples, but use nikki's known assignments.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


pops <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

filename = "analysis/filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
#SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)

print(paste0("The number of SAMPLES in data: ", length(c(SeqArray::seqGetData(gdsin, "sample.id")))))

print(paste0("The number of SNPs in data: ",  length(c(SeqArray::seqGetData(gdsin, "variant.id")))))

# only the atlantic samples from Nikki
keep <- pops$Lab_ID

length(keep)

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate= NaN)  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
dat$Population <- str_extract(dat$IDs, "^[^Tt]+")
dat$Population <- sub("^(FB).*$", "\\1", dat$Population)

# with ids
d <- ggplot(dat, aes(PC1, PC2, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[2], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic Only: PC1, PC2')
d

ggsave("figures/PCA_AtlOnlyKnown_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, label=Population, color=Population)) +
  geom_text(size =3) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic Only: PC1, PC3')
d

ggsave("figures/PCA_AtlOnlyKnown_1_3.png",
       d, w=7, h=5)


SeqArray::seqGetData(gdsin, "$chrom_pos")
SnpLoad <- snpgdsPCASNPLoading(pca.out, gdsin)
pc<-1
hist(sort(abs(SnpLoad$snploading[pc,]),decreasing=T,index.return=T)[[1]],breaks = 30,main="PC loadings: PC1",)

high.snp.ids <- SnpLoad$snp.id[which(abs(SnpLoad$snploading[1,]) >= 0.035)]
SnpLoad



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
### subset to just the bottom left samples
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

samples <- SeqArray::seqGetData(gdsin, "sample.id")
pop.ids <- str_extract(samples, "^[^Tt]+")
pop.ids <- sub("^(FB).*$", "\\1", pop.ids)

sum(dat$PC1 < 0 & dat$PC2 < 0)
#145
keep <- dat$IDs[dat$PC1 < 0 & dat$PC2 < 0]

length(keep)

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate=NaN)  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )

# add in the labels based on nikki assignments
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))
dat$IDs <- pca.out$sample.id
out <- merge(x=dat, y=pops, by.x="IDs", by.y="indiv", all.x=T)
dat<- out


d <- ggplot(dat, aes(PC1, PC2, fill=pop)) +
  geom_point(size=4, pch=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic with unknown: PC1, PC2') +
  scale_fill_manual(values=c("darkgreen","lawngreen","orange3", "red3"))
d

head(dat)
dat[dat$PC2< -0.3,]

# drop the 44Tt samples:

keep2 <- keep[grep("44Tt", keep, invert=T)]
keep <- keep2

snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
                                     start.pos="random", num.thread=1, remove.monosnp = T, 
                                     sample.id = keep,
                                     missing.rate=NaN)  
snpset.id <- unlist(unname(snpset))

pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = T, maf = 0.05,
                               snp.id=snpset.id,
                               sample.id = keep) # filtering for pruned SNPs

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )

# add in the labels based on nikki assignments
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))
dat$IDs <- pca.out$sample.id
out <- merge(x=dat, y=pops, by.x="IDs", by.y="indiv", all.x=T)
dat<- out


d <- ggplot(dat, aes(PC1, PC2, fill=pop)) +
  geom_point(size=4, pch=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic with unknown: PC1, PC2') +
  scale_fill_manual(values=c("darkgreen","lawngreen","orange3", "red3"))
d


ggsave("figures/PCA_AtlWithUnknown_1_2.png",
       d, w=7, h=5)

#pc1 and 3
d <- ggplot(dat, aes(PC1, PC3, fill=pop)) +
  geom_point(size=4, pch=21) +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC3: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic with unknown: PC1, PC3') +
  scale_fill_manual(values=c("darkgreen","lawngreen","orange3", "red3"))

d

ggsave("figures/PCA_AtlWithUnknown_1_3.png",
       d, w=7, h=5)


### save indiv ids to use with admixture.
samples <- SeqArray::seqGetData(gdsin, "sample.id")

write.table(file= "analysis/coastalAtl_unknowns.txt", samples,
            row.names=F, quote = F, col.names = F)

# and subset down the orange population
head(dat)
samples2 <- dat[,c("IDs", "pop")]
samples3 <- samples2[which(samples2$pop != "Orange_84"),]
samples4 <- samples2[is.na(samples2$pop),]
samples5 <- rbind(samples3, samples4)

allOrange <- samples2[which(samples2$pop == "Orange_84"),]
sampOrange <- allOrange[sample(nrow(allOrange), 20, replace=FALSE),]

sampout <- rbind(samples5, sampOrange)
nrow(sampout)
write.table(file= "analysis/coastalAtl_orangeSubset_unknowns.txt", sampout$IDs,
            row.names=F, quote = F, col.names = F)







