
library(vcfR)
library(adegenet)
library(ggpubr)
library(dplyr)
library(poppr)

vcf <- read.vcfR( "data/ATLOnly_filtered_LDthin_ChrCorrected.vcf", verbose = FALSE )

genl<-vcfR2genlight(vcf)

# add population ids
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))

genl@ind.names <- gsub("b", "",genl@ind.names)

ids <- data.frame(IDs = genl@ind.names)
#result <- ids %>%
#  left_join(pops %>% select( Lab.ID.., Pop.Structure.Location), by = c("IDs" = "Lab.ID.."))

#genl@pop <- as.factor(result$Pop.Structure.Location)

## Full data set
# re-running pca, basically to make sure consistent again with plink, etc. 
pca1 <- glPca(genl,center = T, scale = T, nf = 5)

#png("../figures/dapc_eigen.png", h=4, w=4, units="in", res=300)
barplot(100*pca1$eig/sum(pca1$eig), col = heat.colors(50), main="PCA Eigenvalues") # retain first 5 axes, incremental decrease after 2
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
#dev.off()

#proportion of explained variance by first three axes
a1<-pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis
a2<-pca1$eig[2]/sum(pca1$eig) # proportion of variation explained by 2nd axis 
a3<-pca1$eig[3]/sum(pca1$eig) # proportion of variation explained by 3rd axis
pcvar <- data.frame(Axis = c(1:3), Proportion = c(a1,a2,a3))
pcvar

pca1.scores <- as.data.frame(pca1$scores)
pca1.scores$pop <- pop(genl)
pca1.scores$ind <- genl@ind.names

set.seed(89)
num_pops <- length(levels(factor(pca1.scores$pop)))

pcplt <- merge(pca1.scores, pops, by.x="ind", by.y="indiv", all.x=T)
# plot PC 1 and 2
pca1.p<-ggscatter(pcplt, x = "PC1", y = "PC2", color = "pop",
                  ellipse = T, ellipse.level = 0.95, size = 3,
                  xlab = paste0("PC1 (",round(pcvar[1,2]*100,2),"%)"),
                  ylab = paste0("PC2 (",round(pcvar[2,2]*100,2),"%)"),
                  palette =c("darkgreen","lawngreen","orange3", "red3")
)
pca1.p
#ggsave("../figures/dapc_PCA.png",pca1.p, h=5, w=5)



###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------
# run dapc
###------------------------------------------------------------------------
###------------------------------------------------------------------------
###------------------------------------------------------------------------

# follow the k-1 recommendation for the number of PCAs. So k=4
grp <- find.clusters(genl, max.n.clust=40)

# keep k-1: 3 
# lowest BIC is 11 clusters- plateaus around 11. seems... incorrect.
# do 4 instead, based on biology
dapc1 <- dapc(genl, grp$grp, n.pca=20, n.da=4)

#temp <- optim.a.score(dapc1)
col.in <- c("#E69F00","#56B4E9", "#009E73", "#CC79A7", "red")
scatter.dapc(dapc1)
scatter.dapc(dapc1, grp=genl@pop)

#png(file="../figures/dapc.4.png", h=4, w=4, units="in", res=300)
scatter(dapc1, grp=grp$grp, 
        bg="white", pch=c(16,15,18,17, 16), col=col.in,
        cex=2,
        cstar=0,clab=0,
        scree.pca=FALSE,scree.da=FALSE, leg=TRUE, posi.leg="bottomleft")

#dev.off()
