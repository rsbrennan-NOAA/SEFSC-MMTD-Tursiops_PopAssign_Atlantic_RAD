library(tidyverse)
dat <- read.csv("analysis/out.imiss", 
                header=T, sep="\t")

hist(dat$F_MISS, breaks=30)

dat[which(dat$F_MISS > 0.4),]
nrow(dat[which(dat$F_MISS > 0.4),])

write.table(file="analysis/rm_missing.txt",data.frame(d=dat$INDV[which(dat$F_MISS > 0.4)]), col.names=F, 
  row.names=F, quote=F)




###### post filtering

dat <- read.csv("analysis/filtered.6.depth.ldepth.mean", 
                header=T, sep="\t")

hist(dat$MEAN_DEPTH, breaks=40)
hist(dat$MEAN_DEPTH, breaks=60, xlim=c(0, 500))

mean(dat$MEAN_DEPTH)
# 28.34
# *3 = 85.02

quantile(dat$MEAN_DEPTH, probs=.99)
#138
sum(dat$MEAN_DEPTH > 85.02)
dat[which(dat$MEAN_DEPTH > 85.02),]

sum(dat$MEAN_DEPTH > quantile(dat$MEAN_DEPTH, probs=.99))
dat[which(dat$MEAN_DEPTH > quantile(dat$MEAN_DEPTH, probs=.99)),]

# HDplot
library(vcfR)

vcfInput<-read.vcfR("analysis/filtered.6.vcf.gz")

source("scripts/HDplot.R")

HDplotResults<-HDplot(vcfInput)

head(HDplotResults)
mean(HDplotResults$H)

hist(HDplotResults$num_hets)

HDplotResults %>% ggplot()+geom_point(aes(x=H,y=D))

#plot H and ratio
HDplotResults %>% ggplot()+geom_point(aes(x=H,y=ratio))


# sites beyond the distribution are likely paralogs
## D- > than 10 for sure, but maybe even 5?
## H - > than about 0.6


sum((HDplotResults$H > 0.6))
#123
sum((abs(HDplotResults$D) > 6))
#445

# positions to exclude:
datexclude <- HDplotResults[which(HDplotResults$H > 0.6 | abs(HDplotResults$D) > 6),]
posexclude <- datexclude[,1:2]
nrow(posexclude)
#498
write.table(posexclude, file="C:/Users/Reid.Brennan/Documents/projects/spermWhaleRad/HD_exclude.txt",
            quote=F, col.names = FALSE, row.names=FALSE)