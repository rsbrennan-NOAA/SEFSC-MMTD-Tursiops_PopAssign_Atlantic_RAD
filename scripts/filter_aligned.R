### alignment stats
library(ggplot2)
dat <- read.csv("analysis/alignment_stats.csv")
dat$sample <- gsub("\\.merged\\.sorted", "", dat$sample)
dflib <- read.csv("merge_summary.csv", header=F)
colnames(dflib) <- c("sample", "lanes", "dir")
hist(dat$map_percent_q20)

all <- merge(dat, dflib, by="sample")



hist(dat$mapped_q20, breaks=30)


ggplot(all, aes(x = total_reads, fill = as.factor(lanes))) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 40) +
  geom_vline(xintercept = quantile(all$total_mapped, c(0.025)), 
             linetype = "dashed", color = "black", size=2) +
  scale_fill_discrete(name = "Lanes") +
  theme_classic(base_size = 14) +
  labs(x = "Total Reads", 
       y = "Count",
       title = "Total # of Reads by Number of Lanes")

# there are some samples originally included that are the wrong species. we knew this ahead of time. Drop them. 

## 41 are from brazil
all$sample[grep("^41",all$sample)]

all2 <- all[grep("^41",all$sample, invert=T),]
nrow(all2)

# another species
all2$sample[grep("^157",all2$sample)]

all3 <- all2[grep("^157",all2$sample, invert=T),]
nrow(all3)

# drop Tadu
all4 <- all3[grep("^Tadu",all3$sample, invert=T),]
nrow(all4)





length(all$sample[which(all$total_reads < quantile(all$total_mapped, c(0.025)))])
# 11
# 399 remain

write.table(file="scripts/bam.list",paste0("/home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/analysis/merged_bams/",
      all$sample[which(all$total_reads >= quantile(all$total_mapped, c(0.025)))],
      ".merged.sorted.bam"), row.names=F, quote=F, col.names =F)
















