### alignment stats
library(ggplot2)
dat <- read.csv("alignment_stats.csv")
dat$sample <- gsub("\\.merged\\.sorted", "", dat$sample)
dflib <- read.csv("merge_summary.csv", header=F)
colnames(dflib) <- c("sample", "lanes", "dir")
hist(dat$map_percent_q20)

all <- merge(dat, dflib, by="sample")



hist(dat$mapped_q20, breaks=30)


table(dflib$V2)
dflib[which(dflib$V2 == 4),]


ggplot(all, aes(x = total_reads, fill = as.factor(lanes))) +
  geom_histogram(position = "identity", alpha = 0.9, bins = 40) +
  geom_vline(xintercept = quantile(all$total_mapped, c(0.025)), 
             linetype = "dashed", color = "black", size=2) +
  scale_fill_discrete(name = "Lanes") +
  theme_classic(base_size = 14) +
  labs(x = "Total Reads", 
       y = "Count",
       title = "Total # of Reads by Number of Lanes")

length(all$sample[which(all$total_reads < quantile(all$total_mapped, c(0.025)))])
# 11
# 399 remain

write.table(file="scripts/bam.list",paste0("/home/rbrennan/Tursiops-NC-PopulationAssignment-RAD/analysis/merged_bams/",
      all$sample[which(all$total_reads >= quantile(all$total_mapped, c(0.025)))],
      ".merged.sorted.bam"), row.names=F, quote=F)



#------------------------------------------------------------------------------
## make a map to figure out where these samples are:
## in case they're really clustered or something














