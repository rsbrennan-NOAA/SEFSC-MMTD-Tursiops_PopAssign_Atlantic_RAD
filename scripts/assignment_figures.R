library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tidyverse)
library(ggspatial)
library(marmap)
library(RColorBrewer)

library(scatterpie)
library(maps)


micro_19 <- read.csv("MicroAssign_19.csv")
micro_43 <- read.csv("MicroAssign_19.csv")
all_results <- read.table("analysis/structure/assignPOP_replicates.txt", header=T)

rad_table <- all_results %>%
  group_by(ind_id, pred.pop) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

rad_long <- as_tibble(data.frame(ind_id = rad_table$ind_id,
                       X= c(""), 
                       population=rad_table$pred.pop,
                       proportion = rad_table$proportion,
                       group = "RAD"))



micro_19_long <- micro_19 %>%
            select(-lat, -lon) %>%
            pivot_longer( 
                              cols = c(DarkGreen_17, Red_21, Orange_84,LightGreen_14),
                              names_to = "population",
                              values_to = "proportion")
micro_19_long$group <- "micro_19"
micro_43_long <- micro_43 %>%
          select(-lat, -lon) %>%
          pivot_longer(
                              cols = c(DarkGreen_17, Red_21, Orange_84,LightGreen_14 ),
                              names_to = "population",
                              values_to = "proportion")

micro_43_long$group <- "micro_43"

all_dat <- rbind(rbind(rad_long, micro_19_long), micro_43_long)

unique(all_dat$ind_id)

d <- ggplot(all_dat, aes(x = ind_id, y = proportion, fill = population)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Proportion of Assignments") +
  xlab("Individual ID") +
  scale_fill_manual(values = c("LightGreen_14" = "lightgreen",
                               "DarkGreen_17" = "darkgreen",
                               "Red_21" = "firebrick3",
                               "Orange_84" = "orange")) +
  facet_wrap(~group, nrow=3)

d

ggsave("figures/population_assign_Micro_RAD.png", 
       h=5, w=8)


# make plot of mean probability of assignment for each indiv:
d <- ggplot(all_results, aes(x = ind_id, y = assignment_prob, 
                             fill = pred.pop)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Probability of Assignments") +
  xlab("Individual ID") +
  scale_fill_manual(values = c("LightGreen_14" = "lightgreen",
                               "DarkGreen_17" = "darkgreen",
                               "Red_21" = "firebrick3",
                               "Orange_84" = "orange")) +
  ggtitle("Probabilty of assignments: 100 reps")
d

ggsave("figures/population_assign_Probs.png", 
       h=4, w=5)


# 
head(all_results)


summary_stats <- all_results %>%
  group_by(ind_id, pred.pop) %>%
  summarise(
    count = n(),
    median_darkgreen = median(darkgreen),
    median_red = median(red),
    median_orange = median(orange),
    median_lightgreen = median(lightgreen),
    median_prob = median(assignment_prob)
  ) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# find confident assignments:
high_assigns <- summary_stats %>%
  group_by(ind_id) %>%
  slice_max(order_by = proportion, n = 1) %>%
  ungroup()

nrow(high_assigns)
high_assigns$assigned <- FALSE
high_assigns$assigned[which(high_assigns$median_prob > 0.6)] <- TRUE

as.data.frame(high_assigns)



#pltpt2 <- cbind(pltpt, qval)
wide_summary <- summary_table %>%
  pivot_wider(
    id_cols = ind_id,
    names_from = pred.pop,
    values_from = proportion,
    values_fill = 0  # Fill NAs with 0 since no occurrence means 0 proportion
  )

# get coords.

coords <- read.csv("Tursiops_RADseq_Metadata.csv")
coords2 <- data.frame(ID = coords$Lab.ID,
                      lat = coords$Lat,
                      lon = coords$Long)

pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))
nrow(pops)
sum(pops$indiv %in% coords2$ID)

m1 <- merge(coords2, pops, by.x="ID", by.y="indiv")
nrow(m1)
wide_summary$ID <- as.factor(wide_summary$ind_id)

wide_summary$ID[!wide_summary$ID %in% coords2$ID]

unknowns <- merge(coords2, wide_summary, by="ID")
unknowns

unknowns2 <- merge(unknowns, data.frame(ID = high_assigns$ind_id, pred.pop=high_assigns$pred.pop, assigned = high_assigns$assigned), by="ID")
unknowns <- unknowns2 %>%
  select(-ind_id )


m1$lon <- as.numeric(m1$lon)

high_assigns$assigned
world <- ne_countries(scale = "medium", returnclass = "sf")
usa <- st_as_sf(maps::map("state", fill=TRUE, plot =FALSE))



unknowns$lat <- as.numeric(unknowns$lat)
unknowns$lon <- as.numeric(unknowns$lon)
unknowns$pred.pop <- ifelse(unknowns$assigned == FALSE, "No confidence", unknowns$pred.pop)
unknowns$pred.pop <- as.factor((unknowns$pred.pop))
unknowns <- unknowns %>%
  select(-assigned )
unknowns$group <- "RAD"

head(unknowns)
# 
micro_19_2 <- micro_19 %>%
  select(ind_id, lat, lon, DarkGreen_17, Red_21, Orange_84) %>%
  rename(ID = ind_id) %>%
  mutate(pred.pop = case_when(
    DarkGreen_17 >= Red_21 & DarkGreen_17 >= Orange_84 ~ "DarkGreen_17",
    Red_21 >= DarkGreen_17 & Red_21 >= Orange_84 ~ "Red_21",
    Orange_84 >= DarkGreen_17 & Orange_84 >= Red_21 ~ "Orange_84"
  ))

micro_19_2$group <- "micro_19"

micro_43_2 <- micro_43 %>%
  select(ind_id, lat, lon, DarkGreen_17, Red_21, Orange_84) %>%
  rename(ID = ind_id) %>%
  mutate(pred.pop = case_when(
    DarkGreen_17 >= Red_21 & DarkGreen_17 >= Orange_84 ~ "DarkGreen_17",
    Red_21 >= DarkGreen_17 & Red_21 >= Orange_84 ~ "Red_21",
    Orange_84 >= DarkGreen_17 & Orange_84 >= Red_21 ~ "Orange_84"
  ))

micro_43_2$group <- "micro_43"

pltall <- rbind(rbind(unknowns, micro_19_2), micro_43_2)

p <- 
  
  ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey70") +
  geom_sf(data = usa, fill = NA, color = "grey70") +
  geom_point(data =m1, aes(x=lon, y=lat, color=as.factor(pop)),
             shape= 17, size = 1.5,
             alpha=1) +
  geom_scatterpie(data = pltall, 
                  aes(x = lon, y = lat), 
                  cols = c("Orange_84", "Red_21", "DarkGreen_17"), 
                  alpha = 1, 
                  pie_scale=5,
                  color = "black",
                  linewidth = 0.1) +
  scale_fill_manual(values=c("orange3", "red3", "darkgreen")) +
  scale_color_manual(values=c("darkgreen","lawngreen","orange3", "red3")) +
  theme_bw() +
  theme(
    #panel.background = element_rect(fill = "white"),
    panel.grid = element_blank(),
    #axis.text = element_blank(),
    #axis.ticks = element_blank()
    legend.position = "top",
    legend.title=element_blank()) +
  xlab("Longitude")+
  ylab("Latitude") +
  coord_sf(xlim = c(-83, -70), ylim = c(30, 41), expand = FALSE) +
  facet_grid(group ~ pred.pop)

p

ggsave("figures/assignment_map_RAD.png", p, h=8, w=8)



#---------
# micros
### 19

### 43




















# add into PCA
#---------------------------------------------------------------------------------

keep <- c(m1$ID, unknowns$ID)
length(keep)
pops <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

filename = "analysis/filtered.final_ids"
filename.gds = paste0(filename, ".gds")
filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
#SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

gdsin = SeqArray::seqOpen(filename.gds)

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
known <- merge(dat, m1, by.x="IDs", by.y="ID")
unknown_df <- merge(dat, unknowns, by.x="IDs", by.y="ID")
nrow(dat)

unknown_df_lowconf <- unknown_df[unknown_df$assigned == FALSE, ]


d <- ggplot(known, aes(PC1, PC2, color=pop)) +
  geom_point(size=4, pch=17, alpha=0.7) + 
  #geom_point(data=unknown_df,  aes(PC1, PC2, fill=pop),size=4, pch=14, fill="grey") +
  geom_scatterpie(data = unknown_df, 
                  aes(x = PC1, y = PC2), 
                  cols = c("Orange_84", "Red_21", "DarkGreen_17"), 
                  alpha = 1, 
                  pie_scale=3,
                  color = "black",
                  linewidth = 0.1) +
  geom_point(data=unknown_df_lowconf,  aes(PC1, PC2),size=4,stroke=2, pch=4,color="black") +
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Atlantic with unknown: PC1, PC2') +
  scale_color_manual(values=c("darkgreen","lawngreen","orange3", "red3"))+
  scale_fill_manual(values=c("orange3","red3","darkgreen", "lawngreen"))
d


ggsave("figures/PCA_Atl_assignments_1_2.png",
       d, w=7, h=5)
