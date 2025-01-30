# rubias
# https://cran.r-project.org/web/packages/rubias/vignettes/rubias-overview.html
library(rubias)
library(adegenet)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(vcfR)


# read in vcf
vcf <- read.vcfR( "data/ATLOnly_filtered_LDthin_ChrCorrected.vcf", 
                  verbose = FALSE )

gen <- vcfR2genind(vcf)
geno <- genind2df(gen)
data <- matrix(NA, nrow=nInd(gen), ncol=2*nLoc(gen))
colnames(data) <- paste0(rep(locNames(gen), each=2), c("","_2"))
rownames(data) <- indNames(gen)

for(i in 1:ncol(geno)) {
  locus <- geno[,i]
  data[,2*i-1] <- substr(locus, 1, 1)  
  data[,2*i] <- substr(locus, 2, 2)
  # Handle NAs
  na_idx <- is.na(locus)
  data[na_idx, 2*i-1] <- NA
  data[na_idx, 2*i] <- NA
}
data <- as.data.frame(data)

# get populations for the samples
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))

result <- data %>%
  rownames_to_column("indiv") %>%
  left_join(pops, by = "indiv") %>%
  select(indiv, pop)

# change to rubias format.
df <- cbind(data.frame(sample_type = rep(NA, nrow(data)),
                       repunit = as.character(result$pop),
                       collection = as.character(result$pop),
                       indiv = result$indiv),
                      data
)


# pull out reference and mixture samples:

ref_ids <- !is.na(df$repunit)
mix_ids <- is.na(df$repunit)

ref_df <- df[ref_ids,]
mix_df <- df[mix_ids,]
mix_df$collection <- "unknown"
mix_df$sample_type <- "mixture"
ref_df$sample_type <- "reference"

# subsample reference orange indivs:
ref_noOrange <- ref_df[ref_df$collection != "Orange_84",]
ref_OnlyOrange <- ref_df[ref_df$collection == "Orange_84",]
ref_subOrange <- ref_OnlyOrange[sample(1:nrow(ref_OnlyOrange), 20),]
nrow(ref_subOrange)

ref_subsample <- rbind(ref_noOrange, ref_subOrange)

mix_est <- infer_mixture(reference = ref_subsample, 
                         mixture = mix_df, 
                         gen_start_col = 5)

mix_est

# for mixing proportions
rep_mix_ests <- mix_est$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi))  # adding mixing proportions over collections in the repunit

rep_indiv_ests <- mix_est$indiv_posteriors %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))

rep_indiv_ests

# get the maximum-a-posteriori population for each individual
map_rows <- mix_est$indiv_posteriors %>%
  group_by(indiv) %>%
  top_n(1, PofZ) %>%
  ungroup()

table(map_rows$collection)

mix_est$indiv_posteriors


normo <- tibble(z_score = rnorm(1e06))
ggplot(map_rows, aes(x = z_score)) +
  geom_density(colour = "blue") +
  geom_density(data = normo, colour = "black")

# are any individuals low confidence?
# The z_score statistic is most useful as a check for individuals. 
# It is intended to be a quick way to identify aberrant individuals. 
# If you see a z-score to the maximum-a-posteriori population for an individual 
# in your mixture sample that is considerably less than z_scores you saw in the reference, 
# then you might infer that the individual doesn’t actually fit any of the populations in the reference well.

map_rows[order(map_rows$z_score, decreasing=F),]
# A tibble: 19 × 10
#mixture_collection indiv   repunit  collection  PofZ log_likelihood z_score n_non_miss_loci
#<chr>              <chr>   <chr>    <chr>      <dbl>          <dbl>   <dbl>           <int>
# 1 unknown            2Tt605  DarkGre… DarkGreen… 1.00          -2568.  -3.16             3077
# 2 unknown            2Tt485  Red_21   Red_21     1.00          -2215.  -1.83             2575
# 3 unknown            14Tt386 Orange_… Orange_84  1             -2113.  -1.81             2452
# 4 unknown            4Tt830  Orange_… Orange_84  1             -2144.  -1.06             2518

mix_est$indiv_posteriors[mix_est$indiv_posteriors$indiv == "2Tt605",]
mix_est$indiv_posteriors[mix_est$indiv_posteriors$indiv == "FB712",]
mix_est$indiv_posteriors$PofZ[mix_est$indiv_posteriors$indiv == "FB712"]

# the z-scores are based on the sum of squared log-likelihood of each individual's 
   # genotype, given the allele counts in each collection. 

mix_est$indiv_posteriors$log_likelihood[mix_est$indiv_posteriors$indiv == "2Tt605"]





#-------------------------------------------------------------------------------
# assignPOP
library(assignPOP)
library(dartR)

#needs either structure or genepop format.

dfgl <- gl.read.vcf("data/ATLOnly_filtered_LDthin_ChrCorrected.vcf", verbose = NULL)
pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))

known_gl <- dfgl[indNames(dfgl) %in% pops$indiv]
unknown_gl <- dfgl[!indNames(dfgl) %in% pops$indiv]

# Create a named vector for lookup
pop_lookup <- setNames(pops$pop, pops$indiv)

# Match populations while keeping original order
matched_pops <- data.frame(
  indiv = known_gl@ind.names,
  pop = pop_lookup[known_gl@ind.names],
  stringsAsFactors = FALSE
)

known_gl@pop <- as.factor(matched_pops$pop)

gl2structure(known_gl, outfile = "ATLOnly.str",
             addcolumns=known_gl@pop,
             outpath = "data")
gl2structure(unknown_gl, outfile = "unknown.str",
             addcolumns=unknown_gl@pop,
             outpath = "data")

known_dat <- read.Structure("./data/ATLOnly.str", 
                      ploidy=2)
unknown_dat <- read.Structure("./data/unknown.str", 
                            ploidy=2)

assign.X( x1=known_dat, x2=unknown_dat, 
          dir="assignPOP_results/", model="svm",
          mplot=F)

results <- read.table("assignPOP_results/AssignmentResult.txt",
                      header=T, sep=" ")

results <- unique(results)

data_long <- pivot_longer(results, 
                          cols = c(DarkGreen_17, Red_21, Orange_84,LightGreen_14),
                          names_to = "Population",
                          values_to = "Proportion")


ggplot(data_long, aes(x = Ind.ID, y = Proportion, fill = Population)) +
  geom_col() +
  scale_fill_manual(values = c("LightGreen_14" = "lightgreen",
                                "DarkGreen_17" = "darkgreen",
                               "Red_21" = "firebrick3",
                               "Orange_84" = "orange")) +
  theme_bw() +
  labs(x = "Individual ID",
       y = "Probability") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("population assignment- all Orange")

ggsave("figures/population_assign_allOrange.png", 
       h=4, w=5)

#---------------------------------------------------------------------
# subset orange individuals


known_gl <- dfgl[indNames(dfgl) %in% pops$indiv]
unknown_gl <- dfgl[!indNames(dfgl) %in% pops$indiv]

# Create a named vector for lookup
pop_lookup <- setNames(pops$pop, pops$indiv)

# Match populations
matched_pops <- data.frame(
  indiv = known_gl@ind.names,
  pop = pop_lookup[known_gl@ind.names],
  stringsAsFactors = FALSE
)

known_gl@pop <- as.factor(matched_pops$pop)

# Get indices of populations to keep 
keep_idx <- which(known_gl@pop %in% c("DarkGreen_17", "Red_21", "LightGreen_14"))

# Orange pop
orange_idx <- which(known_gl@pop == "Orange_84")

# subsample orange
orange_sample <- sample(orange_idx, 14)

# Combine and subset object
final_idx <- sort(c(keep_idx, orange_sample))
new_gl <- known_gl[final_idx,]


gl2structure(new_gl, outfile = "ATLOnly_sub.str",
             addcolumns=new_gl@pop,
             outpath = "data")
gl2structure(unknown_gl, outfile = "unknown.str",
             addcolumns=unknown_gl@pop,
             outpath = "data")

known_dat <- read.Structure("./data/ATLOnly_sub.str", 
                            ploidy=2)
unknown_dat <- read.Structure("./data/unknown.str", 
                              ploidy=2)
length(unknown_dat$LocusName)
length(known_dat$LocusName)

sum(unknown_dat$LocusName == known_dat$LocusName)

assign.X( x1=known_dat, x2=unknown_dat, 
          dir="assignPOP_results_sub/", model="svm",
          mplot=F)

results <- read.table("assignPOP_results_sub/AssignmentResult.txt",
                      header=T, sep=" ")

results <- unique(results)

data_long <- pivot_longer(results, 
                          cols = c(DarkGreen_17, Red_21, Orange_84,LightGreen_14),
                          names_to = "Population",
                          values_to = "Proportion")


ggplot(data_long, aes(x = Ind.ID, y = Proportion, fill = Population)) +
  geom_col() +
  scale_fill_manual(values = c("LightGreen_14" = "lightgreen",
                               "DarkGreen_17" = "darkgreen",
                               "Red_21" = "firebrick3",
                               "Orange_84" = "orange")) +
  theme_bw() +
  labs(x = "Individual ID",
       y = "Probability") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("population assignment- subsampled Orange Pop")


ggsave("figures/population_assign_subsetOrange.png", 
       h=4, w=5)

#------------------------------
# add pca using these individuals
## to demonstrate how subsetting changes results:
#data_long
# results
keep <- c(known_dat$SampleID,unknown_dat$SampleID)
length(keep)
#pops <- read.csv("scripts/ATL_RADseq_samples.txt", sep="\t")

#filename = "analysis/filtered.final_ids"
#filename.gds = paste0(filename, ".gds")
#filename.vcf.gz = paste0(filename, ".vcf.gz")
# Convert VCF to GDS
#SeqArray::seqVCF2GDS(vcf.fn = filename.vcf.gz, out.fn = filename.gds, storage.option="ZIP_RA")

#gdsin = SeqArray::seqOpen(filename.gds)

length(keep)

#snpset <- SNPRelate::snpgdsLDpruning(gdsin, ld.threshold=0.2, autosome.only = F, 
#                                     start.pos="random", num.thread=1, remove.monosnp = T, 
#                                     sample.id = keep,
#                                     missing.rate= NaN)  
snpset.id <- SeqArray::seqGetData(gdsin, "variant.id")
#snpset.id.2 <- 
  
pca.out = SNPRelate::snpgdsPCA(autosome.only = F, gdsin, num.thread=2, 
                               remove.monosnp = F, maf = 0.001,
                               sample.id = keep,
                               snp.id=snpset.id) # keep all

eig = pca.out$eigenval[!is.na(pca.out$eigenval)]
barplot(100*eig/sum(eig), main="PCA Eigenvalues")

# scale 
eig <- 100*eig/sum(eig)

dat <- as.data.frame(pca.out$eigenvect)

colnames(dat) <- c("PC1", "PC2", "PC3",colnames(dat)[4:ncol(dat)] )
dat$IDs <- pca.out$sample.id
known <- merge(dat, matched_pops, by.x="IDs", by.y="indiv")
unknown_df <- merge(dat, results, by.x="IDs", by.y="Ind.ID")
nrow(dat)
nrow(unknown_df)
library(scatterpie)

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
  xlab(paste0("PC1: ",round(eig[1], 2),"% variance")) +
  ylab(paste0("PC2: ",round(eig[3], 2),"% variance")) +
  theme_bw() +
  ggtitle('Subset Orange: PC1, PC2') +
  scale_color_manual(values=c("darkgreen","lawngreen","orange3", "red3"))+
  scale_fill_manual(values=c("orange3","red3","darkgreen", "lawngreen"))
d

ggsave("figures/PCA_subset_1_2.png",
       d, w=7, h=5)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# loop over and make sure assignments are robust

# Get indices of populations to keep 
keep_idx <- which(known_gl@pop %in% c("DarkGreen_17", "Red_21", "LightGreen_14"))

# Orange pop
orange_idx <- which(known_gl@pop == "Orange_84")

unknown_dat <- read.Structure("./data/unknown.str", 
                              ploidy=2)

all_results <- data.frame(
  iteration = integer(),
  pred.pop = character(),
  ind_id = character(),
  darkgreen = numeric(),
  red = numeric(),
  orange = numeric(),
  lightgreen = numeric(),
  assignment_prob = numeric()
)


nreps <- 100
for(i in 1:nreps){
  
  # subsample orange
  orange_sample <- sample(orange_idx, 14)
  
  # Combine and subset object
  final_idx <- sort(c(keep_idx, orange_sample))
  new_gl <- known_gl[final_idx,]
  
  gl2structure(new_gl, outfile = "ATLOnly_sub.str",
               addcolumns=new_gl@pop,
               outpath = "data")
  known_dat <- read.Structure("./data/ATLOnly_sub.str", 
                              ploidy=2)
  
  assign.X( x1=known_dat, x2=unknown_dat, 
            dir="assignPOP_results_sub/", model="svm",
            mplot=F)
  results <- read.table("assignPOP_results_sub/AssignmentResult.txt",
                        header=T, sep=" ")

  results <- unique(results)  
  # get assignment probability
  results <- results %>%
    mutate(assignment_prob = case_when(
      pred.pop == "Orange_84" ~ Orange_84,
      pred.pop == "Red_21" ~ Red_21,
      pred.pop == "DarkGreen_17" ~ DarkGreen_17,
      pred.pop == "LightGreen_14" ~ LightGreen_14
    ))
  
  
  iteration_results <- data.frame(
    iteration = i,
    ind_id = results$Ind.ID,
    pred.pop = results$pred.pop,
    darkgreen = results$DarkGreen_17,
    red = results$Red_21, 
    orange = results$Orange_84,
    lightgreen = results$LightGreen_14,
    assignment_prob = results$assignment_prob
  )
  all_results <- rbind(all_results, iteration_results)

  }

# save all_results

write.table(all_results, file="analysis/structure/assignPOP_replicates.txt", 
            quote = F, row.names=F, sep="\t")





