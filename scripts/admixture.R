# admixture analysis
library(tidyverse)
library(dplyr)
library(forcats)
library(patchwork)

### what is our most likely k?

# read in CV scores:
cvin <- read.csv("analysis/structure/cv_all.txt", sep=":", header=F)
colnames(cvin) <- c("id", "cv")
# fix the formatting to get K into numeric format
cvin$K <- substr(cvin$id, 4, 4)
cvin <- cvin[1:4,]

# plot the results
p <- ggplot(cvin,aes(x=K,y=cv)) +
  geom_point(size=2)  + geom_line(group=1)+
  theme_bw(base_size = 9)

p
ggsave("figures/cv.png", p, h=3, w=3)

# actual results:

samplelist <- read_delim("analysis/structure/LDthin_numCorrect.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:5){
  data <- read_delim(paste0("analysis/structure/LDthin_numCorrect.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}



pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))
#out <- merge(x=dat, y=pops, by.x="IDs", by.y="indiv", all.x=T)

result <- all_data %>%
  left_join(pops, by = c("sample" = "indiv"))

all_data <- result
head(all_data)

# order samples by population:
sampleorder <- all_data$sample[order(all_data$pop)]
sampleorder_unique <- make.unique(as.character(sampleorder))
all_data$IDs <- factor(all_data$sample, levels=sampleorder_unique)
all_data$sample <- all_data$IDs

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=pop), linewidth=2) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  #scale_color_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",                    
                    labels=c("1","2"))

# within each population, order indivs by q val
new_dat<- all_data[all_data$k == 2 & all_data$Q == "Q1",]
sorted_dat <- new_dat %>%
  arrange(pop, value)
sample_order <- unique(sorted_dat$sample)
new_dat2 <- sorted_dat %>%
  mutate(sample = factor(sample, levels = sample_order))

all_data$sample <- factor(all_data$sample, levels=c(new_dat2$sample))
all_data$k <- as.numeric(all_data$k)

p2 <-  
  all_data %>%
  filter(k < 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=pop),
           linewidth = 4, 
           sides="b") +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(
      size = 5, 
      margin = margin(t = -1)  # Reduce top margin to bring labels closer to axis
    ),
    legend.position = "bottom",
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_color_manual(values=c("darkgreen","lawngreen","orange3", "red3"),
                     name = "Population") +
  scale_fill_brewer(palette = "Set1", guide="none") +
  facet_wrap(~k,ncol=1)
p2


combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
combined_plot

ggsave("figures/Admixture_plot.pdf", combined_plot, width = 8, height = 7, units="in")
ggsave("figures/Admixture_plot.png", combined_plot, width = 8, height = 7, units="in")



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# subset orange indivs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

### what is our most likely k?

# read in CV scores:
cvin <- read.csv("analysis/structure/cv_subset.txt", sep=":", header=F)
colnames(cvin) <- c("id", "cv")
# fix the formatting to get K into numeric format
cvin$K <- substr(cvin$id, 4, 4)
cvin <- cvin[1:4,]

# plot the results
p <- ggplot(cvin,aes(x=K,y=cv)) +
  geom_point(size=2)  + geom_line(group=1) +
  theme_bw(base_size = 9)

p
ggsave("figures/cv_orangeSubset.png", p, h=3, w=3)

# actual results:

samplelist <- read_delim("analysis/structure/ATLOnly_orangeSubset_LDthin_numCorrect.fam",
                         col_names = c("individual", "id2", "a", "b", "c", "d"),
                         delim=" ")

# read in all date, in a loop
## first create an empty dataframe
all_data <- tibble(individual=character(),
                   k=numeric(),
                   Q=character(),
                   value=numeric())

# then add all results to this
for (k in 2:4){
  data <- read_delim(paste0("analysis/structure/ATLOnly_orangeSubset_LDthin_numCorrect.",k,".Q"),
                     col_names = paste0("Q",seq(1:k)),
                     delim=" ")
  data$sample <- samplelist$individual
  data$k <- k
  #This step converts from wide to long.
  data %>% gather(Q, value, -sample,-k) -> data
  all_data <- rbind(all_data,data)
}



pops <- read.csv("NC_individual_ids.txt", sep="\t", col.names=c("indiv", "pop"))
#out <- merge(x=dat, y=pops, by.x="IDs", by.y="indiv", all.x=T)

result <- all_data %>%
  left_join(pops, by = c("sample" = "indiv"))

all_data <- result
head(all_data)

# order samples by population:
sampleorder <- all_data$sample[order(all_data$pop)]
sampleorder_unique <- make.unique(as.character(sampleorder))
all_data$IDs <- factor(all_data$sample, levels=sampleorder_unique)
all_data$sample <- all_data$IDs

all_data %>%
  filter(k == 2) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=pop), linewidth=2) +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw() +
  #scale_color_manual(values=c("#eac435","#557fc3", "#03cea4", "#fb4d3d"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_brewer(palette="Set1",name="K",                    
                    labels=c("1","2"))

# within each population, order indivs by q val
new_dat<- all_data[all_data$k == 2 & all_data$Q == "Q1",]
sorted_dat <- new_dat %>%
  arrange(pop, value)
sample_order <- unique(sorted_dat$sample)
new_dat2 <- sorted_dat %>%
  mutate(sample = factor(sample, levels = sample_order))


all_data$sample <- factor(all_data$sample, levels=c(new_dat2$sample))
all_data$k <- as.numeric(all_data$k)

p2 <-  
  all_data %>%
  filter(k < 5) %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_rug(aes(x=sample, y=value, color=pop),
           linewidth = 4, 
           sides="b") +
  geom_bar(stat="identity",position="stack") +
  xlab("Sample") + ylab("Ancestry") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(
      size = 5, 
      margin = margin(t = -1)  # Reduce top margin to bring labels closer to axis
    ),
    legend.position = "bottom",
    axis.title.x = element_blank()
  ) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_color_manual(values=c("darkgreen","lawngreen","orange3", "red3"),
                     name = "Population") +
  scale_fill_brewer(palette = "Set1", guide="none") +
  facet_wrap(~k,ncol=1)
p2


combined_plot <- wrap_plots(p, p2, heights = c(0.15, 1), ncol=1)
combined_plot

ggsave("figures/Admixture_orangeSubset_plot.pdf", combined_plot, width = 8, height = 7, units="in")
ggsave("figures/Admixture_orangeSubset_plot.png", combined_plot, width = 8, height = 7, units="in")





