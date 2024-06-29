library(tidyverse)
library(phyloseq)
library(stringr)
library(lubridate)
library(compositions)
library(padr)
library(ggrepel)
library(gridExtra)
library(abind)
library(ape)
library(ks)
library(ggridges)
library(shapes)
library(driver)
library(readr)
library(dplyr)

set.seed(123)

setwd("/home/ayden/Documents/Silverman Lab/code/fenrir_paper_code/data/mallard")
source("utils.R")

array.id <- 1

params <- read.table('sensitivity_analysis_params.txt', stringsAsFactors = FALSE)
params <- params[array.id,]

# Load Data ---------------------------------------------------------------

# Read in corrected mapping file 
mapping <- import_qiime_sample_data("2016.03.25MappingFile.MergedPool.txt")

# Read in Run Notes 
runnotes <- read_csv("2017.10.26_ResearcherIdentity.csv")
runnotes <- select(runnotes, X.SampleID, Researcher, Note, Comments) %>% 
  as.data.frame()
rownames(runnotes) <- runnotes$X.SampleID
mapping <- cbind(mapping, select(runnotes[as.character(mapping$X.SampleID), ], -X.SampleID))

# Load Sequencing batch Info and add to metadata
# t1 <- readRDS(file.path(paths[[machine]]$dada2, "seqtab.s1.nochim.rds"))
# t2 <- readRDS(file.path(paths[[machine]]$dada2, "seqtab.s2.nochim.rds"))
# t1 <- rownames(t1)
# t2 <- rownames(t2)
# batch <- ifelse(rownames(mapping) %in% t1, 1, 2)
# mapping$batch <- batch

# Read in Phyloseq Object (saved from dada2 scripts)
ps <- readRDS("phyloseq.rds")
sample_data(ps) <- mapping

# Any Increase in B. Ovatus? -----------------------------------------------

tax_table(ps)[,"Species"][!is.na(tax_table(ps)[,"Species"])] %>% 
  as.data.frame() %>% 
  arrange(Species)

# From this result is is apparent that B. Ovatus was not found by dada2 
# or not assigned a taxonomy at the species level. 

# Preprocessing -----------------------------------------------------------

# Now parse Dates from Sample IDs 
# Note: Day 0 is 2016.11.24
# Daily samples were generally taken at 3pm
sample.ids <- as.character(sample_data(ps)$X.SampleID)
postinnoc <- str_detect(sample.ids, "PostInnoc")
times <- str_split(sample.ids, "V", simplify=TRUE)[,1] %>% 
  str_replace(regex("\\.(Days|Day)"), "") %>% 
  str_replace(regex("(?!(m|d))\\.$"), "d\\.15h\\.") %>% 
  str_replace("00md\\.15h\\.", "") %>% 
  str_replace("T\\.\\.", "2015\\.11\\.") %>% 
  str_replace_all(regex("d|h"), "") %>%
  ymd_h()+days(1) 
times[is.na(times)] <- ymd_h("2015.11.1.15")
# Now set time-series to begin 
#times <- times + days(23) # Removed so that day of month = experimental day for plots
sample_data(ps)$time <- times
sample_data(ps)$postinnoc <- postinnoc



# Investigate Retained by Bioreactor --------------------------------------

ps.tmp.total <- tax_glom(ps, "Family")

first5 <- sample_data(ps.tmp.total) %>% 
  as("data.frame") %>% 
  filter(time <= min(time) + days(5)) %>% 
  pull(X.SampleID) %>% 
  as.character()

last5 <- sample_data(ps.tmp.total) %>% 
  as("data.frame") %>% 
  filter(time >= max(time) - days(5)) %>% 
  pull(X.SampleID) %>% 
  as.character()

cs.first5 <- colSums(otu_table(ps.tmp.total)[first5,])
cs.last5 <- colSums(otu_table(ps.tmp.total)[last5,])

still.there <- names(cs.first5[cs.first5 > 0]) %in% names(cs.last5[cs.last5 > 0])
print("Percentage of Families in first 5 still in last 5 days:")
(sum(still.there)/length(still.there))*100

print("These families represent X percentage of total counts")
tmp.totals <- colSums(otu_table(ps.tmp.total))
sum(tmp.totals[names(cs.first5[!still.there])])/sum(tmp.totals)*100

# Clean up some variables
rm(ps.tmp.total, first5, last5, cs.first5, cs.last5, still.there, tmp.totals)


# Continue PreProcessing --------------------------------------------------

# Investigate Distribution of Counts per sample
hist(sample_sums(ps))
quantile(sample_sums(ps), probs=seq(.01, .1, by=0.01))

# Based on this I will remove all samples with fewer than 5000 reads
total.reads <- sum(sample_sums(ps))
ps <- prune_samples(sample_sums(ps)>5000, ps)
sum(sample_sums(ps))/total.reads


# Now collapse to family level
ps <- tax_glom(ps, "Family")

# Now retain only the high abundance taxa
hist(log(taxa_sums(ps)))

ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.90*length(x)), TRUE)
(remaining.reads <- sum(sample_sums(ps))/total.reads)

# Remove Duplicate Samples (ignoring the Noise estimation samples
# and postinnoculation samples)
# Heather and I are not sure why there are duplicates... but there are
duplicates.to.remove <- subset_samples(ps, Normal_Noise_Sample=="Normal" & postinnoc==FALSE) %>%
  sample_data() %>% 
  .[,c('time', 'Vessel')] %>% 
  cbind(., duplicate=duplicated(.)) %>% 
  rownames_to_column("SampleID") %>% 
  filter(duplicate==TRUE) 
ps <- prune_samples(!(sample_names(ps) %in% duplicates.to.remove$SampleID), ps)


# Replace with Manual Created Tree ----------------------------------------

tree <- ape::read.tree("manual_families.tree")
ape::is.binary.tree(tree) # Tree must be binary for PhILR
ape::is.rooted(tree) # Tree must be rooted for PhILR
tree <- ape::makeNodeLabel(tree, method="number", prefix='n')
phy_tree(ps) <- tree


# Setup Time-Series Data --------------------------------------------------

# Keep all time-series datapoints
Y <- subset_samples(ps, Normal_Noise_Sample=="Normal" & 
                      postinnoc==FALSE) %>% 
  otu_table() %>% 
  as("matrix") %>%
  as.data.frame() %>%
  bind_cols(., as(sample_data(ps)[,c("time","Vessel")][rownames(.),], "data.frame")) %>%
  rename(t=time, series=Vessel) %>% 
  arrange(t, series)

# Index of Daily
tt.hourly <- subset_samples(ps, Normal_Noise_Sample=="Normal" &
                              SampleType=="Hourly" &
                              postinnoc==FALSE) %>%
  subset_samples(time<(ymd("2015-12-19")-days(23))) %>% 
  sample_data() %>% 
  .[["time"]] %>% 
  unique() %>% 
  as.data.frame() %>% 
  pad(interval="hour") %>% 
  .[[1]]

# Index of Hourly 
tt.daily <- Y[hour(Y$t)== 15, ] %>% 
  pad(interval="day", group="series") %>% 
  select(t) %>% 
  .[[1]] %>% 
  unique()

tt.total <- unique(Y$t) %>% 
  as.data.frame() %>% 
  pad(interval="hour") %>% 
  .[[1]]

# Pad Y to hourly
Y <- expand.grid(t = tt.total, series = factor(1:4)) %>% 
  left_join(Y, by=c("t", "series"))


# Set up Replicate Samples from End of Time-series Samples --------------------

# Create Vector of "Treatments" - this is just the vessel
Y_csme <- subset_samples(ps, Normal_Noise_Sample=="Noise_Estimation" &
                           SampleType=="Hourly" &
                           postinnoc==FALSE) %>%
  otu_table() %>%
  as("matrix")

Tr <-  subset_samples(ps, Normal_Noise_Sample=="Noise_Estimation" &
                        SampleType=="Hourly" &
                        postinnoc==FALSE) %>%
  sample_data() %>%
  .[["Vessel"]]
Tr <- as.numeric(Tr)

# Add in Time Info
replicate.datetime <- ymd_h("2015.12.22.15") - days(23)
Y_csme <- data.frame("t" = replicate.datetime, 
                    "series" = Tr, 
                     Y_csme)

# Expand for "missing" replicates
Y_csme <- Y_csme %>% 
  group_by(series) %>% 
  mutate(rep.n = 1:n()) %>% 
  ungroup() %>% 
  left_join(expand.grid("rep.n" = 1:max(.$rep.n), "series" = 1:4), ., 
             by=c("rep.n", "series"))
Y_csme$t[is.na(Y_csme$t)] <- replicate.datetime

n.replicates <- max(Y_csme$rep.n)
Y_csme <- dplyr::select(Y_csme, -rep.n) %>% 
  mutate(series = factor(series))


# Combine Longitudinal and Cross-Sectional --------------------------------
# Combining to one array indexed by sampling point rather than time point 
# (see methods). 

Y <- bind_rows(Y, Y_csme) # COMBINE!!!

tt.total <- filter(Y, series =="1")$t
tt.total.index <- seq_along(tt.total)

Y_indicators <- Y %>% 
  group_by(series) %>% 
  arrange(t) %>% 
  mutate(tt.total.index = 1:n(), 
         rep=duplicated(t)) %>%
  ungroup() %>% 
  mutate(cc = complete.cases(.)) %>% 
  select(series, tt.total.index, cc, rep)

tt.observed.ind <- Y_indicators %>% 
  select(-rep) %>% 
  spread(series, cc) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix()

tt.replicate.ind <- Y_indicators %>% 
  select(-cc) %>% 
  spread(series, rep) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix()

tt.observed <- tt.total[(rowSums(tt.observed.ind) > 0)]
tt.sample <- unique(as.POSIXct(format(c(unique(tt.observed), tt.hourly), tz="UTC", usetz = T), tz = "UTC"))

tt.sample.ind <- Y %>% 
  group_by(series) %>% 
  arrange(t) %>% 
  mutate(tt.total.index = 1:n(), 
         rep=duplicated(t)) %>%
  mutate(sample = !rep & (t %in% tt.sample)) %>% 
  ungroup() %>% 
  select(series, tt.total.index, sample) %>% 
  spread(series, sample) %>% 
  arrange(tt.total.index) %>% 
  select(-tt.total.index) %>% 
  as.matrix() %>% 
  rowSums() %>% 
  (function(x) x>0)

# Following pretty much just for plotting
tt.hourly.index <- match(tt.hourly, tt.total)
tt.daily.index <- match(tt.daily, tt.total)

# Convert to Array --------------------------------------------------------

# Now map Y to array
Y <- split(Y, Y$series) %>% 
  map(~dplyr::select(.x, -series, -t)) %>% 
  abind::abind(along=3) %>% 
  aperm(c(1,3,2)) 

Y_total <- Y
Y_total[is.na(Y_total)] <- 0
Y_total <- aperm(Y_total, c(3, 1, 2))
Y_total <- matrix(Y_total, nrow = 10, ncol = 693 * 4)
write.table(Y_total, "Y_total.csv", row.names = FALSE, col.names = FALSE, sep = ",")

# Remove time points/rows of Y that have no valid observations
Y <- Y[array_any_cases(Y),,]
Y.obs <- array_any_cases(Y, margin=c(1,2))

# Fill in NA with zeros
Y[is.na(Y)] <- 0

data <- list(Y_obs_combined = Y, observed_indices_combined = tt.observed.ind, parmas = params)

saveRDS(data, file="data.rds")
