

rm(list = ls())

library(lme4)
library(lmerTest)
library(ggplot2)
library(reshape2)
library(plotly)
library(ggpubr)
library(ggthemes)
library(RColorBrewer)
library(tidyr)
library(ggforce)
library(cowplot)
library(tibble)
library(dplyr)
library(car)
library(gridExtra)
library(grid)
library(tidyverse)
library(qiime2R)

setwd("~/lozupone_lab/microbe-immune-relationships/scatterplots")
source("/Users/jenniferfouquier/lozupone_lab/rscripts/taxonomy-from-hash.R")

metadata = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/biopsy-input-renamed-for-filtering.txt", sep = "\t")
cytof.data = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/CD4CCRX5-time1.txt", sep = "\t")
# greengenes_taxon <- read.csv("taxonomy-gg-closed-refs.tsv", sep = '\t')

greengenes_taxon <- read_qza("~/lozupone_lab/microbe-immune-relationships/microbial-analysis/merged-tables/closed-ref-2/closed-ref2-taxonomy.qza")$data

# de.novo.seqs.97 = "clustered-0.97-closed-ref-relative-freqs-transposed.tsv"
original.analysis = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/original-analysis-eiko.txt", sep = "\t")

taxa <- as.data.frame(t(read_qza("~/lozupone_lab/microbe-immune-relationships/microbial-analysis/merged-tables/closed-ref-2/relative-filtered-features.qza")$data)) %>% rownames_to_column("SampleID")

# metadata = merge(x=metadata, y=cytof.data, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=original.analysis, by.x="Biopsy_ID", by.y = "Biopsy_ID", all = TRUE)


MetadataForLMEs = function(metadata){
  # make function for making linear mixed effects models make sense
  metadata$HIV_Status[metadata$HIV_Status == "Negative"] <- ".Negative"
  return(metadata)
}

### ORIGINAL ANALYSIS 
metadata = metadata %>% filter(!is.na(Eubacterium.biforme.original))

# this is funky, remove sampleID and add it back to remove cols with all 0s
taxa.samples = metadata$SampleID
taxa %>%
  filter(SampleID %in% taxa.samples) -> taxa

taxa.samples = taxa$SampleID # store sample col after filtering
taxa <- taxa %>% select(-SampleID) # remove col
taxa <- taxa %>% select(which(!colSums(taxa, na.rm=TRUE) %in% 0)) # keep whne
taxalist = colnames(taxa)
taxa$SampleID = taxa.samples # add col back

metadata = merge(x=metadata, y=taxa, by.x='SampleID', by.y='SampleID', all = TRUE)

PlotTaxaLineGraphs <- function(metadata, lme.df, response.var){
  
  print(metadata$SampleID)
  important.taxa.list = c(lme.df$important.taxa.list)
  
  # Other required columns
  # items needed for making plots (TODO some of these should be passed in)
  important.columns.for.graphs = c('PID', 'HIV_Status', response.var)
  
  important.columns = c(important.columns.for.graphs, important.taxa.list)

  # only keep necessary columns; long format is bulky but needed for facet wrap
  metadata.condensed = metadata[important.columns]

  metadata.long.format = gather(metadata.condensed,
                                key = "hash",
                                value = "TaxaAbundance", 
                                all_of(c(important.taxa.list)))
  
  metadata.long.format <- merge(metadata.long.format, greengenes_taxon,
                                by.x = 'hash', by.y = 'Feature.ID')
  
  # make a function
  
  myRound = function(needs.rounding){
    return(round(needs.rounding, 3))
  }
  
  lme.df$text.string.stats = paste0("p=", myRound(lme.df$spear.ps), " q=", 
                              myRound(lme.df$spearman.q.list), " rho=", 
                              myRound(lme.df$rholist))
  
  # only add cols I need; not critical
  metadata.long.format <- merge(metadata.long.format, lme.df,
                                by.x = 'hash', by.y = 'important.taxa.list')
  metadata.long.format %>%
    separate(Taxon, into = c('Kingdom',
                             'Phylum',
                             'Class',
                             'Order',
                             'Family',
                             'Genus',
                             'Species'), sep = ';', remove = FALSE) -> metadata.long.format
  
  metadata.long.format$taxonomy.text = paste(metadata.long.format$Family, 
                                             metadata.long.format$Genus, 
                                             metadata.long.format$Species)
  
  metadata.long.format$taxonomy.text[metadata.long.format$hash == 125219] = "Holdemanella biformis"
  metadata.long.format$taxonomy.text[metadata.long.format$hash == 1133220] = "Mediterraneibacter faecis"
  metadata.long.format$taxonomy.text[metadata.long.format$hash == 369058] = "Desulfovibrio piger"
  metadata.long.format$taxonomy.text[metadata.long.format$hash == 469918] = "Catenibacterium mitsuokai"
  # 
  gg = ggplot(metadata.long.format, cor.coef = TRUE, aes(x = TaxaAbundance, 
                                                         y = get(response.var), 
                                                         color=HIV_Status)) +
    geom_point(aes(color = HIV_Status), size=2.7) +
    geom_smooth(method='lm', se = FALSE, color="darkgrey") +
    ylab("CD4+ CCR5+ %") +
    xlab("Relative Abundance") +
    facet_wrap(hash~taxonomy.text+text.string.stats, scales = "free_x", nrow=3) +
    theme_bw(base_size =12) +
    scale_color_manual(values=c("red", "blue"))

  print(gg)
  return(metadata.long.format)
}

AllTaxaLME = function(taxalist, metadata, full.model.string, 
                      reduced.model.string, pdf.name, response.var){
  
  important.taxa.list = c()

  spear.ps = c() 
  spearman.q.list = c()
  rholist = c()
  
  lm.p.list = c()
  lm.q.list = c()
  lme.df = data.frame()
  
  total.samples = length(metadata$SampleID)
  
  # only look at organisms found in at least a certain number of samples
  for (taxon in taxalist){
    total.zeros = sum(metadata[[taxon]]==0, na.rm = TRUE)
    if (total.zeros >= total.samples/2){
      next
    }

    # spearman rank rho/r, so only using one time point
    spearman.cor.results <- cor.test(metadata[[taxon]], 
                                     metadata$CD4_plus_CCR5_plus.original, 
                                     method = "spearman", exact = FALSE, 
                                     na.action = "na.omit", alternative = "greater")
    spear.p =  spearman.cor.results$p.value
    spearman.rho =  spearman.cor.results$estimate

    lm.results = summary(lm(CD4_plus_CCR5_plus.original ~ get(taxon)*HIV_Status, 
                            data = metadata))
    lm.interact.p = lm.results$coefficients[,4][4][[1]]

    spear.ps = c(spear.ps, spear.p)
    rholist = c(rholist, spearman.rho)
    
    lm.p.list = c(lm.p.list, lm.interact.p)
    important.taxa.list = c(important.taxa.list, taxon)
  }
  
  spearman.q.list = p.adjust(spear.ps, method = "fdr") 
  lm.q.list = p.adjust(lm.p.list, method = "fdr")
  # make a new dataframe that contains important taxa for this reln 
  # p and FDR-adjusted q values
  lme.df = data.frame(important.taxa.list, 
                      spear.ps, spearman.q.list, rholist, 
                      lm.p.list, lm.q.list)
  lme.df = subset(lme.df, spearman.q.list < 0.1)

  return(lme.df)
}

response.var = "CD4_plus_CCR5_plus.original"
lme.df = AllTaxaLME(taxalist, metadata, full.model.string, 
                     reduced.model.string, pdf.name, response.var)

metadata.long.format = PlotTaxaLineGraphs(metadata, lme.df, response.var)