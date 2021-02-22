

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

setwd("~/lozupone_lab/microbe-immune-relationships/scatterplots")
source("/Users/jenniferfouquier/lozupone_lab/rscripts/taxonomy-from-hash.R")

metadata = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/biopsy-input-renamed-for-filtering.txt", sep = "\t")
cytof.data = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/CD4CCRX5-time1.txt", sep = "\t")
greengenes_taxon <- read.csv('taxonomy-gg-de-novo-filtered-seqs.tsv', sep = '\t')

de.novo.seqs.97 = "relative-freq-de-novo-freqs-transposed.tsv"

original.analysis = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/original-analysis-eiko.txt", sep = "\t")

LoadQIIMERepSeq = function(repseqs.transposed) {
  df = read.csv(repseqs.transposed, sep = "\t", header = FALSE,  stringsAsFactors=FALSE)
  colnames(df) = df[2,] # fix colnames
  df <- df[-c(1,2), ] # remove rows and keep all columns
  # colnames(df)[1] = "SampleID" 
  df[,2:length(colnames(df))] <- lapply(df[,2:length(colnames(df))], 
                                        function(x) (as.numeric(as.character(x) )))
  return(df)
}


transposed.feature.table = LoadQIIMERepSeq(de.novo.seqs.97)

taxa = transposed.feature.table
taxa.colnames = colnames(taxa)
taxalist = taxa.colnames[c(2:length(taxa.colnames))]

metadata = merge(x=metadata, y=cytof.data, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=taxa, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=original.analysis, by.x="Biopsy_ID", by.y = "Biopsy_ID", all = TRUE)

# print(length(metadata$SampleID)) # before removing non MSM HIV Negative
# metadata = metadata[grepl("MSM", metadata[["MSM_MSW"]]) | grepl("Positive", metadata['HIV_Status']),]
# print(length(metadata$SampleID)) # after removing non MSM HIV Negative


MetadataForLMEs = function(metadata){
  # make function for making linear mixed effects models make sense
  metadata$HIV_Status[metadata$HIV_Status == "Negative"] <- ".Negative"
  return(metadata)
}



### ORIGINAL ANALYSIS
metadata = metadata %>% filter(!is.na(Eubacterium.biforme.original))

spearman.cor.results <- cor.test(metadata$`35119c7a1bef820619b6f3c4c9a9f172`, 
                                 metadata$CD4pos.CCR5pos, method = "spearman", 
                                 exact = FALSE, na.action = "na.omit")

gg = ggplot(metadata, cor.coef = TRUE, aes(x = metadata$`35119c7a1bef820619b6f3c4c9a9f172`, y = CD4pos.CCR5pos)) +
  geom_point(aes(color = HIV_Status)) +
  geom_smooth(method='lm', se = FALSE) +
  theme_bw(base_size =9)

print(gg)

# metadata = metadata %>% filter(`35119c7a1bef820619b6f3c4c9a9f172` < .1)
# metadata = metadata %>% filter(HIV_Status == "Positive" | MSM_MSW == "MSM")

PlotTaxaLineGraphs <- function(metadata, lme.df, response.var){
  
  print(metadata$SampleID)
  # taxa < 0.2 q-values from FDR adjustment
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

  metadata.long.format %>%
    separate(Taxon, into = c('Kingdom',
                             'Phylum',
                             'Class',
                             'Order',
                             'Family',
                             'Genus',
                             'Species'), sep = ';', remove = FALSE) -> metadata.long.format
  
  gg = ggplot(metadata.long.format, cor.coef = TRUE, aes(x = TaxaAbundance, 
                                                         y = get(response.var))) +
    geom_point(aes(color = HIV_Status)) +
    geom_smooth(method='lm', se = FALSE) +
    ylab(response.var)+
    xlab("Relative Abundance")+
    facet_wrap(hash~Genus+Species, scales = "free_x") +
    theme_bw(base_size =9)
    # scale_color_brewer(palette = 'Dark2')

  # plot.grid = plot_grid(gg, tableGrob(lme.df), nrow = 2)
  print(gg)
  # print(plot.grid)
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
  
  for (taxon in taxalist){
    # spearman rank rho/r, so only using one time point
    spearman.cor.results <- cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, 
                                     method = "spearman", exact = FALSE)
    spear.p =  spearman.cor.results$p.value
    spearman.rho =  spearman.cor.results$estimate
    
    if (spear.p <= 0.4 | spear.p >= 0.1){
      print(taxon)
      print(spearman.cor.results)
      print(spear.p)
    }
    
    # working
    lm.results = summary(lm(CD4pos.CCR5pos ~ get(taxon)*HIV_Status, 
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
  lme.df = data.frame(important.taxa.list, spear.ps, spearman.q.list, rholist, 
                      lm.p.list, lm.q.list)

  lme.df = subset(lme.df, spear.ps < 0.1)

  return(lme.df)
}

response.var = "CD4pos.CCR5pos"

# metadata.for.lmes = MetadataForLMEs(metadata)
lme.df = AllTaxaLME(taxalist, metadata, full.model.string, 
                     reduced.model.string, pdf.name, response.var)

write.csv(lme.df, "collapsed-taxa-pless0.2.csv")

metadata.long.format = PlotTaxaLineGraphs(metadata, lme.df, response.var)
# 
# 
# plot_grid(gg, gg, 
#           rel_heights = c(1, 0.4),
#           ncol = 1)
