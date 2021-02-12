

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
library(dplyr)


setwd("~/lozupone_lab/microbe-immune-relationships/scatterplots")
source("/Users/jenniferfouquier/lozupone_lab/rscripts/taxonomy-from-hash.R")

  
metadata = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/biopsy-input-renamed-for-filtering.txt", sep = "\t")
cytof.data = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/CD4CCRX5-time1.txt", sep = "\t")
greengenes_taxon <- read.csv('taxonomy-gg.tsv', sep = '\t')
feature.table = "rel-freqs-collapsed-species-remove-low-freq-transposed.txt"

response.var = "CD4pos.CCR5pos"

p.value.cutoff = 0.05
pdf.name = "analysis-results/regression-plots-CD4pos-CCR5pos-0.05-HIVposHIVnegMSM.pdf" 
significant.table.file = "analysis-results/collapsed-taxa-pless0.05.csv"


LoadQIIMERepSeq = function(repseqs.transposed) {
  df = read.csv(repseqs.transposed, sep = "\t", header = FALSE,  stringsAsFactors=FALSE)
  colnames(df) = df[2,] # fix colnames
  df <- df[-c(1,2), ] # remove rows and keep all columns
  df[,2:length(colnames(df))] <- lapply(df[,2:length(colnames(df))], 
                                        function(x) (as.numeric(as.character(x) )))
  return(df)
}

transposed.feature.table = LoadQIIMERepSeq(feature.table)

taxa = transposed.feature.table
taxa.colnames = colnames(taxa)
taxalist = taxa.colnames[c(2:length(taxa.colnames))]

metadata = merge(x=metadata, y=cytof.data, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=taxa, by.x='SampleID', by.y='SampleID', all = TRUE)

# HIV + and HIV- MSM for other analysis
# print(length(metadata$SampleID)) # before removing non MSM HIV Negative
# metadata = filter(metadata, MSM_MSW == "MSM" | HIV_Status == "Positive")
# print(length(metadata$SampleID)) # after removing non MSM HIV Negative

# MetadataForLMEs = function(metadata){
#   # make function for making linear mixed effects models make sense
#   metadata$HIV_Status[metadata$HIV_Status == "Negative"] <- ".Negative"
#   return(metadata)
# }

ggplotRegressions = function(metadata.long.format){
  print(response.var)
  gg = ggplot(metadata.long.format, cor.coef = TRUE, aes(x = TaxaAbundance, 
                                                         y = get(response.var), 
                                                         color = HIV_Status, 
                                                         group = HIV_Status)) +
    geom_point() +
    geom_smooth(method='lm', se = FALSE) +
    ylab(response.var)+
    xlab("Relative Abundance")+
    facet_wrap(hash~Genus+Species, scales = "free_x") +
    theme_bw(base_size =9)
  return(gg)
}


PlotTaxaLineGraphs <- function(metadata, lme.df, response.var){
  
  # taxa < 0.2 q-values from FDR adjustment
  important.taxa.list = c(lme.df$important.taxa.list)

  important.columns.for.graphs = c('PID', 'HIV_Status', response.var)
  important.columns = c(important.columns.for.graphs, important.taxa.list)

  # only keep necessary columns; long format is bulky but needed for facet wrap
  metadata.condensed = metadata[important.columns]

  metadata.long.format = gather(metadata.condensed,
                                key = "hash",
                                value = "TaxaAbundance", 
                                all_of(c(important.taxa.list)))
  
  metadata.long.format %>%
    separate(hash, into = c('Kingdom',
                             'Phylum',
                             'Class',
                             'Order',
                             'Family',
                             'Genus',
                             'Species'), sep = ';', remove = FALSE) -> metadata.long.format
  return(metadata.long.format)
}

AllTaxaLME = function(taxalist, metadata, full.model.string, 
                      reduced.model.string, pdf.name, response.var){
  
  important.taxa.list = c()
  plist = c()
  rholist = c()
  lme.df = data.frame()
  
  for (taxon in taxalist){
    pvalue = FALSE
    cor.results <- cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, method = "spearman", na.action = "na.omit", exact = FALSE)
    pvalue = cor.results$p.value
    rho = cor.results$estimate
    
    if (is.na(pvalue)){
    } else if (pvalue <= 0.05){
      print(taxon)
      res = summary(lm(metadata$CD4pos.CCR5pos ~ metadata[[taxon]]*metadata$HIV_Status ))
      summary = summary(lm(CD4pos.CCR5pos ~ get(taxon):HIV_Status, data = metadata))
      print(cor.results)
    } 
    plist = c(plist, pvalue)
    rholist = c(rholist, rho)
    important.taxa.list = c(important.taxa.list, taxon)
  }
  
  qlist = p.adjust(plist, method = "fdr")
  # make a new dataframe that contains important taxa for this reln 
  # p and FDR-adjusted q values
  lme.df = data.frame(important.taxa.list, plist, qlist, rholist)
  lme.df = subset(lme.df, plist < p.value.cutoff)

  return(lme.df)
}


lme.df = AllTaxaLME(taxalist, metadata, full.model.string, 
                     reduced.model.string, pdf.name, response.var)

write.csv(lme.df, significant.table.file)

metadata.long.format = PlotTaxaLineGraphs(metadata, lme.df, response.var)

gg = ggplotRegressions(metadata.long.format)

pdf(pdf.name)
print(gg)
dev.off()
