

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

setwd("~/lozupone_lab/microbe-immune-relationships/scatterplots")
source("/Users/jenniferfouquier/lozupone_lab/rscripts/taxonomy-from-hash.R")

metadata = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/biopsy-input-renamed-for-filtering.txt", sep = "\t")
cytof.data = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/CD4CCRX5-time1.txt", sep = "\t")
greengenes_taxon <- read.csv('taxonomy-gg.tsv', sep = '\t')


LoadQIIMERepSeq = function(repseqs.transposed) {
  df = read.csv(repseqs.transposed, sep = "\t", header = FALSE,  stringsAsFactors=FALSE)
  colnames(df) = df[2,] # fix colnames
  df <- df[-c(1,2), ] # remove rows and keep all columns
  df[,2:length(colnames(df))] <- lapply(df[,2:length(colnames(df))], 
                                        function(x) (as.numeric(as.character(x) )))
  return(df)
}

transposed.feature.table = LoadQIIMERepSeq("feature-table-5232-10freq-22samples.tsv")

taxa = transposed.feature.table
taxa.colnames = colnames(taxa)
taxalist = taxa.colnames[c(2:length(taxa.colnames))]

metadata = merge(x=metadata, y=cytof.data, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=taxa, by.x='SampleID', by.y='SampleID', all = TRUE)

MetadataForLMEs = function(metadata){
  # make function for making linear mixed effects models make sense
  metadata$HIV_Status[metadata$HIV_Status == "Negative"] <- ".Negative"
  return(metadata)
}

PlotTaxaLineGraphs <- function(metadata, lme.df, response.var, 
                               facet.vars, facet.cols, greengenes_taxon){
  
  hash.list = c(lme.df$important.taxa.list)

  important.columns.for.graphs = c('PID', 'HIV_Status', response.var)
  important.columns = c(important.columns.for.graphs, hash.list)

  # only keep necessary columns; long format is bulky but needed for facet wrap
  metadata.condensed = metadata[important.columns]

  metadata.long.format = gather(metadata.condensed,
                                key = "Taxa",
                                value = "TaxaAbundance", 
                                all_of(c(hash.list)))
  
  metadata.long.format <- merge(metadata.long.format, greengenes_taxon,
                                by.x = 'hash.list', by.y = 'Feature.ID')
  
  metadata.long.format %>% 
    separate(Taxon, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family',
                             'Genus', 'Species'), 
             sep = '; ') -> metadata.long.format
  
  gg = ggplot(metadata.long.format, cor.coef = TRUE, aes(x = TaxaAbundance, y = get(response.var), 
                                   color = HIV_Status, 
                                   group = HIV_Status)) +
    ylab(response.var) +
    xlab("Abundance") +
    geom_point() +
    geom_smooth(method='lm', se = FALSE) +
    scale_color_manual(values = c("red","blue")) +
    # # Change the size of the texts
    # sp + geom_text(size=6)
    # # Change vertical and horizontal adjustement
    # sp +  geom_text(hjust=0, vjust=0)
    # # Change fontface. Allowed values : 1(normal),
    # # 2(bold), 3(italic), 4(bold.italic)
    # sp + geom_text(aes(fontface=2))
    facet_wrap(facet.vars, ncol = facet.cols, scales = "free_x")
    # facet_wrap_paginate(facets = c("Taxa"), nrow = 2, ncol = 2, page = 1)
  # Here we add our special class
  # class(p) <- c('gg_multiple', class(p))
  # 
  # # You can now print or save all pages
  # ggsave(p, filename = 'multiple_plot.pdf')
  # print(gg)
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
    # spearman rank rho/r, so only using one time point
    
    cor.results <- cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, method = "spearman", na.action = "na.omit", exact = FALSE)
    # cor.results = cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, method = "spearman", na.action = "na.exclude")
    # Linear Mixed Effects Models method (for more than one time point; likely not used
    # for this project)
    # full.model = lmer(as.formula(full.model.string), data=metadata, REML=FALSE)
    # reduced.model = lmer(as.formula(reduced.model.string), data=metadata, REML=FALSE)
    # model.summary = summary(full.model)
    # results = anova(full.model, reduced.model)
    # pvalue = results$`Pr(>Chisq)`[2]
    pvalue = cor.results$p.value
    rho = cor.results$estimate
    
    if (is.na(pvalue)){
    } else if (pvalue <= 0.05){
      print(taxon)
      res = summary(lm(metadata$CD4pos.CCR5pos ~ metadata[[taxon]]*metadata$HIV_Status ))
      # print(res)
      
      # full.model = lmer(CD4pos.CCR5pos ~ get(taxon)*HIV_Status + (1|HIV_Status), data = metadata, REML = FALSE)
      # reduced.model = lmer(CD4pos.CCR5pos ~ 1 + (1|HIV_Status), data = metadata, REML = FALSE)
      # 
      # anova.results = anova(full.model, reduced.model)
      # print(anova.results)
      # 
      summary = summary(lm(CD4pos.CCR5pos ~ get(taxon):HIV_Status, data = metadata))
      # print(summary)
      # 
      # print(pvalue)
      # print(tstatistic)
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
  # print(lme.df)
  # only look at taxa that are less than a q value
  # TODO FIX THIS VALUE ***************
  # lme.df = subset(lme.df, qlist < 0.8)
  lme.df = subset(lme.df, plist < 0.05)

  return(lme.df)
}


full.model.string = "CD4pos.CCR5pos ~ get(taxon)*HIV_Status + (1|PID)"
reduced.model.string = "CD4pos.CCR5pos ~ 1  + (1|PID)"
pdf.name = "lmes-CD4CCR5-bytaxa-20200128.pdf"
response.var = "CD4pos.CCR5pos"

metadata.for.lmes = MetadataForLMEs(metadata)
lme.df = AllTaxaLME(taxalist, metadata.for.lmes, full.model.string, 
                     reduced.model.string, pdf.name, response.var)

facet.vars = vars(Taxa)
facet.cols = 3
metadata.long.format = PlotTaxaLineGraphs(metadata, lme.df, response.var, 
                                          facet.vars, facet.cols, greengenes_taxon)


response.var = "CD4pos.CCR5pos"

gg = ggplot(metadata.long.format, aes(x = TaxaAbundance, y = get(response.var),
                                 color=HIV_Status)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~hash+Family, scales = 'free_x') +
  theme_bw(base_size = 18) +
  scale_color_brewer(palette = 'Dark2')

print(gg)

plot_grid(gg, gg, 
          rel_heights = c(1, 0.4),
          ncol = 1)
