

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

setwd("~/lozupone_lab/microbe-immune-relationships/scatterplots")
source("/Users/jenniferfouquier/lozupone_lab/rscripts/taxonomy-from-hash.R")

metadata = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/biopsy-input.txt", sep = "\t")
cytof.data = read.csv("~/lozupone_lab/microbe-immune-relationships/input-data/CD4CCRX5.txt", sep = "\t")

LoadQIIMERepSeq = function(repseqs.transposed) {
  df = read.csv(repseqs.transposed, sep = "\t", header = FALSE,  stringsAsFactors=FALSE)
  colnames(df) = df[2,] # fix colnames
  df <- df[-c(1,2), ] # remove rows and keep all columns
  colnames(df)[1] = "SampleID" 
  df[,2:length(colnames(df))] <- lapply(df[,2:length(colnames(df))], 
                                        function(x) (as.integer(as.character(x) )))
  return(df)
}

taxa.0.1.300.filtered = LoadQIIMERepSeq("rep-seqs-12922-filtered-DM-0.1-300-MOCKUP.tsv")
taxa = taxa.0.1.300.filtered
taxa.colnames = colnames(taxa)
taxalist = taxa.colnames[c(2:length(taxa.colnames))]

metadata = merge(x=metadata, y=cytof.data, by.x='SampleID', by.y='SampleID', all = TRUE)
metadata = merge(x=metadata, y=taxa, by.x='SampleID', by.y='SampleID', all = TRUE)

# keep only the rows that do not have .2 in SampleID column (time 2)
# We don't want information after diet modification
metadata = metadata[! grepl(".2", metadata[["SampleID"]]),]

MetadataForLMEs = function(metadata){
  # make function for making linear mixed effects models make sense
  metadata$HIV_Status[metadata$HIV_Status == "Negative"] <- ".Negative"
  return(metadata)
}


PlotTaxaLineGraphs <- function(metadata, lme.df, response.var, 
                               facet.vars, facet.cols){
  
  # taxa < 0.2 q-values from FDR adjustment
  important.taxa.list = c(lme.df$important.taxa.list)
  
  # only look at taxa with at least a specified number of datapoints
  # data.points = 15
  # non.zero.taxa.bool = metadata[important.taxa.list] != 0.0
  # taxa.cols.with.data = (colSums(non.zero.taxa.bool) > data.points) == TRUE
  # important.taxa.list = important.taxa.list[taxa.cols.with.data]
  # 
  # Other required columns
  # items needed for making plots (TODO some of these should be passed in)
  important.columns.for.graphs = c('PID', 'HIV_Status', response.var)
  
  important.columns = c(important.columns.for.graphs, important.taxa.list)

  # only keep necessary columns; long format is bulky but needed for facet wrap
  metadata.condensed = metadata[important.columns]
  
  important.taxa.list.modified = c()
  for (i in important.taxa.list){
    print(i)
    i.new = taxonName(i, "taxonomy.tsv")
    print(i.new)
    names(metadata.condensed)[names(metadata.condensed) == i] <- i.new
    important.taxa.list.modified = c(important.taxa.list.modified, i.new)
  }
  print(important.taxa.list.modified)
  
  metadata.long.format = gather(metadata.condensed,
                                key = "Taxa",
                                value = "TaxaAbundance", 
                                all_of(c(important.taxa.list.modified)))
  
  # expand color palettes
  color.count = length(unique(metadata$PID))
  mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(color.count)

  # ggplot with facet wrap: multiple plots faceted by vars passed to facet_wrap
  # TODO CHECK SCALES

  ggplot(metadata.long.format, cor.coef = TRUE, aes(x = TaxaAbundance, y = get(response.var), 
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
    # facet_wrap_paginate(facets = c("Diet", "Taxa"), nrow = 3, ncol = 6, page = 1)
  # Here we add our special class
  # class(p) <- c('gg_multiple', class(p))
  # 
  # # You can now print or save all pages
  # ggsave(p, filename = 'multiple_plot.pdf')
}
# lme.df

AllTaxaLME = function(taxalist, metadata, full.model.string, 
                      reduced.model.string, pdf.name, response.var){
  
  important.taxa.list = c()
  plist = c()
  lme.df = data.frame()
  
  for (taxon in taxalist){
    pvalue = FALSE
    print(taxon)
    # spearman rank rho/r, so only using one time point
    
    cor.results <- cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, na.action = "na.omit")
    # cor.results = cor.test(metadata[[taxon]], metadata$CD4pos.CCR5pos, method = "spearman", na.action = "na.exclude")
    print(cor.results)
    # Linear Mixed Effects Models method (for more than one time point; likely not used
    # for this project)
    # full.model = lmer(as.formula(full.model.string), data=metadata, REML=FALSE)
    # reduced.model = lmer(as.formula(reduced.model.string), data=metadata, REML=FALSE)
    # model.summary = summary(full.model)
    # results = anova(full.model, reduced.model)
    # pvalue = results$`Pr(>Chisq)`[2]
    pvalue = cor.results$p.value
    print(pvalue)
    
    if (is.na(pvalue)){
    } else if (pvalue <= 0.05){

      print(taxon)
      print(cor.results)
    } 
    plist = c(plist, pvalue)
    important.taxa.list = c(important.taxa.list, taxon)
  }
  
  qlist = p.adjust(plist, method = "fdr")
  # make a new dataframe that contains important taxa for this reln 
  # p and FDR-adjusted q values
  lme.df = data.frame(important.taxa.list, plist, qlist)
  print(lme.df)
  # only look at taxa that are less than a q value
  # TODO FIX THIS VALUE ***************
  lme.df = subset(lme.df, qlist < 0.2)
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
PlotTaxaLineGraphs(metadata, lme.df, response.var, facet.vars, facet.cols)
