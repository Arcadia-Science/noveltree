#!/usr/bin/env Rscript
# Tidyverse should already be included here
library(tidyverse)
library(reshape)
library(cowplot)
library(elbow)

# Read in summaries of MCL Inflation parameter performance.
fpaths <- list.files(pattern = 'summary.tsv')
fs <- list()
for(i in 1:length(fpaths)){
  fs[[i]] <- read.delim(fpaths[[i]])
}
res <- do.call(rbind, fs)

res$NumOGs_GT_4spp <- res$NumOGs_GT_4spp / res$NumOGs
res$NumOGs_All_spp <- res$NumOGs_All_spp / res$NumOGs

res <- melt(res, id.vars = 'InflationParam')
vars <- unique(res$variable)
plts <- list()

ylabs <- 
  c('Number of Orthogroups', '% Orthogroups with\n>= 4 Spp.', 
    '% Orthogroups\nwith All Spp.', 'Mean Copy # Per Spp/Per OG', 
    'Protein Domain Score', '% Genes in OGs', '% Genes in ssOGs', 
    'Mean % Species Overlap')

# Initialize empty vector to hold inflection point results. 
elbows <- c()

for(i in 1:length(vars)){
  var <- vars[i]
  
  if(i %in% c(1:9)){
    inflect <- 
      elbow(res[which(res$variable == vars[[i]]),c(1,3)])$InflationParam_selected
    elbows[i] <- inflect
    
    plts[[i]] <- 
      ggplot(data = res[which(res$variable == var),], 
             aes(x = InflationParam, y = value)) +
      theme_classic() + 
      geom_vline(xintercept = inflect) + 
      geom_point(size = 3) +
      geom_line() +
      ylab(ylabs[i])
  }else{
    plts[[i]] <- 
      ggplot(data = res[which(res$variable == var),], 
             aes(x = InflationParam, y = value)) +
      geom_point(size = 3) +
      geom_line() +
      ylab(ylabs[i]) + 
      theme_classic()
  }
}

og.summs <- 
  plot_grid(plts[[1]], plts[[2]], plts[[3]], 
            plts[[4]], plts[[5]], plts[[6]], 
            plts[[7]], plts[[8]], 
            ncol = 3, nrow = 3)

best.i <- median(elbows)

ggsave(og.summs, filename = 'inflation_summaries.pdf',
       height = 9, width = 9)

# And write out the inflation parameter we selected
sink("best-inflation-param.txt")
best.i
sink()

# And write out package versions
sink("version.txt")
packageVersion("tidyverse")
packageVersion("reshape")
packageVersion("cowplot")
packageVersion("elbow")
sink()
