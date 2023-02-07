#!/usr/bin/env Rscript
# Tidyverse should already be included here
library(tidyverse)
library(reshape)
library(cowplot)
library(elbow)

# Get the minimum number of species for orthogroup inclusion
# (used for plotting)
args = commandArgs(trailingOnly=TRUE)
min_spp <- args[1]

# Read in summaries of MCL Inflation parameter performance.
res <- read.delim('cogeqc_results.tsv')

res <- melt(res, id.vars = 'inflation_param')
vars <- unique(res$variable)
plts <- list()

ylabs <-
  c('Number of Orthogroups', paste0('% Orthogroups with\n>= ', min_spp, ' Species'),
    'Mean Species Copy\nNumber Per-OG', 'InterPro Score', 'OMA Score', 
    '% Genes in ssOGs', 'Mean Number of\nPer-Species ssOGs', 'Mean Pairwise\n% Species Overlap')

# In some cases we want to identify the inflation parameter that is the best
# or most representative "compromise" (e.g. % genes in ssOGs, which typically
# exhibits a pattern of hovering around some value before inflecting
# sharply). For these, we want to identify these inflection points (using elbow).
# For others, we want the value that maximizes some value (e.g InterPro score).
# Initialize empty vector to hold the results.
best <- c()
invariant <- c()

for(i in 1:length(vars)){
  # As a safety, check if the values are constant for all inflation parameters:
  # If so, we'll ignore these
  invariant[i] <- var(res$variable == vars[i]) == 0
  
  # Get the results for this summary stat  
  tmp.res <- res[which(res$variable == vars[i]),]
  # A check to make sure that we are not dealing with missing values only
  if(sum(is.na(tmp.res$value)) != length(tmp.res$value)){
    if(i %in% c(4,5)){
      # A check to make sure that we are not dealing with missing values only
      inflect <-
        elbow(tmp.res[,c(1,3)])$inflation_param_selected
      best <- c(best, inflect)
    
      plts[[i]] <-
        ggplot(data = tmp.res,
               aes(x = inflation_param, y = value)) +
        theme_classic() +
        geom_vline(xintercept = inflect) +
        geom_point(size = 3) +
        geom_line() +
        ylab(ylabs[i])
    }else{      
      plts[[i]] <-
        ggplot(data = tmp.res,
               aes(x = inflation_param, y = value)) +
        geom_point(size = 3) +
        geom_line() +
        ylab(ylabs[i]) +
        theme_classic()
    }
  }else{
      plts[[i]] <-
        ggplot(data = tmp.res,
               aes(x = inflation_param, y = value)) +
        theme_classic() +
        geom_point(size = 3) +
        geom_line() +
        ylab(ylabs[i])
  }
}

og_summs <-
  plot_grid(plts[[1]], plts[[2]], plts[[3]], plts[[4]],
            plts[[5]], plts[[6]], plts[[7]], plts[[8]],
            ncol = 4, nrow = 2)

best_i <- mean(best)

ggsave(og_summs, filename = 'inflation_summaries.pdf',
       height = 6, width = 12)

# And write out the inflation parameter we selected
sink("best_inflation_param.txt")
best_i
sink()

# And write out package versions
sink("version.txt")
packageVersion("tidyverse")
packageVersion("reshape")
packageVersion("cowplot")
packageVersion("elbow")
sink()
