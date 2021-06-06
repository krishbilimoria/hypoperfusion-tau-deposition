# Entorhinal Cortex: Final Script
# 2020-07-03

# Loading libraries
{
library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(scales)
}

# Writing folder directories (if not available)
{
dir.create("output")
dir.create("output/shapiro")
dir.create("output/baseline")
dir.create("output/primary")
dir.create("output/images")
}

# Reading input data files
{
ast6y <- read_excel("input/asl-stats-sheet-tau-6-years.xlsx") %>%
  as.data.frame()
colnames(ast6y) <- tolower(colnames(ast6y))

amyloid <- read_excel("input/asl-stats-sheet-amyloid-6-years.xlsx") %>%
  as.data.frame()

demog <- read_csv("input/ptdemog.csv") %>%
  as.data.frame()
colnames(demog) <- tolower(colnames(demog))

amyloid <- amyloid[,c(9,707:948)]
colnames(amyloid) <- tolower(colnames(amyloid))

ast6y <- right_join(ast6y, amyloid, by = c("imageuid" = "imageuid"))
}

## Subsetting demographics of our patients only
{
ids <- ast6y$rid

demog <- demog[match(ids, demog$rid), ]

ast6y <- right_join(ast6y, demog, by = c("rid" = "rid"))

ast6y$yearadmission <- substr(ast6y$userdate, start = 1, stop = 4) %>%
  as.numeric()
ast6y$age <- ast6y$yearadmission - ast6y$ptdobyy
}

## Subsetting data

### Entorhinal cortex
{
entorhinal.lh                       <- ast6y[,c(1:11,1118,1142,141,767,911)]
entorhinal.rh                       <- ast6y[,c(1:11,1118,1142,453,819,981)]

colnames(entorhinal.lh) <- c("colprot","rid", "viscode.x", "viscode2.x", "diagnosis", "examdate", "version", "loniuid", "imageuid", "rundate", "rawqc", "ptgender", "age", "meancbf", "tau.suvr", "amyloid.suvr")  
colnames(entorhinal.rh) <- c("colprot","rid", "viscode.x", "viscode2.x", "diagnosis", "examdate", "version", "loniuid", "imageuid", "rundate", "rawqc", "ptgender", "age", "meancbf", "tau.suvr", "amyloid.suvr")  

region.names <- c("entorhinal.lh", "entorhinal.rh")
all <- list(entorhinal.lh, entorhinal.rh) %>%
  setNames(region.names)
}

# Baseline Differences

## Differences in age between controls, mci, and ad in left and right hemispheres
{
an.age <- lapply(all, function (x) car::Anova(aov(age ~ diagnosis, x)))  %>%
  setNames(region.names) %>%
  dplyr::bind_rows() %>%
  as.data.frame

an.age <- an.age[1,]
write.csv(an.age, file="output/baseline/baseline-anova-age.csv")
}

## Differences in gender between controls, mci, and ad in left and right hemispheres
{
chisq.gender <- lapply(all, function (x) chisq.test(x$ptgender, x$diagnosis)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame

colnames(chisq.gender) <- c("x.statistic", "p.value")

write.csv(chisq.gender, file="output/baseline/baseline-chisq-gender.csv")
}

# Shapiro-Wilk Test: Assessing Normality

## Take the mean cbf + suvr from each of the brain areas, and do the Shapiro-Wilk test

### Mean CBF in Controls, MCI, and AD
{
shapiro.cn.meancbf <- lapply (all, function (x) subset(x, diagnosis == "CN")) %>%
  lapply (function (x) shapiro.test(x$meancbf)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.cn.meancbf$cohort <- "CN"
colnames(shapiro.cn.meancbf) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.cn.meancbf, "output/shapiro/shapiro-cn-meancbf.csv")

shapiro.mci.meancbf <- lapply (all, function (x) subset(x, diagnosis == "MCI")) %>%
  lapply (function (x) shapiro.test(x$meancbf)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.mci.meancbf$cohort <- "MCI"
colnames(shapiro.mci.meancbf) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.mci.meancbf, "output/shapiro/shapiro-mci-meancbf.csv")

shapiro.ad.meancbf <- lapply (all, function (x) subset(x, diagnosis == "AD")) %>%
  lapply (function (x) shapiro.test(x$meancbf)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.ad.meancbf$cohort <- "AD"
colnames(shapiro.ad.meancbf) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.ad.meancbf, "output/shapiro/shapiro-ad-meancbf.csv")

shapiro.meancbf <- rbind(shapiro.cn.meancbf, shapiro.mci.meancbf, shapiro.ad.meancbf)
shapiro.meancbf <- shapiro.meancbf[order(shapiro.meancbf$pvalue),]
write.csv(shapiro.meancbf, "output/shapiro/shapiro-total-meancbf.csv") #Combined rbinded values

shapiro.meancbf[(shapiro.meancbf$pvalue < 0.05),] %>%
  write.csv("output/shapiro/shapiro-not-normal-meancbf.csv")
}

### Tau in Controls, MCI, and AD
{
shapiro.cn.tau <- lapply (all, function (x) subset(x, diagnosis == "CN")) %>%
  lapply (function (x) shapiro.test(x$tau.suvr)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.cn.tau$cohort <- "CN"
colnames(shapiro.cn.tau) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.cn.tau, "output/shapiro/shapiro-cn-tau.csv")

shapiro.mci.tau <- lapply (all, function (x) subset(x, diagnosis == "MCI")) %>%
  lapply (function (x) shapiro.test(x$tau.suvr)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.mci.tau$cohort <- "MCI"
colnames(shapiro.mci.tau) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.mci.tau, "output/shapiro/shapiro-mci-tau.csv")

#shapiro.ad.tau <- lapply (all, function (x) subset(x, diagnosis == "AD")) %>%
#  lapply (function (x) shapiro.test(x$tau.suvr)) %>%
#  dplyr::bind_rows() %>%
#  t %>%
#  as.data.frame
#shapiro.mci.tau$cohort <- "AD"
#colnames(shapiro.ad.tau) <- c("wstatistic", "pvalue", "cohort") #Insufficient number of cases to run test
#write.csv(shapiro.ad.suvr, "output/shapiro.ad.suvr.csv")

shapiro.tau <- rbind(shapiro.cn.tau, shapiro.mci.tau) #,shapiro.ad.meancbf) 
shapiro.tau <- shapiro.tau[order(shapiro.tau$pvalue),]
write.csv(shapiro.tau, "output/shapiro/shapiro-total-tau.csv") #Combined

shapiro.tau[(shapiro.tau$pvalue < 0.05),] %>%
  write.csv("output/shapiro/shapiro-not-normal-tau.csv")
}

### Amyloid in Controls, MCI, AD
{
shapiro.cn.amyloid <- lapply (all, function (x) subset(x, diagnosis == "CN")) %>%
  lapply (function (x) shapiro.test(x$amyloid.suvr)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.cn.amyloid$cohort <- "CN"
colnames(shapiro.cn.amyloid) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.cn.amyloid, "output/shapiro/shapiro-cn-amyloid.csv")

shapiro.mci.amyloid <- lapply (all, function (x) subset(x, diagnosis == "MCI")) %>%
  lapply (function (x) shapiro.test(x$amyloid.suvr)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.mci.amyloid$cohort <- "MCI"
colnames(shapiro.mci.amyloid) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.mci.amyloid, "output/shapiro/shapiro-mci-amyloid.csv")

#shapiro.ad.amyloid <- lapply (all, function (x) subset(x, diagnosis == "AD")) %>%
#  lapply (function (x) shapiro.test(x$amyloid.suvr)) %>%
#  lapply (function (x) c(x$statistic, x$p.value)) %>%
#  dplyr::bind_rows() %>%
#  t %>%
#  as.data.frame
#shapiro.ad.amyloid$cohort <- "AD"
#colnames(shapiro.ad.amyloid) <- c("wstatistic", "pvalue", "cohort")
#write.csv(shapiro.ad.amyloid, "output/shapiro/shapiro-ad-amyloid.csv") #Insufficient sample size

shapiro.amyloid <- rbind(shapiro.cn.amyloid, shapiro.mci.amyloid) #,shapiro.ad.amyloid) 
shapiro.amyloid <- shapiro.amyloid[order(shapiro.amyloid$pvalue),]
write.csv(shapiro.amyloid, "output/shapiro/shapiro-total-amyloid.csv") #Combined

shapiro.amyloid[(shapiro.amyloid$pvalue < 0.05),] %>%
  write.csv("output/shapiro/shapiro-not-normal-amyloid.csv")
}

## Shapiro-Wilk for Age
{
shapiro.cn.age <- lapply (all, function (x) subset(x, diagnosis == "CN")) %>%
  lapply (function (x) shapiro.test(x$age)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.cn.age$cohort <- "CN"
colnames(shapiro.cn.age) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.cn.age, "output/shapiro/shapiro-cn-age.csv")

shapiro.mci.age <- lapply (all, function (x) subset(x, diagnosis == "MCI")) %>%
  lapply (function (x) shapiro.test(x$age)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.mci.age$cohort <- "MCI"
colnames(shapiro.mci.age) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.mci.age, "output/shapiro/shapiro-mci-age.csv")

shapiro.ad.age <- lapply (all, function (x) subset(x, diagnosis == "AD")) %>%
  lapply (function (x) shapiro.test(x$age)) %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
shapiro.ad.age$cohort <- "AD"
colnames(shapiro.ad.age) <- c("wstatistic", "pvalue", "cohort")
write.csv(shapiro.ad.age, "output/shapiro/shapiro-ad-age.csv")

shapiro.age <- rbind(shapiro.cn.age, shapiro.mci.age, shapiro.ad.age) #,shapiro.ad.meancbf) 
shapiro.age <- shapiro.age[order(shapiro.age$pvalue),]
write.csv(shapiro.age, "output/shapiro/shapiro-total-age.csv") #Combined

shapiro.age[(shapiro.age$pvalue < 0.05),] %>%
  write.csv("output/shapiro/shapiro-not-normal-age.csv")
}

# Primary Analysis

## Kruskal-Wallis test between diagnosis groups
{
kw <- lapply(all, function (x) kruskal.test(meancbf ~ diagnosis, x))  %>%
  lapply (function (x) c(x$statistic, x$p.value)) %>%
  setNames(region.names) %>%
  dplyr::bind_rows() %>%
  t %>%
  as.data.frame
colnames(kw) <- c("kw.statistic", "p.value")
}

### Benjamini-Hochberg Adjustment to p-values for multiple hypothesis testing
{
kw$pval.bh.adjusted <- p.adjust(kw$p.value, "fdr")
write.csv(kw, file="output/primary/primary-analysis-kruskal-wallis.csv")
}

# Plots
{
bp.cbf.lh <-  entorhinal.lh %>%
  ggplot( aes(x=diagnosis, y=meancbf)) + 
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot() + 
  scale_x_discrete(limits=c("CN", "MCI", "AD")) +
  scale_y_continuous(breaks=seq(0,600000,100000), labels = comma) +
  ylim(0,600000) +
  theme_classic() +
  ggtitle(paste0("Entorhinal.lh, mean cerebral blood flow"))

bp.cbf.rh <-  entorhinal.rh %>%
  ggplot( aes(x=diagnosis, y=meancbf)) + 
  stat_boxplot(geom = 'errorbar') +
  geom_boxplot() + 
  scale_x_discrete(limits=c("CN", "MCI", "AD")) +
  scale_y_continuous(breaks=seq(0,600000,100000), labels = comma) +
  ylim(0,600000) +
  theme_classic() +
  ggtitle(paste0("Entorhinal.rh, mean cerebral blood flow"))

save_plot("output/images/bp.cbf.entorhinal.lh.png", bp.cbf.lh, ncol = 1, nrow=1, base_asp = 1.1)
save_plot("output/images/bp.cbf.entorhinal.rh.png", bp.cbf.rh, ncol = 1, nrow=1, base_asp = 1.1)
}