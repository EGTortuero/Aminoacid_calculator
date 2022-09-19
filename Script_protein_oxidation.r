# Importing the libraries
library(pgirmess) # Multiple comparison test between treatments or treatments versus control after Kruskal-Wallis test
library(PMCMRplus) # Performs Nemenyi's non-parametric all-pairs comparison test for Kruskal-type ranked data.
library(Hmisc) # Spearman's rho rank correlation coefficients: is there any kind of trend?

# 1. Analysis of the amino acid composition for all E. coli proteins according to the subcellular location

# 1.a. Importing the whole table
setwd(".")
aa_norm <- read.csv("Table_aa_content_ecoli_wo_unknown_median_norm.csv", header=T, sep="\t", dec=".")

# 1.b. Testing if the distribution of the amino acids is homogeneous among the subcelullar location
# The reason because we run into Kruskal-Wallis test is that the data is not normally distributed.
headers <- names(aa_norm)[-c(1,2)]
results <- list()
for(i in headers){
  print(i) # Print the amino acid
  results[[i]] <- kruskal.test(formula(paste("aa_norm$", i, " ~ aa_norm$Subcellular_Location", sep="")))
}

# However, the p-values should be adjusted accordingly. We used Sequential Bonferroni to do that
p.adjust(c(results$A$p.value, results$C$p.value, results$D$p.value, results$E$p.value,
           results$F$p.value, results$G$p.value, results$H$p.value, results$I$p.value,
           results$K$p.value, results$L$p.value, results$M$p.value, results$N$p.value,
           results$P$p.value, results$Q$p.value, results$R$p.value, results$S$p.value,
           results$T$p.value, results$V$p.value, results$W$p.value, results$Y$p.value), method="BH")

# It should be nice to know, per every amino acid, which comparison(s) are significant.
# In this case, we use as a first approach the post-hoc comparison provided by pgirmess.
results2 <- list()
for(i in headers){
  print(i) # Print the amino acid
  results2[[i]] <- kruskalmc(formula(paste("aa_norm$", i, " ~ aa_norm$Subcellular_Location", sep="")))
}
print(results2) # To know the results...

# Hmmm... Not very helpful. Time to retrieve the p-value according to Nemenyi's test as
# implemented in the PMCMRplus package
results3 <- list()
for(i in headers){
  print(i) # Print the amino acid
  results3[[i]] <- kwAllPairsNemenyiTest(formula(paste(i, " ~ Subcellular_Location", sep="")), data = aa_norm)
}

# But the values are not corrected at all! Time to adjust them using Sequential Bonferroni
# However, Nemenyi test keep in mind the Family Wise correction by default (Sachs 1997)
results4 <- list()
for(i in headers){
  print(i) # Print the amino acid
  results4[[i]] <- p.adjust(kwAllPairsNemenyiTest(formula(paste(i, " ~ Subcellular_Location", sep="")), data = aa_norm)$p.value, method="BH")
}

# 1.c. Creating the tables according to the different subcellular location for further analyses
aminoacids_Cytop <- subset(aa_norm, Subcellular_Location == "1")
aminoacids_InMem <- subset(aa_norm, Subcellular_Location == "2")
aminoacids_Perip <- subset(aa_norm, Subcellular_Location == "3")
aminoacids_OuMem <- subset(aa_norm, Subcellular_Location == "4")
aminoacids_Extra <- subset(aa_norm, Subcellular_Location == "5")

# 2. Analysis per aminoacid

# 2.1. Basic statistics on the distribution of an amino acid
# The use of the median and the mad (median absolute deviation) is due to the fact that
# the distribution of the amino acids in the cell does not follow a normal distribution.
# However, it is more common to describe these values as quartiles (which could be helpful
# for the boxplots)
results_cytop <- list()
results_inmem <- list()
results_perip <- list()
results_oumem <- list()
results_extra <- list()
for(i in headers){
  print(i) # Print the amino acid
  results_cytop[[i]] <- quantile(aminoacids_Cytop[[i]], c(.025, 0.5, 0.975))
  results_inmem[[i]] <- quantile(aminoacids_InMem[[i]], c(.025, 0.5, 0.975))
  results_perip[[i]] <- quantile(aminoacids_Perip[[i]], c(.025, 0.5, 0.975))
  results_oumem[[i]] <- quantile(aminoacids_OuMem[[i]], c(.025, 0.5, 0.975))
  results_extra[[i]] <- quantile(aminoacids_Extra[[i]], c(.025, 0.5, 0.975))
}

# 2.2. Is the distribution of M more or less constant in  all subcellular compartment?
# Kruskal test to re-affirm the hypothesis that the distribution of M is uneven.
resultsKW_GENERAL <- list()
resultsKW_POSTHOC <- list()
for(i in headers){
  print(i) # Print the amino acid
  resultsKW_GENERAL[[i]] <- kruskal.test(formula(paste("aa_norm$", i, " ~ aa_norm$Subcellular_Location", sep="")))
  print(resultsKW_GENERAL[[i]])
  resultsKW_POSTHOC[[i]] <- kruskalmc(formula(paste("aa_norm$", i, " ~ aa_norm$Subcellular_Location", sep="")))
  print(resultsKW_POSTHOC[[i]])
}

# As the data had unequal variances (due to the different size), Student's t test is not 
# reliable for our analyses. As the data does not follow the normality, we cannot apply
# Welch test... So the only choice I have is the Mann-Whitney's U:
results_ci <- list()
results_cp <- list()
results_co <- list()
results_ce <- list()
results_ip <- list()
results_io <- list()
results_ie <- list()
results_po <- list()
results_pe <- list()
results_oe <- list()
results_holm <- list()
for(i in headers){
  print(i) # Print the amino acid
  results_ci[[i]] <- wilcox.test(aminoacids_Cytop[[i]], aminoacids_InMem[[i]], correct=FALSE)
  results_cp[[i]] <- wilcox.test(aminoacids_Cytop[[i]], aminoacids_Perip[[i]], correct=FALSE)
  results_co[[i]] <- wilcox.test(aminoacids_Cytop[[i]], aminoacids_OuMem[[i]], correct=FALSE)
  results_ce[[i]] <- wilcox.test(aminoacids_Cytop[[i]], aminoacids_Extra[[i]], correct=FALSE)
  results_ip[[i]] <- wilcox.test(aminoacids_InMem[[i]], aminoacids_Perip[[i]], correct=FALSE)
  results_io[[i]] <- wilcox.test(aminoacids_InMem[[i]], aminoacids_OuMem[[i]], correct=FALSE)
  results_ie[[i]] <- wilcox.test(aminoacids_InMem[[i]], aminoacids_Extra[[i]], correct=FALSE)
  results_po[[i]] <- wilcox.test(aminoacids_Perip[[i]], aminoacids_OuMem[[i]], correct=FALSE)
  results_pe[[i]] <- wilcox.test(aminoacids_Perip[[i]], aminoacids_Extra[[i]], correct=FALSE)
  results_oe[[i]] <- wilcox.test(aminoacids_OuMem[[i]], aminoacids_Extra[[i]], correct=FALSE)
  results_holm[[i]] <- p.adjust(c(results_ci[[i]][['p.value']], results_cp[[i]][['p.value']],
             results_co[[i]][['p.value']], results_ce[[i]][['p.value']],
             results_ip[[i]][['p.value']], results_io[[i]][['p.value']],
             results_ie[[i]][['p.value']], results_po[[i]][['p.value']],
             results_pe[[i]][['p.value']], results_oe[[i]][['p.value']]), method="holm")
}

# 2.3. Is there any trend in the content of an amino acid for all these proteins?
fit <- list()
for(i in headers) {
  print(i)
  fit[[i]] <- lm(formula(paste(i, " ~ Subcellular_Location", sep="")), data = aa_norm)
  print(coef(fit[[i]]))
  print(summary(fit[[i]]))
  print(cor.test(aa_norm[[i]], aa_norm[['Subcellular_Location']], method = 'spearman'))
  #print(rcorr(aa_norm[[i]], aa_norm[["Subcellular_Location"]], type = "spearman"))
}

# 3. Analysis of the amino acid composition for all E. coli inner membrane proteins according to their hydrophobicity
# 3.1. Importing the whole table
aa_norm_2 <- read.csv("Aminoacid_data_BY_PARTS_MEDIAN_NORM.csv", header=T, sep=",", dec=".")
aminoacids_Cytop <- subset(aa_norm, Part == "1")
aminoacids_Perip <- subset(aa_norm, Part == "3")
aminoacids_Trans <- subset(aa_norm, Part == "2")

# 3.2. Basic statistics on the distribution of an amino acid
# The use of the median and the mad (median absolute deviation) is due to the fact that
# the distribution of the amino acids in the cell does not follow a normal distribution.
# However, it is more common to describe these values as quartiles (which could be helpful
# for the boxplots)
results_2_cytop <- list()
results_2_trans <- list()
results_2_perip <- list()
for(i in headers){
  print(i) # Print the amino acid
  results_2_cytop[[i]] <- quantile(aminoacids_Cytop[[i]], c(.025, 0.5, 0.975))
  results_2_trans[[i]] <- quantile(aminoacids_InMem[[i]], c(.025, 0.5, 0.975))
  results_2_perip[[i]] <- quantile(aminoacids_Perip[[i]], c(.025, 0.5, 0.975))
}

# 3.3. Testing if the distribution of the amino acids is homogeneous among the subcelullar location
# The reason because we run into Kruskal-Wallis test is that the data is not normally distributed.
headers_2 <- names(aa_norm_2)[-c(1,2)]

# 3.4. Testing if the distribution of the amino acids is homogeneous among the subcelullar location
# The reason because we run into Kruskal-Wallis test is that the data is not normally distributed.
results_2_KG <- list()
results_2_KP <- list()
results_2_Nemenyi <- list()
for(i in headers_2){
  print(i) # Print the amino acid
  results_2_KG[[i]] <- kruskal.test(formula(paste("aa_norm_2$", i, " ~ aa_norm_2$Part", sep="")))
  results_2_KP[[i]] <- kruskalmc(formula(paste("aa_norm_2$", i, " ~ aa_norm_2$Part", sep="")))
  results_2_Nemenyi[[i]] <- kwAllPairsNemenyiTest(formula(paste(i, " ~ Part", sep="")), data = aa_norm_2)
}

# 3.5. Is there any trend in the content of an amino acid for all these proteins?
fit_2 <- list()
for(i in headers_2) {
  print(i)
  fit_2[[i]] <- lm(formula(paste(i, " ~ Part", sep="")), data = aa_norm_2)
  print(coef(fit_2[[i]]))
  print(summary(fit_2[[i]]))
  print(cor.test(aa_norm_2[[i]], aa_norm_2[['Part']], method = 'spearman'))
  #print(rcorr(aa_norm_2[[i]], aa_norm_2[["Part"]], type = "spearman"))
}
