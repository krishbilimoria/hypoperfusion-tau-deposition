# Hypoperfusion precedes tau deposition in the entorhinal cortex: a retrospective evaluation of ADNI-2 data

Authors: Anish Kapadia, Krish Bilimoria, Prarthna Desai, James T. Grist, Fulvio Zaccagna

Data: Collected from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database: adni.loni.usc.edu, a multicentre project with approximately 50 medical centres and university sites across the United States and Canada. This study used a subset of subject level data from the ADNI-2 database selecting only those who underwent ASL perfusion studies with post-processing quantification. 

Input:
1. 'asl-stats-sheet-amyloid-6-years.xlsx' = Amyloid uptake data
2. 'asl-stats-sheet-tau-6-years.xlsx' = Tau uptake data
3. 'ptdemog.csv' = Patient demographic information

Output: 
1. Baseline: a) Anova testing of baseline age, b) chi-squared test of baseline gender 
2. Images: Boxplots of mean cerebral blood flow across control, mild cognitive impairment, and alzheimer's disease groups
3. Primary: Kruskal-Wallis statistic comparing differences in mean cerebral blood flow between groups
4. Shapiro: Testing normality of distribution of age and mean cerebral blood flow in control, mild cognitive impairment, and alzheimer's disease groups
