# Clears global environment
rm(list = ls(all.names = TRUE)) 

# If needed, install and load packages
library(readxl) # for data import 
library(table1) # to make tables
library(tidyverse) # for data cleaning and organization
select <- dplyr::select
library(rstatix) # for Dunn's test
library(FSA) # package that contains Dunn test function
library(car) # for ANCOVA
library(multcomp) #for Tukey's post hoc
library(ggplot2) # for graphing
library(viridis) # for graphing

##############################################################
####### 1. Demographic Table
##############################################################

# Import data table
DemographicData <- read.csv("Input Data/2021_06_15 IS Project Demographics.csv")

# Identify factors that will be the columns in the table and specify the order they will appear.
DemographicData$Device <- factor(DemographicData$Device, 
                                 levels=c("NS/NV", "SM", "3rd Gen", "4th Gen"), 
                                 labels=c("NS/NV", "Smoker", "3rd Gen", "4th Gen"))

DemographicData$Sex <- factor(DemographicData$Sex, 
                              levels=c("M", "F"), 
                              labels=c("Male", "Female"))

DemographicData$Race <- factor(DemographicData$Race, 
                               levels=c("W", "B", "API", "MO"), 
                               labels=c("White", "Black", "Asian or Pacific Islander", "Mixed/Other"))

DemographicData$Hispanic. <- factor(DemographicData$Hispanic.,
                                    levels=c("NO", "YES"),
                                    labels=c("No", "Yes"))

# Create cleaner labels for column titles.
label(DemographicData$Hispanic.) <- "Hispanic"
label(DemographicData$Serum.Cotinine..ng.mL..ELISA) <- "Serum Cotinine (ng/mL)"

# Determine whether continuous variabes are normally or non-normally distributed
shapiro.test(DemographicData$Age) # p-value = 1.114e-07 indicates non-normally distributed
shapiro.test(DemographicData$BMI) # p-value = 3.49e-05 indicates non-normally distributed
shapiro.test(DemographicData$Serum.Cotinine..ng.mL..ELISA) # p-value = 5.627e-06 indicates non-normally distributed

# Create function for custom table so that Mean (SD) is shown for continuous variables
my.render.cont.mean.sd <- function(x, name, ...) {
  if (!is.numeric(x)) return(render.categorical.default(x))
  what <- switch(name,
                 Age = "Mean (SD)",
                 BMI  = "Mean (SD)",
                 Serum.Cotinine..ng.mL..ELISA = "Mean (SD)")
  parse.abbrev.render.code(c("", what))(x)
}


# Function for adding p-values to table
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # For numeric variables, perform a non-parametric multiple comparisons test (Kruskal Wallis)
    p <- kruskal.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a Fisher's exact test (or Chi squared if there are expected frequencies >5)
    # Simulate p value = TRUE added because groups are too small
    p <- fisher.test(table(y, g), simulate.p.value = TRUE)$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
}

# Make table. Table was exported using cmd+a on the table in Rstudio viewer -> cmd+c -> 
# paste Special in MS Word -> .html format.
table1(~ Sex + Race + Hispanic. + Age + BMI + Serum.Cotinine..ng.mL..ELISA | Device, 
       data = DemographicData, 
       render = my.render.cont.mean.sd,
       render.missing = NULL,
       overall = NULL,
       extra.col = list(`P-value` = pvalue))

# Between group significance testing. Added manually to html table in Word doc.

AgeSig <- dunnTest(Age ~ Device, data = DemographicData, method = "bh")
AgeSig <- data.frame(AgeSig$res)
AgeSig <- format(AgeSig, scientific = FALSE)

CotSig <- dunnTest(Serum.Cotinine..ng.mL..ELISA ~ Device, data = DemographicData, method = "bh")
CotSig <- data.frame(CotSig$res)
CotSig <- format(CotSig, scientific = FALSE)

##############################################################
####### 2. Sputum Cell Differentials
##############################################################

##############################################################
####### 2.1 Summary Table and Crude Analysis
##############################################################

# Import data table as data frame and make Subject IDs row labels.
DiffData <- data.frame(read_excel("Input Data/2021_09_03 IS Project Sputum Characteristics.xlsx"))
DiffData <- data.frame(DiffData[, -1], row.names = DiffData$Subject_ID)

# Remove rows with NA values for sputum characteristics (indicates no cell differential available)
DiffData <- na.omit(DiffData)

# Import device information for each subject
MediatorData_Devices <- data.frame(read_excel("Input Data/2021_06_29 IS Project Mediators.xlsx"))
MediatorData_Devices <- data.frame(MediatorData_Devices$Device, row.names = MediatorData_Devices$Subject_ID)
colnames(MediatorData_Devices) <- "Device"

# Merge data frames
DiffData <- transform(merge(MediatorData_Devices, DiffData, by = 0), row.names=Row.names, Row.names = NULL)

# Global average for percent squamous (to report if needed)
squamous <- stats.default(DiffData$PercSquam)
mean.squamous <- squamous$MEAN
sem.squamous <- squamous$SD/sqrt(squamous$N)

DiffData_Squam_Counts <- DiffData %>% group_by(Device) %>% 
  summarize(n_gt80 = sum(PercSquam > 80))

# Remove total cell count, select weight as they are dependent on the person processing the sample
# Remove % squamous as it is a marker of sample quality/contamination and not of biological importance
DiffData <- subset(DiffData, select = -c(SampleWeight, TCC, PercSquam))

# Identify factors that will be the columns in your table and specify the order you want them to appear.
DiffData$Device <- factor(DiffData$Device, 
                          levels=c("NS/NV", "SM", "3rd Gen", "4th Gen"), 
                          labels=c("NS/NV", "Smoker", "3rd Gen", "4th Gen"))


# Create row labels for variables of interest
label(DiffData$CellsPerMG) <- "Cells/mg"
label(DiffData$PercPMN) <- "% PMN"
label(DiffData$PMNPerMG) <- "PMN/mg"
label(DiffData$PercMac) <- "% Macrophage"
label(DiffData$MacPerMG) <- "Macrophages/mg"
label(DiffData$PercEos) <- "% Eosinophil"
label(DiffData$EosPerMG) <- "Eosinophils/mg"
label(DiffData$PercLym) <- "% Lymphocyte"
label(DiffData$LymPerMG) <- "Lymphocytes/mg"
label(DiffData$BronchPerMG) <- "Bronchial Cells/mg"
label(DiffData$PercBronch) <- "% Bronchial Cells"

# Test normality of each sputum metric
DiffNorm <- apply(DiffData[2:12], 2, shapiro.test)

# Create data frame to summarize results
DiffNorm <- do.call(rbind.data.frame, DiffNorm)
DiffNorm <- format(DiffNorm, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
DiffNorm$p.value.adj <- p.adjust(DiffNorm$p.value, "BH")

# Add column for normality conclusion
DiffNorm <- DiffNorm %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T))

# Create a function for the summary stats you want on your table.
my.render.cont.mean.sem <- function(x) {
  s <- stats.default(x)
  s$SEM <- with(s, SD/sqrt(N))
  with(stats.apply.rounding(s), c("",
                                  "Mean (SEM)"=sprintf("%s (%s)", MEAN, SEM)))
}

# Make table. This code requires the my.render.cont.mean.sem function defined above.
table1(~ CellsPerMG + MacPerMG + PercMac + PMNPerMG + PercPMN + EosPerMG + PercEos + LymPerMG + PercLym + BronchPerMG + PercBronch | Device, 
       data = DiffData,
       render.continuous = my.render.cont.mean.sem,
       render.missing = NULL,
       extra.col = list(`P-value` = pvalue),
       overall = NULL)

# Dunn's post hoc to determine significant between-group comparisons. 
# Significance markers added manually to table in Word. 

# Create variable and group names
groups <- names(DiffData)[1]
variables_diffs <- names(DiffData) [2:12]

# Perform Dunn's test
DiffDunnTestRes <- lapply(variables_diffs, function(x) {
  rstatix::dunn_test(DiffData, reformulate(groups, x),  
                     p.adjust.method = "BH")
})

# Summarize results
CellsPerMG <- data.frame(DiffDunnTestRes[[1]])
PercPMN <- data.frame(DiffDunnTestRes[[2]])
PMNPerMG <- data.frame(DiffDunnTestRes[[3]])
PercMac <- data.frame(DiffDunnTestRes[[4]])
MacPerMG <- data.frame(DiffDunnTestRes[[5]])
PercEos <- data.frame(DiffDunnTestRes[[6]])
EosPerMG <- data.frame(DiffDunnTestRes[[7]])
PercLym <- data.frame(DiffDunnTestRes[[8]])
LymPerMG <- data.frame(DiffDunnTestRes[[9]])
PercBronch <- data.frame(DiffDunnTestRes[[10]])
BronchPerMG <- data.frame(DiffDunnTestRes[[11]])

DiffDunnRes <- bind_rows(CellsPerMG, PercPMN, PMNPerMG, PercMac, MacPerMG, PercEos, EosPerMG, PercLym, LymPerMG, PercBronch, BronchPerMG)

DiffDunnResTrimmed <- filter(DiffDunnRes, p.adj.signif != "ns")

##############################################################
####### 2.2 ANCOVA
##############################################################

## DATA SETUP

# Imported metadata file. Kept only sex, race, and age as these were significantly different between the device groups.
MetaData <- data.frame(read_excel("Input Data/2021_07_15 IS Project Metadata.xlsx"))
MetaData <- data.frame(MetaData[ , 2:6], row.names = MetaData$Subject_ID)

# These covariates were significantly different between the exposure groups (by Chi squared or ANOVA)
MetaData <- subset(MetaData, select =c ("Sex", "Race", "Age")) 

# Collapsed race into white and non-white because groups not sufficiently large to divide more
MetaData[MetaData == "MO" | MetaData == "API" | MetaData == "B"] <- "NW" 

# Filter metadata to only include subjects with differential data. 
MetaDataDiffs <- MetaData[rownames(MetaData) %in% rownames(DiffData),]

# Merge metadata and differential data.
DiffDataAll <- transform(merge(MetaDataDiffs, DiffData, by = 0), row.names=Row.names, Row.names = NULL)

## NORMALITY ASSESSMENT

## RAW DATA

# Test normality of each sputum metric
DiffNorm <- apply(DiffData[2:12], 2, shapiro.test)

# Create data frame to summarize results
DiffNorm <- do.call(rbind.data.frame, DiffNorm)
DiffNorm <- format(DiffNorm, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
DiffNorm$p.value.adj <- p.adjust(DiffNorm$p.value, "BH")

# Add column for normality conclusion
DiffNorm <- DiffNorm %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) 
# Viewing this output shows that none of our groups are normally distributed

## LOG TRANSFORMED DATA 

# Add one to every data in the matrix and log2 transform to avoid negative values.
DiffDataLog <- log2(DiffData[2:12]+1)

# Test normality of each sputum metric
DiffNormLog <- apply(DiffDataLog, 2, shapiro.test)

# Create data frame to summarize results
DiffNormLog <- do.call(rbind.data.frame, DiffNormLog)
DiffNormLog <- format(DiffNormLog, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
DiffNormLog$p.value.adj <- p.adjust(DiffNormLog$p.value, "BH")

# Add column for normality conclusion
DiffNormLog <- DiffNormLog %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) 
# Viewing this output shows that some groups normally distributed now but not all of them.

# Graphing example histograms from mediators with different distributions. 
dev.off()

par(mfrow = c(2,3))

hist(DiffData$PercMac,
     main = "% Macrophage",
     xlab = "% Macrophage",
     ylab = "Frequency")

hist(DiffData$CellsPerMG,
     main = "Cells/mg",
     xlab = "Cells/mg",
     ylab = "Frequency")

hist(DiffData$PercBronch,
     main = "% Bronchial Cells",
     xlab = "% Bronchial Cells",
     ylab = "Frequency")

hist(DiffDataLog$PercMac,
     main = "% Macrophage Transformed",
     xlab = "% Macrophage",
     ylab = "Frequency")

hist(DiffDataLog$CellsPerMG,
     main = "Cells/mg Transformed",
     xlab = "Cells/mg",
     ylab = "Frequency")

hist(DiffDataLog$PercBronch,
     main = "% Bronchial Cells Transformed",
     xlab = "% Bronchial Cells",
     ylab = "Frequency")

# Graphing example quantile-quantile plots from mediators with different distributions. 
dev.off()

par(mfrow = c(2,3))

qqnorm(DiffData$PercMac,
       main = "% Macrophage")
qqline(DiffData$PercMac)

qqnorm(DiffData$CellsPerMG,
       main = "Cells/mg")
qqline(DiffData$CellsPerMG)

qqnorm(DiffData$PercBronch,
       main = "% Bronchial Cells")
qqline(DiffData$PercBronch)

qqnorm(DiffDataLog$PercMac,
       main = "% Macrophage Transformed")
qqline(DiffDataLog$PercMac)

qqnorm(DiffDataLog$CellsPerMG,
       main = "Cells/mg Transformed")
qqline(DiffDataLog$CellsPerMG)

qqnorm(DiffDataLog$PercBronch,
       main = "% Bronchial Cells Transformed")
qqline(DiffDataLog$PercBronch)

## Log2 transforming data did not make data distributions normal for most of the sputum metrics 
## as indicated by Shapiro-Wilk test, but it did improve some distributions to be closer to normal 
## by Shapiro-Wilk, histograms, and QQ plots, so log2 transformed data will be used for ANCOVA.

## ANCOVA

# Merge back cell differential data with metadata
AllDiffDataLog <- transform(merge(MetaData, DiffDataLog, by = 0), row.names=Row.names, Row.names = NULL)
AllDiffDataLog <- transform(merge(MediatorData_Devices, AllDiffDataLog, by = 0), row.names=Row.names, Row.names = NULL)

# Create empty data frame
ANCOVAres_diffs = data.frame()

# Add row names to data frame so that it will be able to add ANCOVA results
rownames <- c("(Intercept)", "Device", "Sex", "Race", "Age", "Residuals")
ANCOVAres_diffs <- data.frame(cbind(rownames))

# Assign row names 
ANCOVAres_diffs <- data.frame(ANCOVAres_diffs[, -1], row.names = ANCOVAres_diffs$rownames)

# Perform ANCOVA over all columns in 
for (i in 5:ncol(AllDiffDataLog)) {
  
  fit = aov(as.formula(paste0(names(AllDiffDataLog)[i], "~", paste0("`", names(AllDiffDataLog[1:4]), "`", collapse="+"), sep = "")), AllDiffDataLog)
  res <- data.frame(Anova(fit, type="III"))
  res <- subset(res, select = Pr..F.)
  names(res)[names(res) == 'Pr..F.'] <- noquote(paste0(names(AllDiffDataLog[i])))
  
  ANCOVAres_diffs <- transform(merge(ANCOVAres_diffs, res, by = 0), row.names=Row.names, Row.names = NULL)
}

# Put columns in alphabetical order
ANCOVAres_diffs <- ANCOVAres_diffs[, order(names(ANCOVAres_diffs))]

# Transpose for easy viewing
ANCOVAres_diffs <- data.frame(t(ANCOVAres_diffs))

# Delete columns for intercept and residuals
ANCOVAres_diffs <- subset(ANCOVAres_diffs, select = c(Device, Sex, Age, Race))

# Summarize Results
DeviceSigMetrics <- as.vector(row.names(filter(ANCOVAres_diffs, Device < 0.05)))
DeviceSigMetrics

AgeSigMetrics <- as.vector(row.names(filter(ANCOVAres_diffs, Age < 0.05)))
AgeSigMetrics

# Because there are so few significant associations with covariates, crude analysis with 
# Kruskal-Wallis followed by Dunn's multiple comparison test will be used in the paper. 

##############################################################
####### 3. Soluble Mediators
##############################################################

##############################################################
####### 3.1 Soluble Mediator Table and Crude Analysis
##############################################################

# Import data table as data frame and make Subject IDs row labels.
MediatorData <- data.frame(read_excel("Input Data/2021_06_29 IS Project Mediators.xlsx"))
MediatorData_Devices <- data.frame(MediatorData$Device, row.names = MediatorData$Subject_ID)
colnames(MediatorData_Devices) <- "Device"
MediatorData <- data.frame(MediatorData[ , 3:47], row.names = MediatorData$Subject_ID)

# Filter out mediators that were detected in fewer than 25 subjects (~25% of the cohort, since we have four groups).
MediatorData <- MediatorData[, which(as.numeric(colSums(MediatorData != 0)) > 25)]

# Impute zeroes as the sqrt of the lowest value in that column.
MediatorData <- as.data.frame(apply(MediatorData[ , 1:43], 2, function(x) replace(x, x == 0, sqrt(min(x[x>0])))))

# Put column names in alphabetical order
MediatorData <- MediatorData[,order(names(MediatorData))]

# Add groups back in to data frame
MediatorData <- transform(merge(MediatorData_Devices, MediatorData, by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)

# Identify factors that will be the columns in your table and specify the order you want them to appear.
MediatorData$Device <- factor(MediatorData$Device, 
                              levels=c("NS/NV", "SM", "3rd Gen", "4th Gen"), 
                              labels=c("NS/NV", "Smoker", "3rd Gen", "4th Gen"))

# Test normality of each sputum metric
MedNorm <- apply(MediatorData[1:43], 2, shapiro.test)

# Create data frame to summarize results
MedNorm <- do.call(rbind.data.frame, MedNorm)
MedNorm <- format(MedNorm, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
MedNorm$p.value.adj <- p.adjust(MedNorm$p.value, "BH")

# Add column for normality conclusion
MedNorm <- MedNorm %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) # none of our groups are normally distributed

# Copy and paste output of this formula into Table1 function (couldn't get this to work by inserting it into line in function)
as.formula(paste0("~", paste0("`", names(MediatorData[2:ncol(MediatorData)]), "`", collapse="+"), sep = "")) 

# Make table. Cmd+A and copy-paste table into word document. Changed some formatting to make the table more condensed. 
table1(~ Albumin + bFGF + CRP + dsDNA + Eotaxin + Eotaxin3 + Flt1 + GMCSF + 
         IFNg + IL10 + IL12p40 + IL12p70 + IL13 + IL15 + IL16 + IL17 + 
         IL1a + IL1B + IL2 + IL4 + IL5 + IL6 + IL7 + IL8 + IP10 + MCP1 + 
         MIP1a + MIP1B + MMP2 + MMP9 + MPO + NE + PIGF + SAA + sICAM1 + 
         sVCAM1 + TARC + Tie2 + TNFa + Uteroglobin + VEGF + VEGFC + 
         VEGFD | Device, 
       data = MediatorData,
       render.continuous = my.render.cont.mean.sem,
       render.missing = NULL,
       extra.col = list(`P-value` = pvalue),
       overall = NULL)

## DUNN'S NON-PARAMETRIC MULTIPLE COMPARISONS TEST

# Create group variable and variables of interest to be tested
groups <- names(MediatorData)[1]
variables_mediators <- names(MediatorData) [2:44]

# Apply Dunn's Test
MediatorDunnTestRes <- lapply(variables_mediators, function(x) {
  rstatix::dunn_test(MediatorData, reformulate(groups, x),  
                     p.adjust.method = "BH")
})

# Summarize Dunn Test results in one data frame
MediatorDunnRes <- data.frame()

for (i in 1:43) {
  res_df <- data.frame(MediatorDunnTestRes[[i]])
  MediatorDunnRes <- bind_rows(MediatorDunnRes, res_df)
}

# Pull out only comparisons that are significant. Significance symbols added manually in word document. 
MediatorDunnResTrimmed <- filter(MediatorDunnRes, p.adj.signif != "ns")

##############################################################
####### 3.2 ANCOVA
##############################################################

## DATA SETUP

# Mediator data and metadata already imported, but need to merge metadata and device data. 
MetaData <- transform(merge(MediatorData_Devices, MetaData, by=0 ,all=TRUE), row.names=Row.names, Row.names=NULL)

# Removing device column from mediator data.
MediatorData <- MediatorData[2:44]

## NORMALITY ASSESSMENT

# RAW DATA

# Test normality of each sputum metric
MedNorm <- apply(MediatorData, 2, shapiro.test)

# Create data frame to summarize results
MedNorm <- do.call(rbind.data.frame, MedNorm)
MedNorm <- format(MedNorm, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
MedNorm$p.value.adj <- p.adjust(MedNorm$p.value, "BH")

# Add column for normality conclusion. This produces a data frame that summarizes the results of the normality test.
MedNorm <- MedNorm %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) # none of our groups are normally distributed

# LOG TRANSFORM DATA 

# Add one to every data in the matrix and log2 transform to avoid negative values. 
# A pseudocount of 1 is added to all values in the data frame because log(0) gives an error. 
# This also ensures all resulting values are positive.
MediatorDataLog <- log2(MediatorData+1)

# Test normality of each sputum metric
MedNormLog <- apply(MediatorDataLog, 2, shapiro.test)

# Create data frame to summarize results
MedNormLog <- do.call(rbind.data.frame, MedNormLog)
MedNormLog <- format(MedNormLog, scientific = FALSE)

# Add column to adjust for multiple hypothesis testing
MedNormLog$p.value.adj <- p.adjust(MedNormLog$p.value, "BH")

# Add column for normality conclusion. This produces a data frame that summarizes the results of the normality test.
MedNormLog <- MedNormLog %>% mutate(normal = ifelse(p.value.adj < 0.05, F, T)) # Some groups normally distributed now but not all of them.

# Graphing example histograms from mediators with different distributions. 
dev.off()

par(mfrow = c(2,3))

hist(MediatorData$IL12p40,
     main = "IL12p40",
     xlab = "Concentration (pg/mL)",
     ylab = "Frequency")

hist(MediatorData$Uteroglobin,
     main = "Uteroglobin",
     xlab = "Concentration (pg/mL)",
     ylab = "Frequency")

hist(MediatorData$MMP2,
     main = "MMP2",
     xlab = "Concentration (pg/mL)",
     ylab = "Frequency")

hist(MediatorDataLog$IL12p40,
     main = "IL12p40",
     xlab = "Psudeo Log2 Concentration (pg/mL)",
     ylab = "Frequency")

hist(MediatorDataLog$Uteroglobin,
     main = "Uteroglobin",
     xlab = "Concentration (pg/mL)",
     ylab = "Frequency")

hist(MediatorDataLog$MMP2,
     main = "MMP2",
     xlab = "Psudeo Log2 Concentration (pg/mL)",
     ylab = "Frequency")

# Graphing example QQ plots from mediators with different distributions. 
dev.off()

par(mfrow = c(2,3))

qqnorm(MediatorData$IL12p40,
       main = "IL12p40")
qqline(MediatorData$IL12p40)

qqnorm(MediatorData$Uteroglobin,
       main = "Uteroglobin")
qqline(MediatorData$Uteroglobin)

qqnorm(MediatorData$MMP2,
       main = "MMP2")
qqline(MediatorData$MMP2)

qqnorm(MediatorDataLog$IL12p40,
       main = "IL12p40 Transformed")
qqline(MediatorDataLog$IL12p40)

qqnorm(MediatorDataLog$Uteroglobin,
       main = "Uteroglobin Transformed")
qqline(MediatorDataLog$Uteroglobin)

qqnorm(MediatorDataLog$MMP2,
       main = "MMP2 Transformed")
qqline(MediatorDataLog$MMP2)

## Log2 transforming data did not make data distributions normal for most of the mediators as indicated by Shapiro-Wilk test, 
## but it did improve many distributions to be closer to normal by Shapiro-Wilk, histograms, and QQ plots, so log2 transformed data will be used for ANCOVA

## ANCOVA

# Merge data
AllData <- transform(merge(MetaData, MediatorData, by = 0), row.names=Row.names, Row.names = NULL)
AllDataLog <- transform(merge(MetaData, MediatorDataLog, by = 0), row.names=Row.names, Row.names = NULL)

# Create empty data frame
ANCOVAres_mediators = data.frame()

# Add row names to data frame so that it will be able to be merged with ANCOVA results
rownames <- c("(Intercept)", "Device", "Sex", "Race", "Age", "Residuals")
ANCOVAres_mediators <- data.frame(cbind(rownames))

# Assign row names and revert back to empty data frame
ANCOVAres_mediators <- data.frame(ANCOVAres_mediators[, -1], row.names = ANCOVAres_mediators$rownames)

# Perform ANCOVA over all columns in your data frame. 
# AllDataLog[1:4] specifies the metadata columns in the input and may need to be changed depending on how many metadata columsn you have. 
# I specified 5:ncol because my mediator data started in column 5.
for (i in 5:ncol(AllDataLog)) {
  
  fit = aov(as.formula(paste0(names(AllDataLog)[i], "~", paste0("`", names(AllDataLog[1:4]), "`", collapse="+"), sep = "")), AllDataLog)
  res <- data.frame(Anova(fit, type="III"))
  res <- subset(res, select = Pr..F.)
  names(res)[names(res) == 'Pr..F.'] <- noquote(paste0(names(AllDataLog[i])))
  
  ANCOVAres_mediators <- transform(merge(ANCOVAres_mediators, res, by = 0), row.names=Row.names, Row.names = NULL)
}

# Put columns in alphabetical order
ANCOVAres_mediators <- ANCOVAres_mediators[, order(names(ANCOVAres_mediators))]

# Transpose for easy viewing
ANCOVAres_mediators <- data.frame(t(ANCOVAres_mediators))

# Delete columns for intercept and residuals so there is a summary table of just data of interest.
# This table is reported in the paper.
ANCOVAres_mediators <- subset(ANCOVAres_mediators, select = c(Device, Sex, Age, Race))
ANCOVAres_mediators <- round(ANCOVAres_mediators, 4)

# Summarize Results
DeviceSigMediators <- as.vector(row.names(filter(ANCOVAres_mediators, Device < 0.05)))
RaceSigMediators <- as.vector(row.names(filter(ANCOVAres_mediators, Race < 0.05)))
AgeSigMediators <- as.vector(row.names(filter(ANCOVAres_mediators, Age < 0.05)))
SexSigMediators <- as.vector(row.names(filter(ANCOVAres_mediators, Sex < 0.05)))

## ANCOVA POST-HOC TESTING 

# SETUP

# Filter data to only those significant for Device in ANCOVA
DeviceMediatorsLog <- select(AllDataLog, all_of(DeviceSigMediators))

# Merge demographic data back in
DeviceMediatorsLog <- transform(merge(MetaData, DeviceMediatorsLog, by = 0), row.names=Row.names, Row.names = NULL)

# Set devices as factors
DeviceMediatorsLog$Device <- factor(DeviceMediatorsLog$Device, levels = c("NS/NV", "SM", "3rd Gen", "4th Gen"))

# DUNNETT'S WITH CONTROL

# Create empty data frame
DunnettRes <- data.frame()

# Add row names to data frame so that it will be able to add results
rownames_dunnett <- c("SM - NS/NV", "3rd Gen - NS/NV", "4th Gen - NS/NV")
DunnettRes <- data.frame(cbind(rownames_dunnett))

# Assign row names
DunnettRes <- data.frame(DunnettRes[, -1], row.names = DunnettRes$rownames_dunnett)

# Perform Dunnett's Test
for (i in 5:ncol(DeviceMediatorsLog)) {
  
  fit = aov(as.formula(paste0(names(DeviceMediatorsLog)[i], "~", paste0("`", names(DeviceMediatorsLog[1:4]), "`", collapse="+"), sep = "")), 
            DeviceMediatorsLog)
  
  posthoc <- summary(glht(fit, linfct = mcp(Device = "Dunnett")), test = adjusted("none"))
  
  res <- summary(posthoc)$test
  
  res_df <- data.frame(cbind (res$coefficients, res$sigma, res$tstat, res$pvalues))
  colnames(res_df) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)")
  res_df <- res_df[4]
  names(res_df)[names(res_df) == 'Pr(>|t|)'] <- noquote(paste0(names(DeviceMediatorsLog[i])))
  
  DunnettRes <- transform(merge(DunnettRes, res_df, by = 0), row.names=Row.names, Row.names = NULL)
}

# Optional code for applying p-value adjustment (I did not use this)
# DunnettResAdj = data.frame(apply(DunnettRes, 2, p.adjust, method = "BH"))

# TUKEY'S ALL COMPARISONS

# Create empty data frame
TukeyRes <- data.frame()

# Add row names to data frame so that it will be able to add ANCOVA results
rownames_tukey <- c("SM - NS/NV", "3rd Gen - NS/NV", "4th Gen - NS/NV", "3rd Gen - SM", "4th Gen - SM", "4th Gen - 3rd Gen")
TukeyRes <- data.frame(cbind(rownames_tukey))

# Assign row names 
TukeyRes <- data.frame(TukeyRes[, -1], row.names = TukeyRes$rownames_tukey)

# Perform Tukey test
for (i in 5:ncol(DeviceMediatorsLog)) {
  
  fit = aov(as.formula(paste0(names(DeviceMediatorsLog)[i], "~", paste0("`", names(DeviceMediatorsLog[1:4]), "`", collapse="+"), sep = "")), 
            DeviceMediatorsLog)
  
  posthoc <- summary(glht(fit, linfct = mcp(Device = "Tukey")), test = adjusted("none"))
  
  res <- summary(posthoc)$test
  
  res_df <- data.frame(cbind (res$coefficients, res$sigma, res$tstat, res$pvalues))
  colnames(res_df) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)")
  res_df <- res_df[4]
  names(res_df)[names(res_df) == 'Pr(>|t|)'] <- noquote(paste0(names(DeviceMediatorsLog[i])))
  
  TukeyRes <- transform(merge(TukeyRes, res_df, by = 0), row.names=Row.names, Row.names = NULL)
}

# Examined p values only for 3rd versus 4th gen users only.

# GRAPHING DEVICE SIGNIFICANT MEDIATORS

# Convert data into a form that can be used as input for faceting
DeviceMediatorsLogLong <- gather(DeviceMediatorsLog, key = "Mediator", value = "Log2Concentration", DeviceSigMediators)

# Create panel of graphs using fact_wrap
dev.off()

theme_set(theme_bw())

pdf("Output Figures/SolubleMediatorFigurePanelJitter.pdf",
    colormodel = "cmyk",
    bg = "transparent",
    width = 7.3,
    height = 8.1)

figurepanel <- ggplot(DeviceMediatorsLogLong, aes (x = Device, y = Log2Concentration, fill = Device)) +
  geom_boxplot(color = "black",
               outlier.shape = NA) +
  geom_jitter(position = position_jitter(0.15),
              shape = 20) +
  labs(y = "Log2 (Concentration (pg/mL))") +
  scale_y_continuous(expand = expansion(mult = c(0.25, 0.7))) +
  facet_wrap(~Mediator, scales = "free_y", nrow = 4) +
  theme(strip.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_viridis(begin = 0.25, end = 1, discrete = TRUE)

figurepanel

dev.off()