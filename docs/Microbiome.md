## Code used for shotgun sequencing analysis of microbiome samples of the Davar et al (2024) paper

For shotgun sequencing analysis of the microbiome samples, the [JAMS package](https://github.com/johnmcculloch/JAMS_BW) was used for both analysis *within* and *between* samples. Although the JAMS package is a *single* software package written mostly in R, the two phases of shotgun sequencing analyis, namely, analysis **within** a sample and comparison **between** samples is run independently in JAMS.

Briefly, fastqs for each sample are first put through a pipeline called [JAMSalpha](https://github.com/johnmcculloch/JAMS_BW/wiki/JAMSalpha), which yields a _single_ file with extension _.jams_ (called a jamsfile).

These jamsfiles (one for each sample) are then used, together with metadata, to obtain features-by-sample matrices through the [JAMSbeta](https://github.com/johnmcculloch/JAMS_BW/wiki/JAMSbeta) pipeline, detailed below.

### Pre-analysis of samples starting from fastqs.
In order to run fastqs for a single sample through JAMSalpha, use the JAMSalpha command on bash, with one line per sample, as in the example below:
```bash
JAMSalpha -f /path/to/forward/Sample_R1.fastq -r /path/to/reverse/Sample_R2.fastq -d /path/to/JAMSdb/JAMSdbApr2020_96Gbk2db -A metagenome -p SamplePrefix
```

The database used for this paper was built in December 2022 and can be found [here](https://hpc.nih.gov/~mccullochja/JAMSdb202212.tar.gz).

### Comparison between samples
#### Building a feature table using the outputs from JAMSalpha
The JAMSbeta pipeline will look for sample prefixes present within a column named "Sample" of the metadata file supplied and select the corresponding jamsfiles present within the folder specified following the `-y` argument:

```bash
JAMSbeta -p DavarCpG -t metadata_file.tsv -y /path/to/folder/with/jamsfiles -n LKT,Product,ECNumber,GO,Interpro -u
```

This will create an R-session image (.RData) containing [SummarizedExperiment](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html) objects which can be used with the JAMS package plotting functions.

Note the use of the "-u" option (stratify function by taxonomy) when running JAMSbeta. This will ensure that non-taxonomic (functional) SummarizedExperiment objects will also contain information as to which LKTs contribute to each functional feature count within each sample. See the code for Supplemental figure 7c-g.

## Analysis of fecal samples obtained from humans

The JAMSbeta script will have generated an object of class list called _expvec_. This named list contains as each element, a JAMS-compatible SummarizedExperiment object for all samples by all features within a given analysis space. Taxonomic data can be accessed with expvec$LKT; primary gene annotation at expvec$Product; Enzyme Commission number annotation with expvec$ECNumber; gene ontology annotation with expvec$GO and [Interpro](https://www.ebi.ac.uk/interpro/) signatures with expvec$Interpro.


We encourage you to browse through the standard [JAMS beta tutorial](https://github.com/johnmcculloch/JAMS_BW/wiki/Tutorial_JAMSbeta) in order to become familiar to plotting syntax in the JAMS package.

All code below was run in R after loading the JAMS package and code for Figures 7a, 7b, 7c, 7d, 7f, 8b, S6, S7b, S7c, and S7d was written and developed by John A. McCulloch. 


For consistency, obtain a vector with taxonomic features which survive "moderate" filtering, i.e. LKTs surviving at least 250 PPM in 15% of samples and having an estimated genome completeness of at least 10% in 5% of samples.
```R
library(JAMS)
currobj <- filter_experiment(ExpObj = expvec$LKT, featcutoff = c(250, 15),asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = c(10, 5))

moderateLKTs <- rownames(currobj)
```

Define vectors with the names of samples at specific timepoints. This will make it more convenient for plotting subsequent plots (see below).

```R
Samples_Intermediate <- c("PD1UOP004W3CpG202205", "PD1UOP005W7CpG202205", "PD1UOP007W3CpG202205", "PD1UOP008ACpG202112", "PD1UOP009CpG202112", "PD1UOP010BCpG202112", "PD1UOP011BCpG202112", "PD1UOP013BCpG202112", "PD1UOP014BCpG202112", "PD1UOP015BCpG202112", "PD1UOP016BCpG202112", "PD1UOP018BCpG202112", "PD1UOP019BCpG202112", "PD1UOP021BCpG202112", "PD1UOP022ECpG202112", "PD1UOP023BCpG202112", "PD1UOP024BCpG202112", "PD1UOP025BCpG202112", "PD1UOP026W3CpG202205", "PD1UOP027BCpG202112", "PD1UOP030BCpG202112", "PD1UOP032W3CpG202205", "PD1UOP033BCpG202112", "PD1UOP034W3CpG202205", "PD1UOP035BCpG202112", "PD1UOP036W3CpG202205", "PD1UOP001W3CpG202205", "PD1UOP002W3CpG202205", "PD1UOP037W3CpG202205")

Samples_4_pre_surgical_timepoints <- c("PD1UOP004W1CpG202205", "PD1UOP007W1CpG202205", "PD1UOP008W1CpG202205", "PD1UOP009CpG202112", "PD1UOP010ACpG202112", "PD1UOP013ACpG202112", "PD1UOP014ACpG202112", "PD1UOP015ACpG202112", "PD1UOP016ACpG202112", "PD1UOP018ACpG202112", "PD1UOP019ACpG202112", "PD1UOP021ACpG202112", "PD1UOP022W1CpG202205", "PD1UOP023ACpG202112", "PD1UOP024ACpG202112", "PD1UOP025ACpG202112", "PD1UOP026W1CpG202205", "PD1UOP027ACpG202112", "PD1UOP030ACpG202112", "PD1UOP032ACpG202112", "PD1UOP033ACpG202112", "PD1UOP034ACpG202112", "PD1UOP035ACpG202112", "PD1UOP036ACpG202112", "PD1UOP001W1CpG202205", "PD1UOP005W7CpG202205", "PD1UOP011ACpG202112", "PD1UOP002W3CpG202205", "PD1UOP037W3CpG202205", "PD1UOP008ACpG202112", "PD1UOP010BCpG202112", "PD1UOP014BCpG202112", "PD1UOP015BCpG202112", "PD1UOP016BCpG202112", "PD1UOP018BCpG202112", "PD1UOP019BCpG202112", "PD1UOP021BCpG202112", "PD1UOP022ECpG202112", "PD1UOP024BCpG202112", "PD1UOP025BCpG202112", "PD1UOP027BCpG202112", "PD1UOP033BCpG202112", "PD1UOP035BCpG202112", "PD1UOP001W3CpG202205", "PD1UOP004W3CpG202205", "PD1UOP007W3CpG202205", "PD1UOP026W3CpG202205", "PD1UOP032W3CpG202205", "PD1UOP034W3CpG202205", "PD1UOP036W3CpG202205", "PD1UOP007ACpG202112", "PD1UOP008BCpG202112", "PD1UOP010CCpG202112", "PD1UOP011BCpG202112", "PD1UOP013BCpG202112", "PD1UOP014CCpG202112", "PD1UOP015CCpG202112", "PD1UOP016CCpG202112", "PD1UOP018CCpG202112", "PD1UOP019CCpG202112", "PD1UOP021CCpG202112", "PD1UOP022FCpG202112", "PD1UOP023BCpG202112", "PD1UOP024CCpG202112", "PD1UOP025CCpG202112", "PD1UOP026CpG202112", "PD1UOP027CCpG202112", "PD1UOP030BCpG202112", "PD1UOP032BCpG202112", "PD1UOP035CCpG202112", "PD1UOP036BCpG202112", "PD1UOP004W7CpG202205", "PD1UOPO33W7CpG202205", "PD1UOP005ACpG202112", "PD1UOP007BCpG202112", "PD1UOP008CCpG202112", "PD1UOP010DCpG202112", "PD1UOP011CCpG202112", "PD1UOP015DCpG202112", "PD1UOP016DCpG202112", "PD1UOP018DCpG202112", "PD1UOP019DCpG202112", "PD1UOP022GCpG202112", "PD1UOP035DCpG202112")
```


### Code used for Figure 7

## Figure 7 panel a
```R
# Code for Figure 7a
# Use all samples
pdf("Figure_7_a.pdf", paper = "a4r")
print(plot_Ordination(ExpObj = expvec[["LKT"]], samplesToKeep = NULL, algorithm = "tUMAP",
    distmethod = "bray", compareby = "Patient_ID", colourby = NULL, shapeby = "Patient_ID", 
    connectby = "Patient_ID", connection_orderby = "Timepoint_days_labels",
    textby = "Timepoint_days_labels", dotsize = 1.3, grid = FALSE, forceaspectratio = 1, permanova = FALSE))
dev.off()

```

## Figure 7 panel b
```R
# Code for Figure 7b
pdf("Figure_7_b.pdf", paper = "a4r")
print(plot_Ordination(ExpObj = expvec[["LKT"]], samplesToKeep = Samples_Intermediate,
    featuresToKeep = moderateLKTs,  applyfilters = NULL, algorithm = "tUMAP", 
    distmethod = "bray", compareby = "MPR", colourby = "MPR",
    dotsize = 1, plotcentroids = TRUE, highlight_centroids = TRUE, grid = FALSE,
    forceaspectratio = 1, permanova = TRUE, permanova_permutations = 10000))
dev.off()

```

## Figure 7 panel c
```R
# Code for Figure 7c
pdf("Figure_7_c.pdf", paper = "a4r")

plot_relabund_heatmap(samplesToKeep = Samples_Intermediate, featuresToKeep = moderateLKTs, 
    ExpObj = expvec[["LKT"]], hmtype = "comparative", compareby = "MPR", applyfilters = NULL, 
    featcutoff = NULL, cdict = cdict, colcategories = c("MPR", "Relapse"), 
    showonlypbelow = 0.05, adj_pval_for_threshold = FALSE, cluster_rows = TRUE, 
    splitcolsby = "MPR", scaled = TRUE, minl2fc = 1, invertbinaryorder = TRUE,
    showGram = TRUE, no_underscores = TRUE, label_samples = FALSE)
```

## Figure 7 panel d
```R
# Code for Figure 7d
pdf("Figure_7_d.pdf", paper = "a4r")

print(plot_Ordination(ExpObj = expvec[["LKT"]], samplesToKeep = Samples_4_pre_surgical_timepoints, 
    featuresToKeep = moderateLKTs, applyfilters = NULL, algorithm = "tUMAP", 
    distmethod = "bray", compareby = "MPR", colourby = "MPR",
    dotsize = 1, plotcentroids = TRUE, highlight_centroids = TRUE, 
    grid = FALSE, forceaspectratio = 1, permanova = TRUE, permanova_permutations = 10000))
dev.off()

```

## Figure 7 panel f
```R
# Code for Figure 7f
pdf("Figure_7_f.pdf", paper = "a4r")

plot_relabund_features(ExpObj = expvec[["LKT"]], glomby = "Gram", 
    samplesToKeep = Samples_Intermediate, featuresToKeep = c("negative", "positive", "undetermined"),
    aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", 
    compareby = "MPR", compareby_order = c("Non-MPR", "MPR"), 
    colourby = "MPR", signiflabel = "p.format", ignoreunclassified = TRUE, y_axis_range = 1000000)

dev.off()

```

## Figure 8 panel a
Export LKT relative abundances in parts per million (PPM)
```R
export_expvec_to_XL(expvec = expvec, usefulexp = "LKT", filename = "CpG_Relabund_PPM_moderate_LKTs.xlsx", applyfilters = "moderate",  asPPM = TRUE)
```
Cladogram by Richard Rodrigues


## Figure 8 panel b
```R
# Code for Figure 8b
pdf("Figure_8_b.pdf", paper = "a4r")
for (anal in c("Product", "ECNumber", "GO", "Interpro")){
    print(plot_Ordination(ExpObj = expvec[[anal]], samplesToKeep = Samples_Intermediate,
    featuresToKeep = NULL, applyfilters = NULL, algorithm = "tUMAP", distmethod = "bray",
    compareby = "MPR", colourby = "MPR", dotsize = 1, plotcentroids = TRUE,
    highlight_centroids = TRUE, grid = FALSE, forceaspectratio = 1, permanova = TRUE,
    permanova_permutations = 10000))
}
dev.off()

```

## Supplementary Figure 6
```R
# Code for Supplementary Figure 6
pdf("Figure_S6.pdf", paper = "a4r")
for(anal in c("Product", "ECNumber", "GO", "Interpro")){
plot_relabund_heatmap(samplesToKeep = Samples_Intermediate, featuresToKeep = NULL,
    ExpObj = expvec[[anal]], hmtype = "comparative", compareby = "MPR", applyfilters = NULL,
    featcutoff = NULL, cdict = cdict, colcategories = c("MPR", "Relapse"), 
    showonlypbelow = 0.01, adj_pval_for_threshold = FALSE, cluster_rows = TRUE, 
    splitcolsby = "MPR", scaled = TRUE, minl2fc = 1, invertbinaryorder = TRUE,
    showGram = TRUE, no_underscores = TRUE, label_samples = FALSE)
}
```

## Supplementary Figure 7b

First, get a vector of the relevant LPS accessions we want to plot from the "Product" SummarizedExperiment object.
```R
Product_ftt <- as.data.frame(rowData(expvec$Product))
LPS_genes_Product_df <- Product_ftt[grep("lpx", Product_ftt$Description), ]
#Add the one whose gene symbol does not begin with "lpx"
LPS_genes_Product_df <- rbind(LPS_genes_Product_df, Product_ftt["3-deoxy-D-manno-octulosonic acid transferase", ])
LPS_genes_Product <- LPS_genes_Product_df$Accession

Relevant_LPS_accessions <- subset(LPS_genes_Product_df, Description %in% c("lpxT", "lpxP", "lpxM", "waaA", "lpxL", "lpxC", "lpxB", "lpxA", "lpxD", "lpxK", "lpxH"))[]$Accession
#Remove duplicate entry
Relevant_LPS_accessions <- Relevant_LPS_accessions[Relevant_LPS_accessions != "UDP-3-O-(3-hydroxymyristoyl)glucosamine N-acyltransferase"]
```

Plot a heatmap, scaled and non-scaled of these gene products within the "Intermediate" timepoint samples only.

```R
# Code for Supplementary Figure 7b
pdf("Figure_S7b.pdf", paper = "a4r")

for (sca in c(TRUE, FALSE)){
    plot_relabund_heatmap(samplesToKeep = Samples_Intermediate, featuresToKeep = Relevant_LPS_accessions, 
    ExpObj = expvec[["Product"]], hmtype = "comparative", compareby = "MPR", 
    applyfilters = NULL, featcutoff = NULL, cdict = NULL, colcategories = c("MPR", "Regression"), 
    showonlypbelow = NULL, adj_pval_for_threshold = FALSE, textby = "Timepoint_days_labels", 
    subsetby = NULL, max_rows_in_heatmap = 100, cluster_rows = FALSE,
    row_order = Relevant_LPS_accessions, ordercolsby = NULL, 
    splitcolsby = "MPR", column_split_group_order = c("MPR", "No_MPR"),
    scaled = sca, minl2fc = NULL, maxnumheatmaps = 2, invertbinaryorder = FALSE, 
    showGram = FALSE, cluster_samples_per_heatmap = TRUE, cluster_features_per_heatmap = TRUE,
    addtit = paste(ts, "- Only LPS related Genes"), no_underscores = TRUE, 
    label_samples = TRUE, returnstats = FALSE)
}
dev.off()

```
## Supplementary Figure 7c to g
Plot boxplots for relevant LPS signatures representing their taxonomically independent relative abundance within each sample for each group. Additionally, stratify the relative abundance of each LPS feature of interest by its contributing taxa.

```R
# Code for Supplementary Figure 7 c to g
pdf("Figure_S7c-g.pdf", paper = "a4r")

print(plot_relabund_features(ExpObj = expvec[["Product"]], glomby = NULL, 
    samplesToKeep = Samples_Intermediate, featuresToKeep = Relevant_LPS_accessions, 
    aggregatefeatures = FALSE, compareby = "MPR", compareby_order = c("MPR", "No_MPR"),
    colourby = "MPR", applyfilters = NULL, log2tran_main_plot = FALSE, log2tran_strat_plot = TRUE,
    statsonlog = FALSE, stratify_by_taxlevel = "LKT", maxnumtaxa = 20, plot_points_on_taxonomy = TRUE))
dev.off()

```

## Machine Learning section
## Supplementary Figure 8
By Jonathan Badger
```R

library(JAMS) # https://github.com/johnmcculloch/JAMS_BW
library(pROC)
library(gridExtra)
library(grid)
library(caret)
library(DESeq2)
library(edgeR)

# code based on Amiran's to make a median ROC from replicate runs stored in rds
median_roc <- function(rds, offset = 0, good="Good", addtit="") {
  data <- readRDS(rds)
  sensitivity <- vector()
  specificity <- vector()
  for (i in 1:100) {
    actual_labels <- as.character(data[[5+offset]][i][[1]])
    predicted_probabilities <- as.data.frame(data[[7+offset]][i])[[good]]
    if (is.null(predicted_probabilities)) {
      predicted_probabilities <- as.numeric(data[[7+offset]][i][[1]])
    }
    roc_curve_obj <- roc(actual_labels, predicted_probabilities)
    roc_data <- coords(roc_curve_obj)
    roc_data[order(roc_data$sensitivity), ]
    roc_data$specificity <- 1-roc_data$specificity
    sensitivity <- cbind(sensitivity, roc_data$sensitivity)
    specificity <- cbind(specificity, roc_data$specificity)
  }
  specificity.median <- apply(specificity[, -1], 1, median)
  sensitivity.median <- apply(sensitivity[, -1], 1, median)
  
  # calculate quantile
  upper_ci <- apply(sensitivity, 1, quantile, probs = 0.25)
  lower_ci <- apply(sensitivity, 1, quantile, probs = 0.75)
  z <-  data.frame(specificity.median, sensitivity.median, 
                   upper_ci, lower_ci)
  z <- z[order(z$sensitivity.median), ]
  z <- rbind(rep(0, 4), z)
  simple_auc <- function(specificity, sensitivity) {
    TPR <- sensitivity
    FPR <- 1 - specificity
    # inputs already sorted, best scores first 
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    sum(TPR * dFPR) + sum(dTPR * dFPR)/2
  }
  median_auc <- round(simple_auc(specificity.median, sensitivity.median),2)
  
  plot_line <- ggplot(z, aes(specificity.median, sensitivity.median)) +
    geom_line(color = "blue", size = 1)  +
    labs(title = str_c(addtit, "AUC=",median_auc)) +xlim(0,1) + ylim(0,1)
  
  # Add confidence intervals using geom_ribbon
  
  plot_with_intervals <- plot_line +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), 
                fill = "gray", alpha = 0.5)
  # Add a diagonal line
  plot_with_diagonal <- plot_with_intervals +
    geom_abline(intercept = 0, slope = 1, color = "red")+
    #xlim(-1, 2) +  # Set x-axis limits
    #ylim(-1, 2) +
    xlim(0, 1) +  # Set x-axis limits
    ylim(0, 1) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0,1)) + 
    theme_classic() + ylab("sensitivity") +
    xlab("1 - specificity") + theme(aspect.ratio=1) 
  
  
  plot_with_diagonal
}

# plot AUC vs # of top features
feature_plot <- function(rds) {
  data <- readRDS(rds)
  aucs <- unlist(lapply(data, function(x) 
    if(is.null(x$results$AUC)) NA else x$results$AUC))
  positions <- which(!is.na(aucs))
  aucs <- aucs[positions]
  df <- data.frame(Features=positions, AUC=aucs)
  ggplot(data=df, aes(x=Features, y=AUC)) +
    geom_line() +
    geom_point() + theme_minimal() + 
    scale_x_continuous(breaks=seq(0,200,10)) +
    theme(aspect.ratio=1) + ylim(0.5, 1)
}  

# plot rocs from top features
feature_rocs <- function(rds) {
  data <- readRDS(rds)
  aucs <- unlist(lapply(data, function(x) 
    if(is.null(x$results$AUC)) NA else x$results$AUC))
  rocs <- lapply(data, function(x) x$roc)
  positions <- which(!is.na(aucs))
  rocs <- rocs[positions]
  names(rocs) <- sapply(positions, function(x) 
    str_c(x, ": ", round(aucs[x], 2)))
  ggroc(rocs)+theme_minimal()+geom_abline(intercept = 1) +
    theme(aspect.ratio=1) 
}

# feature select using limma
feature_select <- function(exp, cutoff = 50, top=500, plot=FALSE) {
  if (top>nrow(exp)) {
    return(exp)
  }
  data <- data.frame(assay(exp))
  group <- factor(exp$Response_at_surgery_dichot)
  d0 <- DGEList(data)
  d0 <- calcNormFactors(d0)
  mm <- model.matrix(~0 + group)
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  if (length(drop)==0) {
    d <- d0[-drop,]
  } else {
    d <- d0
  }
  y <- voom(d, mm, plot = plot)
  fit <- lmFit(y, mm)
  contr <- makeContrasts(groupGood - groupBad, 
                         levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  top.table <- top.table[1:min(nrow(exp), top),]
  exp[row.names(top.table),]
}

# wrapper around caret to run ML on a SummarizedExperiment object
runML <- function(exp, response, positive, index, method="rf",
                  seed=15, filterFeature=FALSE) {
  set.seed(seed)
  train_control <- trainControl(method="cv", number=10, 
                                classProbs=T, savePredictions = T)
  results <- data.frame()
  models <- list()
  smp <- NULL
  while(is.null(smp) || length(unique(exp[,smp][[response]]))==1 ||
        length(unique(exp[,-smp][[response]]))==1) {
    smp <- sample(index, length(index), TRUE)
  }
  if (filterFeature) {
    exp <- feature_select(exp)
  }
  training <- exp[, smp]
  testing <- exp[, -smp] %>% assay() %>% t()
  model <- train(t(assay(training)), factor(training[[response]]), 
                 trControl = train_control, verbosity = 0,
                 method = method, preProcess=c("center", "scale"))
  obs <- factor(exp[,row.names(testing)][[response]])
  predictions <- factor(predict(model, testing))
  probs <- predict(model, testing, type = "prob")[,positive]
  matrix <- confusionMatrix(predictions, obs,
                            positive = positive)
  roc_obj <- roc(as.numeric(obs), probs)
  auc <- auc(roc_obj)
  stats <- matrix$overall
  stats <- append(stats, c(AUC=auc))
  results <- rbind(results, stats)
  colnames(results) <- names(stats)
  list(results=results, models=model,roc=roc_obj, matrix=matrix)
}

# run trimmed ML method on top features
runTopFeatures <- function(exp, response, positive, imp, top = 10, method="rf", trainf=0.7, seed=15, boot=101) {
  results <- list()
  aucs <- c()
  set.seed(seed)
  for(i in 1:boot) {
    printf("top:%d, boot:%i\n", top, i)
    exp <- exp[imp[1:top,1],]
    index <- sample(ncol(exp), round(ncol(exp)*trainf))
    results[[i]] <- suppressMessages(suppressWarnings(runML(exp, response, positive, index, method, seed)))
    aucs[i] <- results[[i]]$results$AUC
  }
  med <- median(aucs)
  return(results[[which(aucs==med)[1]]])
}

# run runTopFeatures over a number of feature nummbers

runMultiTopFeatures <- function(exp, response, positive, imp, top) {
  results <- list()
  for(t in top) {
    results[[t]] <- runTopFeatures(exp, response, positive, imp, top=t)
  }
  return(results)
}

# create data frame of the top scoring features across replicates
summarize_imp <- function(rds, offset = 0, top = 20, description = NULL) {
  data <- readRDS(rds)
  tot <- NULL
  for(i in 1:length(data)) {
    imp <- data[[4+offset]][[1]]
    imp <- cbind(Feature=row.names(imp),imp)
    if (is.null(tot)) {
      tot <- imp
    } else {
      shared <- intersect(row.names(imp), row.names(tot))
      di <- intersect(row.names(imp), row.names(tot))
      tot <- rbind(tot, imp[di,])
      tot[shared,2] <- tot[shared,2] + imp[shared, 2]
    }
  }
  tot <- tot[order(tot[,2], decreasing = TRUE),]
  tot <- tot[!duplicated(tot$Feature),]
  row.names(tot) <- tot$Feature
  colnames(tot)[2] <- "Score"
  tot$Score <- round(tot$Score, 1)
  if (!is.null(description)) {
    tot$Description <- description[rownames(tot),1]
  }
  head(tot, top)
}

# make a figure with roc, features, and feature selection for a given rds
plot_ml_fig <- function(rds, offset = 0, description = NULL, title = "", good="Good",
                        topFeatures=NULL) {
  tbl <- summarize_imp(rds, offset = offset, description = description)
  colnames(tbl)[3] <- "Direction"
  if (good=="Good") {
    tbl[tbl$Direction=="Good",]$Direction="MPR"
    if ("Bad" %in% tbl$Direction) {
      tbl[tbl$Direction=="Bad",]$Direction="Non MPR"
    }
  } else {
    tbl[tbl$Direction=="pR",]$Direction="MPR"
    tbl[tbl$Direction=="pNR",]$Direction="Non MPR"
  }
  ss <- tableGrob(tbl, rows=NULL)
  roc <- median_roc(rds, offset = offset,  good = good)
  if(!is.null(topFeatures)) {
    featg <- feature_plot(topFeatures)
    featr <- feature_rocs(topFeatures)
    grid.arrange(roc, ss, featg, featr,
                 top = textGrob(title,
                                gp=gpar(fontsize=20)))
  } else {
    grid.arrange(roc, ss,
                 top = textGrob(title,
                                gp=gpar(fontsize=20)))
  }
}

expvec <- readRDS("expvec_CpG_GitHub.rds")
samples <- readRDS("Samples_Intermediate.rds")

# limit to intermediate samples
for(experiment in names(expvec)) {
  expvec[[experiment]] <- expvec[[experiment]][,samples] %>%
    filter_experiment(featcutoff = c(5, 5))
}

rowData(expvec$LKT)$Description <- 
  str_c("Gram ", rowData(expvec$LKT)$Gram)

# loop to run analyses
NUM <- 100
for(experiment in c("ECNumber", "Product","Interpro", "GO", "LKT")) {
  exp <- expvec[[experiment]]
  matrix <- list()
  accuracy <- c()
  auc <- c()
  obs <- list()
  preds <- list()
  probs <- list()
  imp <- list()
  pvalue <- c()
  method <- "rf"
  filename <- str_c(experiment, "_", method, "_fs_70_data.rds")
  if (!file.exists(filename)) {
    for(seed in 1:NUM) {
      set.seed(seed)
      index <- sample(ncol(exp), round(ncol(exp)*0.7))
      if (length(unique(exp[,index]$Response_at_surgery_dichot))==1 ||
          length(unique(exp[,-index]$Response_at_surgery_dichot))==1) {
        index <- sample(ncol(exp), round(ncol(exp)*0.71))
      }
      flog.info(seed)
      results <-
        suppressMessages(suppressWarnings(runML(exp, "Response_at_surgery_dichot", 
                                                "Good", index, method, filterFeature = TRUE)))
      remaining <- exp[,-index]
      preds[[seed]] <- predict(results$models, t(assay(remaining)))
      probs[[seed]] <- predict(results$models, t(assay(remaining)), type="prob")
      obs[[seed]] <-
        factor(remaining$Response_at_surgery_dichot, levels = levels(preds[[seed]]))
      matrix[[seed]] <- confusionMatrix(preds[[seed]], obs[[seed]], positive = "Good")
      pvalue[seed] <- matrix[[seed]]$overall[[6]]
      accuracy[seed] <-  matrix[[seed]]$overall[[1]]
      auc[seed] <- as.numeric(roc(response=obs[[seed]], predictor=probs[[seed]][,"Good"])$auc)
      if (method=="rf") {
        imp[[seed]] <- varImp(results$models)$importance %>% arrange(-Overall) %>% 
          data.frame()
      } else if (method=="svmPoly") {
        imp[[seed]] <- varImp(results$models)$importance %>% arrange(-Good) %>% 
          data.frame()
      }
      prCounts <- rowMeans(assay(exp)[row.names(imp[[seed]]),exp$Response_at_surgery_dichot=="Good"])
      pnrCounts <-  rowMeans(assay(exp)[row.names(imp[[seed]]),exp$Response_at_surgery_dichot=="Bad"])
      imp[[seed]]$direction <- ifelse(prCounts > pnrCounts, "Good", "Bad")
      names(obs[[seed]]) <- remaining$Patient_ID
      names(preds[[seed]]) <- remaining$Patient_ID
    }
    saveRDS(list(accuracy, auc, matrix, imp, obs, preds, probs, pvalue),
            filename)
    topfile <- str_c(experiment, "_topFeatures.rds")
    if (!file.exists(topfile)) {
      tops <- runMultiTopFeatures(exp, response = "Response_at_surgery_dichot", positive = "Good", 
                            imp = summarize_imp(str_c(experiment, "_", method, "_fs_70_data.rds"), top = 200), 
                            top = c(3,5,10,20,30,50,70,100,150,200))
      saveRDS(tops, topfile)
    }
  }
}

fname <- str_c("ml_results_", format(Sys.time(), "%Y-%m-%d"), ".pdf")

pdf(fname, width = 16, height = 15)
for(experiment in c("ECNumber", "Product","Interpro", "GO", "LKT")) {
  plot_ml_fig(str_c(experiment, "_rf_fs_70_data.rds"), 
              description=rowData(expvec[[experiment]])[,"Description", drop=FALSE],
              good="Good", title = str_c(experiment, " 70/30 split"),
              topFeatures = str_c(experiment, "_topFeatures.rds"))
}
dev.off()
```
