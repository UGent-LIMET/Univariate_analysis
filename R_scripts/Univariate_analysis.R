# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part II: Univariate analysis




##########R Pipeline - Part II: univariate analysis##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part II: univariate analysis - start!")
# Part II: statistical analysis

## data_loading
setwd(path_data_in)
COLLUMN_NR_COMPARISON1 <- COLLUMN_NR_LAST_BEFORE_COMPARISONS + 1        #column number of the column 'Comparison1' for pairwise comparison (OPLSDA)
COLLUMN_NR_MULTIPLE_COMPARISON1 <- COLLUMN_NR_COMPARISON1 + AMOUNT_OF_COMPARISONS         #column number of the column 'MultipleComparison1' for multiple comparison (PLSDA, LIMMA)
COLLUMN_NR_PROJECTION1 <- COLLUMN_NR_MULTIPLE_COMPARISON1 + AMOUNT_OF_MULTIPLE_COMPARISONS
COLLUMN_NR_START_VARIABLES <- COLLUMN_NR_PROJECTION1 + AMOUNT_OF_PROJECTIONS

sampleMetadata <- load_sampleMetadata(INPUT_SAMPLES)
#Add "X" to sampleNames
sampleMetadata[,1] <- paste(replicate(nrow(sampleMetadata),"X"), sampleMetadata[,1], sep="")
check_nrow_SM <- nrow(sampleMetadata)

if (INPUT_VARIABLES == VARIABLEMETADATA_FROM_PIPELINE){
  INPUT_VARIABLES <- paste(name_project, '_variableMetadata.txt', sep="")
  if(exists("COLLUMN_NR_START_SAMPLES") == FALSE){ 		#if after merge, will be given value 21, so do not change
    COLLUMN_NR_START_SAMPLES <- 20  #always 20 (auto and manual must be same format)
  }
}
variableMetadata <- load_variableMetadata(INPUT_VARIABLES)
#note: automatically adds X to samples at this point


## set directory to output
setwd(path_data_out)


## merge variables from variableMetadata into the sampleMetadata
### make sampleMatrix
library(data.table)
variableMetadata_from_start_samples <- subset(variableMetadata, select = -c(2:(COLLUMN_NR_START_SAMPLES-1))) #remove info exept CompID
sampleMatrix <- transpose(variableMetadata_from_start_samples)
colnames(sampleMatrix) <- variableMetadata_from_start_samples[ ,1]
sampleMatrix <- as.data.frame(sapply(sampleMatrix, as.numeric))
sampleMatrix$SampleName <- colnames(variableMetadata_from_start_samples)
sampleMatrix <- sampleMatrix[-1,] #remove first row (colnames compids + samplename; so 1 extra col than no of variables), is captured in colnames
#write_dataframe_as_txt_file(sampleMatrix, 'sampleMatrix.txt')

#merge variable intensities from compIDs to correct sample (if order not same), output order is sorted by samplename
sampleMetadata <- merge_accoding_to_SampleName(sampleMetadata, sampleMatrix)
write_dataframe_as_txt_file(sampleMetadata, 'sampleMetadata_variableMatrix_merged.txt')


##Check sampleMatrx 
#see if same number of rows (amount of sampleNames) as when loaded prev, so merge was succesfull
if(nrow(sampleMetadata) != check_nrow_SM){
  stop("ERROR: univariate analysis stopped because SampleMetadata and VariableMetadata are incompatible to merge into correct sampleMatrix.")
}

##Check Metabolitenames are present
if(!("MetaboliteName" %in% colnames(variableMetadata)) == TRUE){ #if this column not present, stop script
  stop("ERROR: univariate analysis stopped because no column 'MetaboliteName' present in variableMetadata.")
}

## hierarchical clustering heat map
# heatmap based on hierarchical clustering (w all samples, wo QCs)
samples_metadata <- sampleMetadata[sampleMetadata$Type == 'Sample',] 
samples_matrix <- from_df_to_matrix(samples_metadata)

# add names from samples and variables (nice name if present) + check ok
if("Name" %in% colnames(samples_metadata)){ 		#if alternative name for samples exists, use these on heat map
  stopifnot(length(unique(samples_metadata$Name)) == nrow(samples_metadata)) #must be unique name in Name column
  rownames(samples_matrix) <- samples_metadata$Name
}else(rownames(samples_matrix) <- samples_metadata$SampleName)
colnames(samples_matrix) <- variableMetadata$MetaboliteName


## heatmap plot with minimal distance
#https://stackoverflow.com/questions/25528059/cluster-data-in-heat-map-in-r-ggplot
#clustering samples (rows) order
ord <- hclust(dist(samples_matrix))$order #info about order in object ord$order
#default clustering: Cluster method: complete; Distance: euclidean see doc: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
#ord
dendrogram <-as.dendrogram(hclust(dist(samples_matrix)))
png(paste(name_project, "_dendogram_clust", ".png", sep=""), width=10, height=5, units="in", res=300)
par(cex=.2)
plot(dendrogram)
dev.off()

#2nd clustering for the variables (not used on HM)
dendrogram2 <-as.dendrogram(hclust(dist(t(samples_matrix))))
png(paste(name_project, "_dendogram_clust_variables", ".png", sep=""), width=7, height=5, units="in", res=150)
par(cex=.5)
plot(dendrogram2)
dev.off()

#long format for ggplot
library("reshape") 
samples_matrix$SampleName <- rownames(samples_matrix) #depending on which name (system or nic name), as chosen above
data_melt <- melt(samples_matrix, id = c("SampleName")) 

#info from dendogram, for order in ggplot
comp_names2 <- NULL
comp_names2$SampleName <- rownames(samples_matrix) #idem accoding to name chosen above
comp_names2 <- data.frame(comp_names2)
comp_names2$ord <- ord

#merge info long format + order clustering accoding to name chosen above
data_melt <- merge(data_melt, comp_names2, by = "SampleName")
data_melt$order_plot <- data_melt$ord #replace order according to dist
data_melt$order_plot <- as.numeric(data_melt$order_plot)

#factor in order of hdist clustering to SampleName, see in plot function
#plot using info clust
heatmap_hclust <- plot_heatmap_hclust(data_melt)
ggsave(filename = paste(name_project, "_heatmap_clust", ".png", sep=""), heatmap_hclust, width = 7, height = 7, dpi = 300)
ggsave(filename = paste(name_project, "_heatmap_clust", ".pdf", sep=""), heatmap_hclust, width = 7, height = 7, dpi = 300)

#log heatmap based on hierarchical clustering (only on 1 comparison)
data_melt_log <- data_melt
data_melt_log$value <- log10(data_melt$value + 1) #log10(x+1)

heatmap_hclust <- plot_heatmap_hclust(data_melt_log)
ggsave(filename = paste(name_project, "_heatmap_clust", "LOG.png", sep=""), heatmap_hclust, width = 7, height = 7, dpi = 300)
ggsave(filename = paste(name_project, "_heatmap_clust", "LOG.pdf", sep=""), heatmap_hclust, width = 7, height = 7, dpi = 300)



########## univariate analysis comparison ##########
## for each comparison

nsamples_metadata <- sampleMetadata

#report t-test p-values for each pairwise comparison and each variable
name_report_ttest <- paste(name_project,'_pvalues_Comparisons.txt', sep="")
line0 <- paste(c("Pairwise comparison", "Standard", "Hypothesis test", "Group Nrs" ,"Amount per group", "Mean per group", "Median per group", "SD per group", "p-value", "Fold change"), collapse= '\t')
append_result_to_report(line0, name_report_ttest)
q_values <- NULL #adjusted p-value, calculated per comparison but written to report after loop

if(AMOUNT_OF_COMPARISONS >= 1){
  for(pairwise_comparison in 1:AMOUNT_OF_COMPARISONS){
    
    #pairwise_comparison <- 1 #for testing
    print(paste("start calculation comparison ", pairwise_comparison, sep=""))
    
    
    ### prepare matrix    
    comp_ <- nsamples_metadata[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1]
    stopifnot(class(comp_)=='integer') #need only 2 possible {-1,1} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    variable.names2 <- variableMetadata$MetaboliteName 

    #make sure no NAs present in matrix  
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x))) #NAN -> NA
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))] #del variables with all 0's
    colnames(samples_matrix_comp_no0) <- variable.names2 #no 'x', "."
    
    #remove if all values are 0 for each metabolite
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,apply(samples_matrix_comp_no0, 2, function(x) sum(x) != 0)]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1])
    
    
    
    ### Heat maps (ggplot)
    
    #wo transpose for ggplot heatmap:
    Tsamples_matrix_comp_no0_selection <- samples_matrix_comp_no0 #matrix with comp retained, <lod as NA and colnames = standards
    
    # add names from samples and variables (nice name if present)
    if("Name" %in% colnames(samples_metadata_comp)){ 		#if alternative name for samples exists, use these on heat map
      rownames(Tsamples_matrix_comp_no0_selection) <- samples_metadata_comp$Name
    }else(rownames(Tsamples_matrix_comp_no0_selection) <- samples_metadata_comp$SampleName)
    colnames(Tsamples_matrix_comp_no0_selection) <- colnames(samples_matrix_comp_no0)
    
    
    library("reshape") 
    #head(Tsamples_matrix_comp_no0_selection)
    Tsamples_matrix_comp_no0_selection$SampleName <- rownames(Tsamples_matrix_comp_no0_selection)
    data_melt <- melt(Tsamples_matrix_comp_no0_selection, id = c("SampleName")) 
    
    #write copy for 'manual' editing
    name_df <- paste(name_project,'_TableHeatMap_comparison', pairwise_comparison, '.txt', sep="")
    write_dataframe_as_txt_file(Tsamples_matrix_comp_no0_selection, name_df) #samplenames present in last column
    
    
    #order accord to group for heatmap in ggplot
    comp <- as.character(samples_metadata_comp[,COLLUMN_NR_COMPARISON1+pairwise_comparison-1])

    if("Name" %in% colnames(samples_metadata_comp)){ 		#if alternative name for samples exists, use these on heat map
      sn <- as.character(samples_metadata_comp$Name)
    }else(sn <- as.character(samples_metadata_comp$SampleName))
    comp_names <- as.data.frame(cbind(sn, comp))
    comp_names <- comp_names[order(comp_names$comp, decreasing = FALSE),]
    SampleName_reorder <- comp_names$sn #variable needed for ggplot correct order samplenames
    SampleName_reorder <- as.character(SampleName_reorder)

    comp_names$SampleName <- comp_names$sn #no SampleName as variable, but with this name for merge
    comp_names$SampleName <- as.character(comp_names$SampleName)
    comp_names$order_plot <- c(1:(nrow(comp_names))) #for correct HM order!! must start w 1 or shift
    
    data_melt <- merge(data_melt, comp_names, by = "SampleName")
    #data_melt$order_plot <- as.numeric(data_melt$order_plot) #ipv int?
    
    
    #heatmap plot with grouping comparison above on plot
    heatmap_comp <- plot_heatmap_w_group(data_melt)
    ggsave(filename = paste(name_project, "_heatmap_group_comparison" , pairwise_comparison, ".png", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    ggsave(filename = paste(name_project, "_heatmap_group_comparison" , pairwise_comparison, ".pdf", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    
    
    ##log heatmap plot
    data_melt_log <- data_melt
    data_melt_log$value <- log10(data_melt$value + 1) #log10(x+1)
    
    heatmap_comp <- plot_heatmap_w_group(data_melt_log)
    ggsave(filename = paste(name_project, "_heatmap_group_comparison" , pairwise_comparison, "LOG.png", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    ggsave(filename = paste(name_project, "_heatmap_group_comparison" , pairwise_comparison, "LOG.pdf", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    
    
    
    ### Boxplots for each variable 
    #https://statisticsglobe.com/boxplot-in-r
    #Boxplots are a popular type of graphic that visualize the minimum non-outlier, the first quartile, the median, the third quartile, and the maximum non-outlier of numeric data in a single plot
    
    standards <- ncol(samples_matrix_comp_no0)
    for(standard_nr in 1:standards){
      #standard_nr <- 1 #for testing 1, 6, 8, 21, 31 in comp2 (3 group), 12&13 (same input table: ok, same in vm), 108&109 in comp2
      #print(standard_nr)
      
      name_standard <- colnames(samples_matrix_comp_no0[standard_nr])
      
      #also name to write boxplot.png (so no "/" eg)
      library("stringr") 
      name_standardWOsymbol <- colnames(samples_matrix_comp_no0[standard_nr])
      name_standardWOsymbol <- str_replace_all(name_standardWOsymbol, "[^[:alnum:]]", "")    # Delete non-alphanumeric
      
      ### hypothesis testing 
      
      ## check normality for each group in comparison
      comp_variable <- cbind(samples_matrix_comp_no0[standard_nr], comp)
      colnames(comp_variable) <- c("variable", "comp")
      
      #also from here extract addit info for table
      names_groups <- NULL
      amount_groups <- NULL
      mean_groups <- NULL
      median_groups <- NULL
      sd_groups <- NULL
      
      #calculate normality per group in comp
      normality_variables <- NULL
      for(nr in sort(unique(comp))){ #run over groups in comp -1 and +1
        #print(nr)
        
        y1 <- comp_variable[comp_variable$comp == nr,] 
        y1 <- y1[,1]
        #hist(y1)
        
        #The Shapiro-Wilk Test tests the null hypothesis that the samples come from a normal distribution vs. the alternative hypothesis that the samples do not come from a normal distribution. 
        #In this case, the p-value of the test is 0.005999, which is less than the alpha level of 0.05. This suggests that the samples do not come a normal distribution.
        #https://www.statology.org/anova-assumptions/
		    #shapiro.test(rnorm(100, mean = 5, sd = 3))
	    	#shapiro.test(runif(100, min = 2, max = 4))

        #! if keep normal => error if all intensit are SAME in 1 group (all both groups 0 are not present in samples_matrix_comp_no0)
        #if 1 group all are SAME value: ! "Error in shapiro.test(test) : all 'x' values are identical" => normal OK (very narrow distribution inf small)
        tryCatch(
          expr = {
            y1_normality_variable <- shapiro.test(y1)   #p < 0.05 --> NO normal distribution
            normality_variable <- as.numeric(unlist(y1_normality_variable[2]))
            #print(normality_variable)
          },
          error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
            print(paste0("cannot calculate normality for ", name_standard, " in group ", nr, " of comparison ", pairwise_comparison))
            normality_variable <- 10000000 #if all values same => NO normal distribution, use non-parametric test
          }
        )
        
        #append all normality values
        normality_variables <- c(normality_variables, normality_variable)
        #print(normality_variables)
        
        #add info per group
        names_groups <- c(names_groups, nr)
        amount_groups <- c(amount_groups, length(y1))
        mean_groups <- c(mean_groups, format(round(mean(y1),digits=3),nsmall=3))
        median_groups <- c(median_groups, format(round(median(y1),digits=3),nsmall=3))
        sd_groups <- c(sd_groups, format(round(sd(y1),digits=3),nsmall=3))
      }
      
      #info in 1 output cell
      names_groups <- paste(unlist(names_groups), collapse='; ')
      amount_groups <- paste(unlist(amount_groups), collapse='; ')
      mean_groups <- paste(unlist(mean_groups), collapse='; ')
      median_groups <- paste(unlist(median_groups), collapse='; ')
      sd_groups <- paste(unlist(sd_groups), collapse='; ')
      fold_changes <- ( (as.numeric(strsplit(mean_groups,"; ", fixed=T)[[1]][1]) +0.00000001) / (as.numeric(strsplit(mean_groups,"; ", fixed=T)[[1]][2])+0.00000001) ) #mean ration fold change
        
      #manual with +/-1 select  #todo make work if needed for other label 
      y1 <- comp_variable[comp_variable$comp == sort(unique(comp))[1],]   #group "-1"
      y1 <- y1[,1]
      
      y2 <- comp_variable[comp_variable$comp == sort(unique(comp))[2],]   #group "+1"
      y2 <- y2[,1]
      
      
    
      ## ttest if normal distribution OK
      if(all(normality_variables >= 0.05)){

        ttest_univar_comp <- t.test(x=y1, y=y2,
                                    alternative="two.sided", #2 sided since diff can be either higher/lower
                                    mu=0,
                                    paired = FALSE,   #unpaired set
                                    var.equal = FALSE, #auto select Welch instead of checking => see boxplots is diff var
                                    conf.level = 0.95) 
        # p<0.05 --> significant verwerpt H0 (diff between 2 groups)
        p_value <- unlist(ttest_univar_comp[3])
        
        #name_test
        name_test <- unlist(ttest_univar_comp$method)
        
        #report t-test p-values for each pairwise comparison and each variable (add here)
        line_ <- paste(c(pairwise_comparison, name_standard, name_test, names_groups, amount_groups, mean_groups, median_groups, sd_groups, p_value, fold_changes), collapse='\t')
        append_result_to_report(line_, name_report_ttest)
        
        if(is.na(p_value)){
          p_value <- 10000000 #if intensity in both gropus = 0 => means = 0, no diff in means so p-value is inf
          #actually not needed anymore since rm from matrix_no_0 but leave in code
        }
        if(p_value < 0.05){
          boxplot_comp <- plot_duoboxplot(samples_matrix_comp_no0)
          
          png(paste(name_project, "_boxplot_comp" , pairwise_comparison, "_std_", name_standardWOsymbol, "_SIGNIF.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
        if(p_value >= 0.05){
          ## boxplots also NS stds plot
          boxplot_comp <- plot_duoboxplot(samples_matrix_comp_no0)
          png(paste(name_project, "_boxplot_comp" , pairwise_comparison, "_std_", name_standardWOsymbol, "_NS.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
      }
      
      ## Mann-Whitney U test if normal distribution NOK
      if(any(normality_variables < 0.05)){
        ##Wilcoxon Rank Sum Test == Mann-Whitney U test
        wicox_univar_comp <- wilcox.test(x=y1, y=y2,
                                         alternative="two.sided", #2 sided since diff can be either higher/lower
                                         mu=0,
                                         paired = FALSE,   #unpaired set
                                         exact = NULL, correct = TRUE,
                                         conf.level = 0.95) 
        # p<0.05 --> significant verwerpt H0 (diff between 2 groups)
        p_value <- unlist(wicox_univar_comp[3])
        
        #name_test
        name_test <- unlist(wicox_univar_comp$method)
        
        #report wilcoxon p-values for each pairwise comparison and each variable (add here)
        line_ <- paste(c(pairwise_comparison, name_standard, name_test, names_groups, amount_groups, mean_groups, median_groups, sd_groups, p_value, fold_changes), collapse='\t')
        append_result_to_report(line_, name_report_ttest)
        
        if(is.na(p_value)){
          p_value <- 10000000 #if intensity in both gropus = 0 => means = 0, no diff in means so p-value is inf
          #actually not needed anymore since rm from matrix_no_0 but leave in code
        }
        if(p_value < 0.05){
          boxplot_comp <- plot_duoboxplot(samples_matrix_comp_no0)
          
          png(paste(name_project, "_boxplot_comp" , pairwise_comparison, "_std_", name_standardWOsymbol, "_SIGNIF.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
        if(p_value >= 0.05){
          ## boxplots also NS stds plot
          boxplot_comp <- plot_duoboxplot(samples_matrix_comp_no0)
          png(paste(name_project, "_boxplot_comp" , pairwise_comparison, "_std_", name_standardWOsymbol, "_NS.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
      }
    }
    
    ### Volcano plot (using mean ratio fold change) + calc adj p-value
    volcano_df <- read.table(name_report_ttest, header=TRUE, sep="\t") 
    volcano_df <- volcano_df[volcano_df$Pairwise.comparison == pairwise_comparison,]
    q <- p.adjust(volcano_df$p.value, method = "fdr")
    q_values <- c(q_values, q)
    volcano_df$Adjusted_p_value <- q
    
    plot_volcano_comp <- plot_volcano(volcano_df)
    png(paste(name_project, "_volcano_comp" , pairwise_comparison, ".png", sep=""), width=5, height=5, units="in", res=150)
    plot(plot_volcano_comp)
    dev.off() 
    
  }
  
  #### After all pairwise comparisons done, write adj p-values to report
  report_pairwise <- read.table(name_report_ttest, header=TRUE, sep="\t")
  report_pairwise$Adjusted_p_value <- q_values
  name_report_ttest2 <- paste(name_project,'_adjpvalues_Comparisons.txt', sep="")
  write_dataframe_as_txt_file(report_pairwise, name_report_ttest2) 
  
}



###################




########## univariate analysis MultiComparison ##########

nsamples_metadata <- sampleMetadata

#report ANOVA p-values for each multiple comparison and each variable
name_report_ANOVA <- paste(name_project,'_pvalues_MultipleComparisons.txt', sep="")
line0 <- paste(c("Multiple comparison", "Standard", "Hypothesis test", "Groups Nrs" ,"Amount per group", "Mean per group", "Median per group", "SD per group", "p-value"), collapse = '\t')
append_result_to_report(line0, name_report_ANOVA)
q_values <- NULL #adjusted p-value, calculated per comparison but written to report after loop

#report posthoc after anova if signif
name_report_posthoc <- paste(name_project,'_TukeyHSD_MultipleComparisons.txt', sep="")
line0 = paste(c("Multiple comparison", "Standard", "comp", "diff", "lwr", "upr", "p adj"), collapse= '\t')
append_result_to_report(line0, name_report_posthoc)

#report posthoc after Krustal if signif
name_report_posthocDunn <- paste(name_project,'_Dunn_MultipleComparisons.txt', sep="")
line0 = paste(c("Multiple comparison", "Standard", "comp", "chi2", "Z", "altP", "altP.adjusted"), collapse= '\t')
append_result_to_report(line0, name_report_posthocDunn)


## for each multiple comparison
if(AMOUNT_OF_MULTIPLE_COMPARISONS >= 1){
  for(multiple_comparison in 1:AMOUNT_OF_MULTIPLE_COMPARISONS){
    
    #multiple_comparison <- 1
    
    print(paste("start calculation multiple comparison ", multiple_comparison, sep=""))
    
    ### remove non-essential data for comparison
    #remove samples that are not included in comp (Na instead of 0/1)
 
    comp_ <- nsamples_metadata[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1]
    stopifnot(class(comp_)=='integer') #need only 2 possible {-1,1} arguments, rest is NA and not ""!
    samples_metadata_comp <- nsamples_metadata[!is.na(comp_), ]    
    
    #remove variables that do not occur in comparison (not detected in both groups); value NaN after scaling
    samples_matrix_comp <- from_df_to_matrix(samples_metadata_comp)
    variable.names2 <- variableMetadata$MetaboliteName 
    
    #make sure no NAs present in matrix  
    samples_matrix_comp_no0 <- data.frame(sapply(samples_matrix_comp, function(x) ifelse(is.nan(x), NA, x))) #NAN -> NA
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,which(unlist(lapply(samples_matrix_comp_no0, function(x)!all(is.na(x)))))] #del variables with all 0's
    colnames(samples_matrix_comp_no0) <- variable.names2 #no 'x', "."
    
    #remove if all values are 0 for each metabolite
    samples_matrix_comp_no0 <- samples_matrix_comp_no0[,apply(samples_matrix_comp_no0, 2, function(x) sum(x) != 0)]
    
    comp <- as.factor(samples_metadata_comp[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1])
    
    
    
    ### Heat maps (ggplot)
    #wo transpose for ggplot heatmap:
    Tsamples_matrix_comp_no0_selection <- samples_matrix_comp_no0 #matrix with comp retained, <lod as NA and colnames = standards
    
    # add names from samples and variables (nice name if present)
    if("Name" %in% colnames(samples_metadata_comp)){ 		#if alternative name for samples exists, use these on heat map
      rownames(Tsamples_matrix_comp_no0_selection) <- samples_metadata_comp$Name
    }else(rownames(Tsamples_matrix_comp_no0_selection) <- samples_metadata_comp$SampleName)
    colnames(Tsamples_matrix_comp_no0_selection) <- colnames(samples_matrix_comp_no0)
    
    
    library("reshape") 
    #head(Tsamples_matrix_comp_no0_selection)
    Tsamples_matrix_comp_no0_selection$SampleName <- rownames(Tsamples_matrix_comp_no0_selection)
    data_melt <- melt(Tsamples_matrix_comp_no0_selection, id = c("SampleName")) 
    
    #write copy for 'manual' editing
    name_df <- paste(name_project,'_TableHeatMap_Mcomparison', multiple_comparison, '.txt', sep="")
    write_dataframe_as_txt_file(Tsamples_matrix_comp_no0_selection, name_df) #samplenames present in last column
    
    
    #order accord to group for heatmap in ggplot
    comp <- as.character(samples_metadata_comp[,COLLUMN_NR_MULTIPLE_COMPARISON1+multiple_comparison-1])
    
    if("Name" %in% colnames(samples_metadata_comp)){ 		#if alternative name for samples exists, use these on heat map
      sn <- as.character(samples_metadata_comp$Name)
    }else(sn <- as.character(samples_metadata_comp$SampleName))
    comp_names <- as.data.frame(cbind(sn, comp))
    comp_names <- comp_names[order(comp_names$comp, decreasing = FALSE),]
    SampleName_reorder <- comp_names$sn #variable needed for ggplot correct order samplenames
    SampleName_reorder <- as.character(SampleName_reorder)
    
    comp_names$SampleName <- comp_names$sn #no SampleName as variable, but with this name for merge
    comp_names$SampleName <- as.character(comp_names$SampleName)
    comp_names$order_plot <- c(1:(nrow(comp_names))) #for correct HM order!! must start w 1 or shift
    
    data_melt <- merge(data_melt, comp_names, by = "SampleName")
    
    ## heatmap according to group
    heatmap_comp <- plot_heatmap_w_group(data_melt)
    ggsave(filename = paste(name_project, "_heatmap_group_Mcomparison" , multiple_comparison, ".png", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    ggsave(filename = paste(name_project, "_heatmap_group_Mcomparison" , multiple_comparison, ".pdf", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    
    ##log heatmaps plots
    data_melt_log <- data_melt
    data_melt_log$value <- log10(data_melt$value +1) #log10(x +1)
    
    heatmap_comp <- plot_heatmap_w_group(data_melt_log)
    ggsave(filename = paste(name_project, "_heatmap_group_Mcomparison" , multiple_comparison, "LOG.png", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    ggsave(filename = paste(name_project, "_heatmap_group_Mcomparison" , multiple_comparison, "LOG.pdf", sep=""), heatmap_comp, width = 7, height = 7, dpi = 300)
    
    

    
    ## Boxplot for each variable 
    #https://statisticsglobe.com/boxplot-in-r
    #Boxplots are a popular type of graphic that visualize the minimum non-outlier, the first quartile, the median, the third quartile, and the maximum non-outlier of numeric data in a single plot
    
    standards <- ncol(samples_matrix_comp_no0)
    for(standard_nr in 1:standards){
      #standard_nr <- 1 #for testing 1, 6
      #print(standard_nr)
      
      name_standard <- colnames(samples_matrix_comp_no0[standard_nr])
      
      #also name to write boxplot.png (so no "/" eg)
      library("stringr") 
      name_standardWOsymbol <- colnames(samples_matrix_comp_no0[standard_nr])
      name_standardWOsymbol <- str_replace_all(name_standardWOsymbol, "[^[:alnum:]]", "")    # Delete non-alphanumeric

      
      ## check normality for each group in comparison
      comp_variable <- cbind(samples_matrix_comp_no0[standard_nr], comp)
      
      #also from here extract addit info for table
      names_groups <- NULL
      amount_groups <- NULL
      mean_groups <- NULL
      median_groups <- NULL
      sd_groups <- NULL
      
      #calculate normality per group in comp
      normality_variables <- NULL
      for(nr in sort(unique(comp))){ #run over groups in Mcomp 1,2,3,...
        #print(nr)
        
        y1 <- comp_variable[comp_variable$comp == nr,] 
        y1 <- y1[,1]
        #hist(y1)
        
        #The Shapiro-Wilk Test tests the null hypothesis that the samples come from a normal distribution vs. the alternative hypothesis that the samples do not come from a normal distribution. 
        #In this case, the p-value of the test is 0.005999, which is less than the alpha level of 0.05. This suggests that the samples do not come a normal distribution.
        #https://www.statology.org/anova-assumptions/
        
        #! if keep normal => error if all intensit are SAME in 1 group (all both groups 0 are not present in samples_matrix_comp_no0)
        #if 1 group all are SAME value: ! "Error in shapiro.test(test) : all 'x' values are identical" => normal OK (very narrow distribution inf small)
        tryCatch(
          expr = {
            y1_normality_variable <- shapiro.test(y1)   #p < 0.05 --> NO normal distribution
            normality_variable <- as.numeric(unlist(y1_normality_variable[2]))
            #print(normality_variable)
          },
          error = function(e){ 
            # (Optional)
            # Do this if an error is caught...
            print(paste0("cannot calculate normality for ", name_standard, " in group ", nr, " of multiple comparison ", multiple_comparison))
            normality_variable <- 10000000 #if all values same => NO normal distribution, use non-parametric test
          }
        )

        #append all normality values
        normality_variables <- c(normality_variables, normality_variable)
        #print(normality_variables)
        
        #add info per group
        names_groups <- c(names_groups, nr)
        amount_groups <- c(amount_groups, length(y1))
        mean_groups <- c(mean_groups, format(round(mean(y1),digits=3),nsmall=3))
        median_groups <- c(median_groups, format(round(median(y1),digits=3),nsmall=3))
        sd_groups <- c(sd_groups, format(round(sd(y1),digits=3),nsmall=3))
      }
      
      #info in 1 output cell
      names_groups <- paste(unlist(names_groups), collapse='; ')
      amount_groups <- paste(unlist(amount_groups), collapse='; ')
      mean_groups <- paste(unlist(mean_groups), collapse='; ')
      median_groups <- paste(unlist(median_groups), collapse='; ')
      sd_groups <- paste(unlist(sd_groups), collapse='; ')
      
      
      ## one way ANOVA test if normal distribution OK
      if(all(normality_variables >= 0.05)){
        ## one way ANOVA
        # https://www.scribbr.com/statistics/anova-in-r/
        # one-way since 1 independ vaiable (= groep VS intensity)
        
        #One-way ANOVA example
        #In the one-way ANOVA, we test the effects of 3 types of fertilizer on crop yield.
        #Two-way ANOVA example
        #In the two-way ANOVA, we add an additional independent variable: planting density. We test the effects of 3 types of fertilizer and 2 different planting densities on crop yield.
        
        comp_variable <- cbind(samples_matrix_comp_no0[standard_nr], comp)
        colnames(comp_variable) <- c("variable", "comp")
        
        model <- aov(formula = variable ~ comp, data = comp_variable) 
        #summary(model)
        sum_test = unlist(summary(model))
        p_value <- sum_test["Pr(>F)1"]
        
        #name_test
        name_test <- "One-way ANOVA"
        
        #report ANOVA p-values for each multiple comparison and each variable (add here)
        line_ <- paste(c(multiple_comparison, name_standard, name_test, names_groups, amount_groups, mean_groups, median_groups, sd_groups, p_value), collapse = '\t')
        append_result_to_report(line_, name_report_ANOVA)
        
        
        if(is.na(p_value)){
          p_value <- 10000000 #if intensity in both gropus = 0 => means = 0, no diff in means so p-value is inf
          #actually not needed anymore since rm from matrix_no_0 but leave in code
        }
        if(p_value < 0.05){
          ##boxplot
          comp <- as.numeric(comp) #to nr (only works from char->num, no factors!)
          comp <- sprintf("%02d",comp) #to txt with 01 so no issue order if >9groups
          
          boxplot_comp <- plot_multiboxplot(samples_matrix_comp_no0)
          
          png(paste(name_project, "_boxplot_Mcomp" , multiple_comparison, "_std ", name_standardWOsymbol, "_SIGNIF.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
          
          
          ##post-hoc test 
          #  => to retrieve which group is responsible for sign diff ANOVA
          #  https://www.r-bloggers.com/2020/10/anova-in-r/
          #  https://www.rdocumentation.org/packages/stats/versions/3.3/topics/TukeyHSD
          #  https://www.r-bloggers.com/2018/09/tukeys-test-for-post-hoc-analysis/
          #  use Tukey HSD test => is used to compare all groups to each other (so all possible comparisons of 2 groups), 
          #  because in pipeline no info 'ref' group to compare against (e.g. Dunnett's test)
          
          # Tukey HSD test:
          post_test <- TukeyHSD(x=model, which="comp", ordered = FALSE, conf.level = 0.95)
          
          #report posthoc p-adjusted values (+ all info) for each multiple comparison and each variable (add here)
          linePH <- as.data.frame(post_test$comp)
          linePH$multicomp <- multiple_comparison
          linePH$standard <- name_standard
          linePH$comp <- rownames(linePH) 
          linePH <- linePH[ , c("multicomp", "standard", "comp", "diff", "lwr", "upr", "p adj")]     #reorer cols
          append_dataframe_to_report(linePH, name_report_posthoc)
  
          #plot visual Tukey Honest Significant Differences
          #(ku hier wel zeggen: mean sign diff! want same disctrib: normal, zie onder)
          png(paste(name_project, "_boxplot_Mcomp" , multiple_comparison, "_std ", name_standardWOsymbol, "_SIGNIF_posthoc.png", sep=""), width=7, height=5, units="in", res=150)
          plot(TukeyHSD(model, "comp"))
          dev.off()
  
        }
        if(p_value >= 0.05){
          comp <- as.numeric(comp) #to nr (only works from char->num, no factors!)
          comp <- sprintf("%02d",comp) #to txt with 01 so no issue order if >9groups
          
          ## boxplots also NS stds plot
          boxplot_comp <- plot_multiboxplot(samples_matrix_comp_no0)
          
          png(paste(name_project, "_boxplot_Mcomp" , multiple_comparison, "_std_", name_standardWOsymbol, "_NS.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
      }
        
      
      ## Kruskal-Wallis Rank Sum Test if normal distribution NOK
      if(any(normality_variables < 0.05)){
        
        comp_variable <- cbind(samples_matrix_comp_no0[standard_nr], comp)
        colnames(comp_variable) <- c("variable", "comp")
        
        #Kruskal-Wallis Rank Sum Test
        kruskal_univar_comp <- kruskal.test(formula = variable ~ comp, data = comp_variable)

        # p<0.05 --> significant verwerpt H0 (diff between 2 groups)
        p_value <- unlist(kruskal_univar_comp[3])
        
        #name_test
        name_test <- unlist(kruskal_univar_comp$method)
        
        #report Kruskal p-values for each pairwise comparison and each variable (add here)
        line_ <- paste(c(multiple_comparison, name_standard, name_test, names_groups, amount_groups, mean_groups, median_groups, sd_groups, p_value), collapse='\t')
        append_result_to_report(line_, name_report_ANOVA)
        
        if(is.na(p_value)){
          p_value <- 10000000 #if intensity in both gropus = 0 => means = 0, no diff in means so p-value is inf
          #actually not needed anymore since rm from matrix_no_0 but leave in code
        }
        if(p_value < 0.05){
          comp <- as.numeric(comp) #to nr (only works from char->num, no factors!)
          comp <- sprintf("%02d",comp) #to txt with 01 so no issue order if >9groups
          
          boxplot_comp <- plot_multiboxplot(samples_matrix_comp_no0)
          
          png(paste(name_project, "_boxplot_Mcomp" , multiple_comparison, "_std_", name_standardWOsymbol, "_SIGNIF.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
          
          ##post-hoc Dunn test 
          #  => to retrieve which group is responsible for sign diff kruskal
          #https://stackoverflow.com/questions/31434166/which-post-hoc-test-should-i-conduct-after-doing-a-kruskal-wallis
          #https://stackoverflow.com/questions/44120904/difference-between-posthoc-kruskal-dunn-test-and-dunn-test
          #https://www.rdocumentation.org/packages/dunn.test/versions/1.3.5/topics/dunn.test
          library(dunn.test)
          invisible(capture.output({ #wo print output
            if(length(sort(unique(comp))) >2){ #if run mcomp with 2 groups, run wo bh correction
              post_test <- dunn.test(x=comp_variable$variable, g= comp_variable$comp, method = 'bh', altp = TRUE, alpha=0.05) 
            }else(post_test <- dunn.test(x=comp_variable$variable, g= comp_variable$comp, alpha=0.05))
            #with Benjamini-Hochberg correction 
            #(why bh: with >>feat, takes into account nr of FDR
            #Bonferroni "punishes" all input p-values equally, whereas Benjamini-Hochberg (as a way to control the FDR) "punishes" p-values accordingly to their ranking. If you have few tests, it should not make too much of a difference, but if you have many tests (as for example when you test all 20,000 human protein coding genes or even larger sets), it will make a difference. In this case Bonferroni will produce false negatives, in other words it will discard significant observations. Therefore usually I use Benjamini-Hochberg.
            #https://www.researchgate.net/post/What_is_your_prefered_p-value_correction_for_multiple_tests)
          }))
          
          #report posthoc p-adjusted values (+ all info) for each multiple comparison and each variable (add here)
          linePH <- as.data.frame(post_test)
          linePH$multicomp <- multiple_comparison
          linePH$standard <- name_standard
          if(length(sort(unique(comp))) == 2){ #if run mcomp with 2 groups, run wo bh correction
            linePH$altP <- "No BH correction performed"
            linePH$altP.adjusted <- "No BH correction performed"
          }
          linePH <- linePH[ , c("multicomp", "standard", "comparisons", "chi2", "Z", "altP", "altP.adjusted")]     #reorer cols
          
          append_dataframe_to_report(linePH, name_report_posthocDunn)
          
          #plot visual Dunn test Significant Differences 
          #=> NO, not based of diff in mean, since not-normal, only eg distrib possible but same as boxplot
          #https://www.researchgate.net/post/Kruskal-wallis-and-Dunns-test-interpretation     #ku hier niet zeggen mean verschillend, n-param
          
        }
        if(p_value >= 0.05){
          comp <- as.numeric(comp) #to nr (only works from char->num, no factors!)
          comp <- sprintf("%02d",comp) #to txt with 01 so no issue order if >9groups
          
          ## boxplots also NS stds plot
          boxplot_comp <- plot_multiboxplot(samples_matrix_comp_no0)
          png(paste(name_project, "_boxplot_Mcomp" , multiple_comparison, "_std_", name_standardWOsymbol, "_NS.png", sep=""), width=7, height=5, units="in", res=150)
          plot(boxplot_comp)
          dev.off()
        }
      }
    }
    
    ### Calc adj p-value
    df <- read.table(name_report_ANOVA, header=TRUE, sep="\t") 
    df <- df[df$Multiple.comparison == multiple_comparison,]
    q <- p.adjust(df$p.value, method = "fdr")
    q_values <- c(q_values, q)
    
  }  
  
  #### After all multiple comparisons done, write adj p-values to report
  report_multiple <- read.table(name_report_ANOVA, header=TRUE, sep="\t")
  report_multiple$Adjusted_p_value <- q_values
  name_report_ANOVA2 <- paste(name_project,'_adjpvalues_MultipleComparisons.txt', sep="")
  write_dataframe_as_txt_file(report_multiple, name_report_ANOVA2)
  
}



###################


print("R pipeline - Part II: univariate analysis - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################

