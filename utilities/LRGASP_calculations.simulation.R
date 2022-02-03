### New version by Fran. Jan 2021

LRGASP_calculations <- function (NAME, class.file, junc.file, out.dir, functions.dir, sim_prefix, onlyReport=TRUE) {
  # Get functions and spike-ins IDs
  setwd(functions.dir)
  source("LRGASP_functions.R")
  
  novel_sim_file = paste(sim_prefix, ".novel_isoforms.tsv", sep = '')
  sim_counts_file = paste(sim_prefix, ".counts.tsv", sep = '')

  # NEW
  sim_list = read.table(novel_sim_file, header = F)$V1

  tot.sim_list = read.table(sim_counts_file, header = F, sep ="\t")

  tot.0.sim_list = setdiff(as.character(tot.sim_list[,1]),sim_list)
  tot.1.sim_list = setdiff(as.character(tot.sim_list[which(tot.sim_list[,3] >= 1), 1]),sim_list)
  tot.5.sim_list = setdiff(as.character(tot.sim_list[which(tot.sim_list[,3] >= 5), 1]),sim_list)
  tot.sim_list = as.character(tot.sim_list[,1])
  
  # identify files in directory
  cat("Evaluation script has being run.\nData used for ", NAME, " pipeline are \n", class.file , "\n", junc.file , "\n")
  sqanti_data=read.table(class.file , sep = "\t", as.is = T, header = T)
  sqanti_data.junc=read.table(junc.file, sep = "\t", as.is = T, header = T)
  
  if (onlyReport==FALSE){
  
  if(all(is.na(sqanti_data$iso_exp))){
    sqanti_data$iso_exp <- 0
  }
  # change names of structural categories
  cat_levels <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
  cat_labels <- c("FSM", "ISM", "NIC", "NNC", "Genic_Genomic",  "Antisense", "Fusion","Intergenic", "Genic_Intron")
  sqanti_data$structural_category = factor(sqanti_data$structural_category,
                                           labels = cat_labels,
                                           levels = cat_levels,
                                           ordered=TRUE)
  
  ### Add LRGASP_id
  
  iso_tags <- isoformTags(sqanti_data.junc)
  sqanti_data <- merge(sqanti_data, iso_tags, by="isoform" , all.x=TRUE)
  sqanti_data$LRGASP_id <- apply(sqanti_data, 1, monoexon_tag)
  sqanti_data <- addSC(sqanti_data)
  
  ### rewrite initial table with LRGASP_id
  
  write.table(sqanti_data, class.file, quote=F, sep = "\t", row.names = FALSE)
  
  }
  ### Evaluation of FSM
  #####################
  print ("FSM evaluation")
  sqanti_data_FSM=subset(sqanti_data, structural_category=="FSM") # should we filter out mono-exons?
  
  # FSM with both 3' and 5'end at less than 50bp of the TSS/TTS associated to the reference match
  sqanti_data_FSM$TP=apply(sqanti_data_FSM,1,TP_function)
  FSM_TPR_TP_abs=length(which(sqanti_data_FSM$TP==TRUE))
  FSM_TPR_TP=FSM_TPR_TP_abs*100/dim(sqanti_data_FSM)[1]
  
  # FSM with both 3' and 5'end at less than 50bp of any TSS/TTS annotated for that gene
  sqanti_data_FSM$TP_gene=apply(sqanti_data_FSM,1,TP_gene_function)
  FSM_TPR_TP_gene_abs=length(which(sqanti_data_FSM$TP_gene==TRUE))
  FSM_TPR_TP_gene=FSM_TPR_TP_gene_abs*100/dim(sqanti_data_FSM)[1]
  
  # Analysis 5' and 3' ends
  sqanti_data_FSM$TP_ref5=apply(sqanti_data_FSM,1,ref5TP_function)  # 5' end matches the reference
  sqanti_data_FSM$TP_ref5_gene=apply(sqanti_data_FSM,1,ref5TP_gene_function) #5' end matches any TSS of the same gene in the reference
  sqanti_data_FSM$TP_ref3=apply(sqanti_data_FSM,1,ref3TP_function)  # 3' end matches the reference
  sqanti_data_FSM$TP_ref3_gene=apply(sqanti_data_FSM,1,ref3TP_gene_function) #3' end matches any TTS of the same gene in the reference
  sqanti_data_FSM$TP_5prime=apply(sqanti_data_FSM,1,fiveTP_function) # 5' end matches a CAGE peak
  sqanti_data_FSM$TP_3prime=apply(sqanti_data_FSM,1,threeTP_function) # 3' end has polyA motif
  
  FSM_TPR_ref5_abs=length(which(sqanti_data_FSM$TP_ref5==TRUE))
  FSM_TPR_ref5=FSM_TPR_ref5_abs*100/dim(sqanti_data_FSM)[1]
  FSM_TPR_ref5_gene_abs=length(which(sqanti_data_FSM$TP_ref5_gene==TRUE))
  FSM_TPR_ref5_gene=FSM_TPR_ref5_gene_abs*100/dim(sqanti_data_FSM)[1]
  FSM_TPR_ref3_abs=length(which(sqanti_data_FSM$TP_ref3==TRUE))
  FSM_TPR_ref3=FSM_TPR_ref3_abs*100/dim(sqanti_data_FSM)[1]
  FSM_TPR_ref3_gene_abs=length(which(sqanti_data_FSM$TP_ref3_gene==TRUE))
  FSM_TPR_ref3_gene=FSM_TPR_ref3_gene_abs*100/dim(sqanti_data_FSM)[1]
  
  FSM_TPR_5primeTP_abs=length(which(sqanti_data_FSM$TP_5prime==TRUE))
  FSM_TPR_5primeTP=FSM_TPR_5primeTP_abs*100/dim(sqanti_data_FSM)[1] # rate for 5'end matching CAGE
  FSM_TPR_3primeTP_abs=length(which(sqanti_data_FSM$TP_3prime==TRUE))
  FSM_TPR_3primeTP=FSM_TPR_3primeTP_abs*100/dim(sqanti_data_FSM)[1] # rate for 3'end with polyA
  
  # "All TP" have any support (by reference gene/transcript or CAGE/polyA) at both 5' and 3' end
  sqanti_data_FSM$TP_all=apply(sqanti_data_FSM,1,allTP_function)
  FSM_TPR_allTP_abs=length(which(sqanti_data_FSM$TP_all==TRUE))
  FSM_TPR_allTP=FSM_TPR_allTP_abs*100/dim(sqanti_data_FSM)[1]
  
  ##### NEEDED?
  # Normalized values
  normalized_FSM_TP=sqanti_data_FSM[,c("associated_transcript","TP_ref5","TP_ref3", "TP", "TP_5prime","TP_3prime","polyA_motif","pos_cage_peak")]
  normalized_FSM_TP=unique(normalized_FSM_TP)
  normalized_FSM_TPR_TP=length(which(normalized_FSM_TP$TP==TRUE))*100/dim(normalized_FSM_TP)[1]
  normalized_FSM_TPR_3primeTP=length(which(normalized_FSM_TP$TP_3prime==TRUE))*100/dim(normalized_FSM_TP)[1]
  normalized_FSM_TPR_5primeTP=length(which(normalized_FSM_TP$TP_5prime==TRUE))*100/dim(normalized_FSM_TP)[1]
  normalized_FSM_TP$allTP=apply(normalized_FSM_TP,1,allTP_norm)
  normalized_FSM_TP[which(is.na(normalized_FSM_TP$allTP)),"allTP"]=TRUE
  normalized_FSM_TPR_allTP=length(which(normalized_FSM_TP$allTP==TRUE))*100/dim(normalized_FSM_TP)[1]
  ######
  
  # Redundancy
  FSM_reference_redundancy=dim(sqanti_data_FSM)[1]/length(unique(sqanti_data_FSM$associated_transcript))
  
  # SJ mean coverage
  sqanti_data_FSM$mean_all_coverage=apply(sqanti_data_FSM,1, mean_cov_all,sqanti_data.junc)
  
  # Write out results
  a.FSM_results=data.frame(row.names = c("Number of isoforms","Reference Match", "5' reference supported (transcript)", "3' reference supported (transcript)",
                                         "5' reference supported (gene)", "3' reference supported (gene)",
                                         "Supported Reference Transcript Model (SRTM)", "Reference redundancy Level"))
  a.FSM_results[,"Absolute value"]="-"
  a.FSM_results[,"Relative value (%)"]="-"
  
  a.FSM_results["Number of isoforms","Absolute value"]=as.integer(dim(sqanti_data_FSM)[1])
  a.FSM_results["Reference Match","Absolute value"]= FSM_TPR_TP_abs
  a.FSM_results["Reference Match","Relative value (%)"]= round(FSM_TPR_TP, digits = 2)
  a.FSM_results["5' reference supported (transcript)","Absolute value"]= FSM_TPR_ref5_abs
  a.FSM_results["5' reference supported (transcript)","Relative value (%)"]= round(FSM_TPR_ref5, digits = 2)
  a.FSM_results["3' reference supported (transcript)","Absolute value"]= FSM_TPR_ref3_abs
  a.FSM_results["3' reference supported (transcript)","Relative value (%)"]= round(FSM_TPR_ref3, digits = 2)
  a.FSM_results["5' reference supported (gene)","Absolute value"]= FSM_TPR_ref5_gene_abs
  a.FSM_results["5' reference supported (gene)","Relative value (%)"]= round(FSM_TPR_ref5_gene, digits = 2)
  a.FSM_results["3' reference supported (gene)","Absolute value"]= FSM_TPR_ref3_gene_abs
  a.FSM_results["3' reference supported (gene)","Relative value (%)"]= round(FSM_TPR_ref3_gene, digits = 2)
  a.FSM_results["Supported Reference Transcript Model (SRTM)","Absolute value"]= FSM_TPR_allTP_abs
  a.FSM_results["Supported Reference Transcript Model (SRTM)","Relative value (%)"]= round(FSM_TPR_allTP, digits = 2)
  a.FSM_results["Reference redundancy Level","Absolute value"]= round(FSM_reference_redundancy, digits = 2)
  
  ### Evaluation of ISM
  ######################
  print ("ISM evaluation")
  sqanti_data_ISM=subset(sqanti_data, structural_category=="ISM")
  
  # ISM with both 3' and 5'end at less than 50bp of the TSS/TTS associated to the reference match
  sqanti_data_ISM$TP=apply(sqanti_data_ISM,1,TP_function)
  ISM_TPR_TP_abs=length(which(sqanti_data_ISM$TP==TRUE))
  ISM_TPR_TP=ISM_TPR_TP_abs*100/dim(sqanti_data_ISM)[1]
  
  # ISM with both 3' and 5'end at less than 50bp of any TSS/TTS annotated for that gene
  sqanti_data_ISM$TP_gene=apply(sqanti_data_ISM,1,TP_gene_function)
  ISM_TPR_TP_gene_abs=length(which(sqanti_data_ISM$TP_gene==TRUE))
  ISM_TPR_TP_gene=ISM_TPR_TP_gene_abs*100/dim(sqanti_data_ISM)[1]
  
  # Analysis 5' and 3' ends
  sqanti_data_ISM$TP_ref5=apply(sqanti_data_ISM,1,ref5TP_function)
  sqanti_data_ISM$TP_ref5_gene=apply(sqanti_data_ISM,1,ref5TP_gene_function) #5' end matches any TSS of the same gene in the reference
  sqanti_data_ISM$TP_ref3=apply(sqanti_data_ISM,1,ref3TP_function)
  sqanti_data_ISM$TP_ref3_gene=apply(sqanti_data_ISM,1,ref3TP_gene_function) #3' end matches any TTS of the same gene in the reference
  sqanti_data_ISM$TP_5prime=apply(sqanti_data_ISM,1,fiveTP_function)
  sqanti_data_ISM$TP_3prime=apply(sqanti_data_ISM,1,threeTP_function)
  
  # Calculate missing exons
  sqanti_data_ISM$missing_exons=apply(sqanti_data_ISM , 1, missing_exons_function)
  sqanti_data_ISM$missing_exons_perc=apply(sqanti_data_ISM,1, function(x){ as.integer(x["missing_exons"])/as.integer(x["ref_exons"])})
  
  ISM_TPR_5primeTP_abs=length(which(sqanti_data_ISM$TP_5prime==TRUE))
  ISM_TPR_5primeTP=ISM_TPR_5primeTP_abs*100/dim(sqanti_data_ISM)[1]
  ISM_TPR_3primeTP_abs=length(which(sqanti_data_ISM$TP_3prime==TRUE))
  ISM_TPR_3primeTP=ISM_TPR_3primeTP_abs*100/dim(sqanti_data_ISM)[1]
  ISM_TPR_ref5_abs=length(which(sqanti_data_ISM$TP_ref5==TRUE))
  ISM_TPR_ref5=ISM_TPR_ref5_abs*100/dim(sqanti_data_ISM)[1]
  ISM_TPR_ref5_gene_abs=length(which(sqanti_data_ISM$TP_ref5_gene==TRUE))
  ISM_TPR_ref5_gene=ISM_TPR_ref5_gene_abs*100/dim(sqanti_data_ISM)[1]
  ISM_TPR_ref3_abs=length(which(sqanti_data_ISM$TP_ref3==TRUE))
  ISM_TPR_ref3=ISM_TPR_ref3_abs*100/dim(sqanti_data_ISM)[1]
  ISM_TPR_ref3_gene_abs=length(which(sqanti_data_ISM$TP_ref3_gene==TRUE))
  ISM_TPR_ref3_gene=ISM_TPR_ref3_gene_abs*100/dim(sqanti_data_ISM)[1]
  
  # All TP, have any support at either 5' or 3' end
  sqanti_data_ISM$TP_all=apply(sqanti_data_ISM,1,allTP_function)
  ISM_TPR_allTP_abs=length(which(sqanti_data_ISM$TP_all==TRUE))
  ISM_TPR_allTP=ISM_TPR_allTP_abs*100/dim(sqanti_data_ISM)[1]
  
  # Normalized values ### NEEDED??
  normalized_ISM_TP=sqanti_data_ISM[ , c("associated_transcript","TP_ref5","TP_ref3", "TP","TP_5prime","TP_3prime","polyA_motif","pos_cage_peak")]
  normalized_ISM_TP=unique(normalized_ISM_TP)
  normalized_ISM_TPR_TP=length(which(normalized_ISM_TP$TP==TRUE))*100/dim(normalized_ISM_TP)[1]
  normalized_ISM_TPR_3primeTP=length(which(normalized_ISM_TP$TP_3prime==TRUE))*100/dim(normalized_ISM_TP)[1]
  normalized_ISM_TPR_5primeTP=length(which(normalized_ISM_TP$TP_5prime==TRUE))*100/dim(normalized_ISM_TP)[1]
  normalized_ISM_TP$allTP=apply(normalized_ISM_TP,1,allTP_norm)
  normalized_ISM_TP[which(is.na(normalized_ISM_TP$allTP)),"allTP"]=TRUE
  normalized_ISM_TPR_allTP=length(which(normalized_ISM_TP$allTP==TRUE))*100/dim(normalized_ISM_TP)[1]
  
  # Redundancy ### NEEDED??
  ISM_reference_redundancy=dim(sqanti_data_ISM)[1]/length(unique(sqanti_data_ISM$associated_transcript))
  
  # SJ mean coverage
  sqanti_data_ISM$mean_all_coverage=apply(sqanti_data_ISM,1, mean_cov_all,sqanti_data.junc)
  
  # Write out results
  b.ISM_results=data.frame(row.names = c("Number of isoforms", "5' reference supported (transcript)", "3' reference supported (transcript)",
                                         "5' and 3' reference supported (gene)", "5' reference supported (gene)", "3' reference supported (gene)",
                                         "Supported Reference Transcript Model (SRTM)", "Reference redundancy Level"))
  b.ISM_results[,"Absolute value"]="-"
  b.ISM_results[,"Relative value (%)"]="-"
  b.ISM_results["Number of isoforms","Absolute value"]=as.integer(dim(sqanti_data_ISM)[1])
  b.ISM_results["5' reference supported (transcript)","Absolute value"]=ISM_TPR_ref5_abs
  b.ISM_results["5' reference supported (transcript)","Relative value (%)"]=round(ISM_TPR_ref5, digits = 2)
  b.ISM_results["3' reference supported (transcript)","Absolute value"]=ISM_TPR_ref3_abs
  b.ISM_results["3' reference supported (transcript)","Relative value (%)"]=round(ISM_TPR_ref3, digits = 2)
  b.ISM_results["5' and 3' reference supported (gene)", "Absolute value"]=ISM_TPR_TP_gene_abs
  b.ISM_results["5' and 3' reference supported (gene)","Relative value (%)"]=round(ISM_TPR_TP_gene, digits = 2)
  b.ISM_results["5' reference supported (gene)","Absolute value"]=ISM_TPR_ref5_gene_abs
  b.ISM_results["5' reference supported (gene)","Relative value (%)"]=round(ISM_TPR_ref5_gene_abs, digits = 2)
  b.ISM_results["3' reference supported (gene)","Absolute value"]=ISM_TPR_ref3_gene_abs
  b.ISM_results["3' reference supported (gene)","Relative value (%)"]=round(ISM_TPR_ref3_gene, digits = 2)
  b.ISM_results["Supported Reference Transcript Model (SRTM)","Absolute value"]=ISM_TPR_allTP_abs
  b.ISM_results["Supported Reference Transcript Model (SRTM)","Relative value (%)"]=round(ISM_TPR_allTP, digits = 2)
  b.ISM_results["Reference redundancy Level","Absolute value"]=round(ISM_reference_redundancy, digits = 2)
  
  ### Evaluation of NIC
  ########################
  print ("NIC evaluation")
  sqanti_data_NIC=subset(sqanti_data, structural_category=="NIC")
  if (nrow(sqanti_data_NIC) > 0 ) {
    sqanti_data_NIC$mean_novel_coverage=apply(sqanti_data_NIC,1, mean_cov_novel, sqanti_data.junc)
    sqanti_data_NIC$mean_known_coverage=apply(sqanti_data_NIC,1, mean_cov_known,sqanti_data.junc)
    sqanti_data_NIC$SJ_wo_cov=apply(sqanti_data_NIC,1,SJ_wo_cov,sqanti_data.junc)
    sqanti_data_NIC$novel_SJ=apply(sqanti_data_NIC,1,novel_SJ_isof,sqanti_data.junc)
    sqanti_data_NIC$TP_ref5_gene=apply(sqanti_data_NIC,1, ref5TP_gene_function)
    sqanti_data_NIC$TP_ref3_gene=apply(sqanti_data_NIC,1, ref3TP_gene_function)
    sqanti_data_NIC$TP_gene=apply(sqanti_data_NIC,1,TP_gene_function)
    sqanti_data_NIC$TP_5prime=apply(sqanti_data_NIC,1,fiveTP_function)
    sqanti_data_NIC$TP_3prime=apply(sqanti_data_NIC,1,threeTP_function)
    
    
    sqanti_data_NIC$TP_all=apply(sqanti_data_NIC,1,allTP_function_novel)
    
    subcat_levels=c("combination_of_known_junctions", "combination_of_known_splicesites" , "intron_retention")
    subcat_labels=c("Comb. known SJ", "Comb. known splice sites", "IR")
    
    sqanti_data_NIC$subcategory=factor(sqanti_data_NIC$subcategory, labels=subcat_labels,
                                       levels=subcat_levels, ordered=TRUE)
    
    NIC_TPR_TP_gene_abs=length(which(sqanti_data_NIC$TP_gene==TRUE))
    NIC_TPR_TP_gene=NIC_TPR_TP_gene_abs*100/dim(sqanti_data_NIC)[1]
    NIC_TPR_TP_ref5_gene_abs=length(which(sqanti_data_NIC$TP_ref5_gene==TRUE))
    NIC_TPR_TP_ref5_gene=NIC_TPR_TP_ref5_gene_abs*100/dim(sqanti_data_NIC)[1]
    NIC_TPR_TP_ref3_gene_abs=length(which(sqanti_data_NIC$TP_ref3_gene==TRUE))
    NIC_TPR_TP_ref3_gene=NIC_TPR_TP_ref3_gene_abs*100/dim(sqanti_data_NIC)[1]
    NIC_TPR_5primeTP_abs=length(which(sqanti_data_NIC$TP_5prime==TRUE))
    NIC_TPR_5primeTP=NIC_TPR_5primeTP_abs*100/dim(sqanti_data_NIC)[1]
    NIC_TPR_3primeTP_abs=length(which(sqanti_data_NIC$TP_3prime==TRUE))
    NIC_TPR_3primeTP=NIC_TPR_3primeTP_abs*100/dim(sqanti_data_NIC)[1]
    NIC_TPR_allTP_abs=length(which(sqanti_data_NIC$TP_all==TRUE))
    NIC_TPR_allTP=NIC_TPR_allTP_abs*100/dim(sqanti_data_NIC)[1]
    NIC_IR_incidence_abs=length(which(sqanti_data_NIC$subcategory=="IR"))
    NIC_IR_incidence=NIC_IR_incidence_abs*100/dim(sqanti_data_NIC)[1] 
    
    ## Write results
    c.NIC_results=data.frame(row.names = c("Number of isoforms", "5' and 3' reference supported (gene)", 
                                           "5' reference supported (gene)", "3' reference supported (gene)","Intron retention incidence"))
    c.NIC_results[,"Absolute value"]="-"
    c.NIC_results[,"Relative value (%)"]="-"
    c.NIC_results["Number of isoforms","Absolute value"]=as.integer(dim(sqanti_data_NIC)[1])
    c.NIC_results["5' and 3' reference supported (gene)","Absolute value"]=NIC_TPR_TP_gene_abs
    c.NIC_results["5' and 3' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_gene, digits = 2)
    c.NIC_results["5' reference supported (gene)","Absolute value"]=NIC_TPR_TP_ref5_gene_abs
    c.NIC_results["5' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_ref5_gene, digits = 2)
    c.NIC_results["3' reference supported (gene)","Absolute value"]=NIC_TPR_TP_ref3_gene_abs
    c.NIC_results["3' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_ref3_gene, digits = 2)
    c.NIC_results["Intron retention incidence","Absolute value"]=NIC_IR_incidence_abs
    c.NIC_results["Intron retention incidence","Relative value (%)"]=round(NIC_IR_incidence, digits = 2)
  } else {
    c.NIC_results=data.frame(row.names = c("Number of isoforms", "5' and 3' reference supported (gene)",
                                           "5' reference supported (gene)", "3' reference supported (gene)",
                                           "5' CAGE supported", "3' polyA motif supported",
                                           "Supported Novel Transcript Model (SNTM)", "Intron retention incidence"))
    c.NIC_results[,"Absolute value"]="0"
    c.NIC_results[,"Relative value (%)"]="0"
    
  }
  
  
  ### Evaluation of NNC
  ########################
  print ("NNC evaluation")
  sqanti_data_NNC=subset(sqanti_data, structural_category=="NNC")
  if (nrow(sqanti_data_NNC) > 0) {
    sqanti_data_NNC$mean_novel_coverage=apply(sqanti_data_NNC,1, mean_cov_novel,sqanti_data.junc)
    sqanti_data_NNC$mean_known_coverage=apply(sqanti_data_NNC,1, mean_cov_known,sqanti_data.junc)
    sqanti_data_NNC$SJ_wo_cov=apply(sqanti_data_NNC,1,SJ_wo_cov,sqanti_data.junc)
    sqanti_data_NNC$Illumina_SJ_support=apply(sqanti_data_NNC,1,SJ_w_cov_perc)
    sqanti_data_NNC$novel_SJ=apply(sqanti_data_NNC,1,novel_SJ_isof,sqanti_data.junc)
    sqanti_data_NNC$novel_SJ_perc=apply(sqanti_data_NNC,1,novel_SJ_isof_perc,sqanti_data.junc)
    sqanti_data_NNC$TP_ref5_gene=apply(sqanti_data_NNC,1, ref5TP_gene_function)
    sqanti_data_NNC$TP_ref3_gene=apply(sqanti_data_NNC,1, ref3TP_gene_function)
    sqanti_data_NNC$TP_gene=apply(sqanti_data_NNC,1,TP_gene_function)
    sqanti_data_NNC$TP_5prime=apply(sqanti_data_NNC,1,fiveTP_function)
    sqanti_data_NNC$TP_3prime=apply(sqanti_data_NNC,1,threeTP_function)
    sqanti_data_NNC$TP_all=apply(sqanti_data_NNC,1,allTP_function_novel)
    
    sqanti_data_NNC$SJ_non_canonical=apply(sqanti_data_NNC,1,non_canonical_SJ,sqanti_data.junc)
    
    subcat_levels=c("at_least_one_novel_splicesite", "intron_retention")
    subcat_labels=c("At least 1 novel SJ", "IR")
    
    sqanti_data_NNC$subcategory=factor(sqanti_data_NNC$subcategory,
                                       labels=subcat_labels,
                                       levels=subcat_levels, 
                                       ordered=TRUE)
    
    NNC_TPR_TP_gene_abs=length(which(sqanti_data_NNC$TP_gene==TRUE))
    NNC_TPR_TP_gene=NNC_TPR_TP_gene_abs*100/dim(sqanti_data_NNC)[1]
    NNC_TPR_TP_ref5_gene_abs=length(which(sqanti_data_NNC$TP_ref5_gene==TRUE))
    NNC_TPR_TP_ref5_gene=NNC_TPR_TP_ref5_gene_abs*100/dim(sqanti_data_NNC)[1]
    NNC_TPR_TP_ref3_gene_abs=length(which(sqanti_data_NNC$TP_ref3_gene==TRUE))
    NNC_TPR_TP_ref3_gene=NNC_TPR_TP_ref3_gene_abs*100/dim(sqanti_data_NNC)[1]
    NNC_TPR_5primeTP_abs=length(which(sqanti_data_NNC$TP_5prime==TRUE))
    NNC_TPR_5primeTP=NNC_TPR_5primeTP_abs*100/dim(sqanti_data_NNC)[1]
    NNC_TPR_3primeTP_abs=length(which(sqanti_data_NNC$TP_3prime==TRUE))
    NNC_TPR_3primeTP=NNC_TPR_3primeTP_abs*100/dim(sqanti_data_NNC)[1]
    NNC_TPR_allTP_abs=length(which(sqanti_data_NNC$TP_all==TRUE))
    NNC_TPR_allTP=NNC_TPR_allTP_abs*100/dim(sqanti_data_NNC)[1]
    NNC_full_Illumina_SJ_support_abs=length(which(sqanti_data_NNC$Illumina_SJ_support==100))
    NNC_full_Illumina_SJ_support=NNC_full_Illumina_SJ_support_abs*100/dim(sqanti_data_NNC)[1]
    NNC_non_canonical_incidence_abs=length(which(sqanti_data_NNC$SJ_non_canonical>0))
    NNC_non_canonical_incidence=NNC_non_canonical_incidence_abs*100/dim(sqanti_data_NNC)[1]
    NNC_RT_switching_incidence_abs=length(which(sqanti_data_NNC$RTS_stage==TRUE))
    NNC_RT_switching_incidence=NNC_RT_switching_incidence_abs*100/dim(sqanti_data_NNC)[1]
    
    
    ## Write results
    d.NNC_results=data.frame(row.names = c("Number of isoforms", "5' and 3' reference supported (gene)", 
                                           "5' reference supported (gene)", "3' reference supported (gene)", "Non-canonical SJ incidence",
                                           "Full Illumina SJ support", "RT-switching incidence"))
    d.NNC_results[,"Absolute value"]="-"
    d.NNC_results[,"Relative value (%)"]="-"
    d.NNC_results["Number of isoforms","Absolute value"]=as.integer(dim(sqanti_data_NNC)[1])
    d.NNC_results["5' and 3' reference supported (gene)","Absolute value"]=NNC_TPR_TP_gene_abs
    d.NNC_results["5' and 3' reference supported (gene)","Relative value (%)"]=round(NNC_TPR_TP_gene, digits = 2)
    d.NNC_results["5' reference supported (gene)","Absolute value"]=NNC_TPR_TP_ref5_gene_abs
    d.NNC_results["5' reference supported (gene)","Relative value (%)"]=round(NNC_TPR_TP_ref5_gene, digits = 2)
    d.NNC_results["3' reference supported (gene)","Absolute value"]=NNC_TPR_TP_ref3_gene_abs
    d.NNC_results["3' reference supported (gene)","Relative value (%)"]=round(NNC_TPR_TP_ref3_gene, digits = 2)
    d.NNC_results["Full Illumina SJ support","Absolute value"]=NNC_full_Illumina_SJ_support_abs
    d.NNC_results["Full Illumina SJ support","Relative value (%)"]=round(NNC_full_Illumina_SJ_support, digits = 2)
    d.NNC_results["Non-canonical SJ incidence","Absolute value"]=NNC_non_canonical_incidence_abs
    d.NNC_results["Non-canonical SJ incidence","Relative value (%)"]=round(NNC_non_canonical_incidence, digits = 2)
    d.NNC_results["RT-switching incidence","Absolute value"]=NNC_RT_switching_incidence_abs
    d.NNC_results["RT-switching incidence","Relative value (%)"]=round(NNC_RT_switching_incidence, digits = 2)
  } else {
    d.NNC_results=data.frame(row.names = c("Number of isoforms", "5' and 3' reference supported (gene)",
                                           "5' reference supported (gene)", "3' reference supported (gene)",
                                           "5' CAGE supported", "3' polyA motif supported",
                                           "Supported Novel Transcript Model (SNTM)", "Non-canonical SJ incidence",
                                           "Full Illumina SJ support", "RT-switching incidence"))
    
    d.NNC_results[,"Absolute value"]="0"
    d.NNC_results[,"Relative value (%)"]="0"
  }
  
  ### Evaluation of simulation
  ##############################################
  
  ### Evaluation of tot.sim
  #####################
  print ("Evaluation using all simulated transcripts")
  num_simulated = length(tot.sim_list)
  sqanti_data$TP=apply(sqanti_data,1,TP_function)
  tot.sim_transcripts=as.integer(length(sqanti_data$isoform))
  tot.sim_called=sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list & 
                                       sqanti_data$TP==TRUE),"associated_transcript"] %>% unique()
  TP=length(tot.sim_called)
  RM_isoforms=sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list & 
                                  sqanti_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  tot.sim_transcripts_incomplete=sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list & 
                                                       sqanti_data$TP==FALSE),"isoform"]
  
  tot.sim_called_wrong_ends=sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list &
                                                  sqanti_data$TP==FALSE), "associated_transcript"] %>%  unique()
  PTP=length(tot.sim_called_wrong_ends)
  tot.sim_not_detected=setdiff(tot.sim_list,sqanti_data$associated_transcript)
  FN=length(tot.sim_not_detected)
  FP_tot.sim_detected=sqanti_data[-which(sqanti_data$associated_transcript %in% tot.sim_list),"isoform"]
  FP=length(FP_tot.sim_detected)
  
  tot.sim_redundancy=length(sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list),"isoform"])/length(unique(sqanti_data[which(sqanti_data$associated_transcript %in% tot.sim_list),"associated_transcript"]))
  
  
  # Write out results
  da.tot.sim_results=data.frame(row.names = c("Number of isoforms simulated", "True Positive detections (TP)", "Number of transcripts associated to TP (Reference Match)",
                                               "Partial True Positive detections (PTP)", "Number of transcripts associated to PTP",
                                               "False Negative (FN)", "False Positive (FP)", 
                                               "Sensitivity", "Precision",
                                               "Non Redundant Precision","Positive Detection Rate",
                                               "False Discovery Rate", "False Detection Rate", "Redundancy"))
  da.tot.sim_results[,"Value"]="-"
  da.tot.sim_results["Number of isoforms simulated","Value"]=num_simulated
  da.tot.sim_results["True Positive detections (TP)","Value"]=as.integer(TP)
  da.tot.sim_results["Number of transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  da.tot.sim_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  da.tot.sim_results["Number of transcripts associated to PTP","Value"]=as.integer(length(tot.sim_transcripts_incomplete))
  da.tot.sim_results["False Negative (FN)","Value"]=as.integer(FN)
  da.tot.sim_results["False Positive (FP)","Value"]=as.integer(FP)
  da.tot.sim_results["Sensitivity","Value"]=round(TP/length(tot.sim_list), digits = 2)
  da.tot.sim_results["Precision","Value"]=round(RM/tot.sim_transcripts, digits = 2)
  da.tot.sim_results["Non Redundant Precision","Value"]=round(TP/tot.sim_transcripts, digits = 2)
  da.tot.sim_results["Positive Detection Rate", "Value"]=round(length(unique(c(tot.sim_called,tot.sim_called_wrong_ends)))/length(tot.sim_list), digits = 2)
  da.tot.sim_results["False Discovery Rate","Value"]=round((FP + PTP)/tot.sim_transcripts, digits = 2)
  da.tot.sim_results["False Detection Rate","Value"]=round((FP)/tot.sim_transcripts, digits = 2)
  da.tot.sim_results["Redundancy","Value"]=round(tot.sim_redundancy, digits = 2)
  
  ### Evaluation of tot.0.sim
  #####################
  print ("Evaluation using all GENCODE simulated transcripts")
  num_simulated = length(tot.0.sim_list)
  tot.0.sim_transcripts=as.integer(length(sqanti_data$isoform))
  tot.0.sim_called=sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list & 
                                                    sqanti_data$TP==TRUE),"associated_transcript"] %>% unique()
  TP=length(tot.0.sim_called)
  RM_isoforms=sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list & 
                                     sqanti_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  tot.0.sim_transcripts_incomplete=sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list & 
                                                          sqanti_data$TP==FALSE),"isoform"]
  
  tot.0.sim_called_wrong_ends=sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list &
                                                  sqanti_data$TP==FALSE), "associated_transcript"] %>%  unique()
  PTP=length(tot.0.sim_called_wrong_ends)
  tot.0.sim_not_detected=setdiff(tot.0.sim_list,sqanti_data$associated_transcript)
  FN=length(tot.0.sim_not_detected)
  FP_tot.0.sim_detected=sqanti_data[-which(sqanti_data$associated_transcript %in% tot.0.sim_list),"isoform"]
  FP=length(FP_tot.0.sim_detected)
  
  tot.0.sim_redundancy=length(sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list),"isoform"])/length(unique(sqanti_data[which(sqanti_data$associated_transcript %in% tot.0.sim_list),"associated_transcript"]))
  
  
  # Write out results
  e.tot.0.sim_results=data.frame(row.names = c("Number of isoforms simulated", "True Positive detections (TP)", "Number of transcripts associated to TP (Reference Match)",
                                               "Partial True Positive detections (PTP)", "Number of transcripts associated to PTP",
                                               "False Negative (FN)", "False Positive (FP)", 
                                               "Sensitivity", "Positive Detection Rate",
                                               "Redundancy"))
  e.tot.0.sim_results[,"Value"]="-"
  e.tot.0.sim_results["Number of isoforms simulated","Value"]=num_simulated
  e.tot.0.sim_results["True Positive detections (TP)","Value"]=as.integer(TP)
  e.tot.0.sim_results["Number of transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  e.tot.0.sim_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  e.tot.0.sim_results["Number of transcripts associated to PTP","Value"]=as.integer(length(tot.0.sim_transcripts_incomplete))
  e.tot.0.sim_results["False Negative (FN)","Value"]=as.integer(FN)
  e.tot.0.sim_results["False Positive (FP)","Value"]=as.integer(FP)
  e.tot.0.sim_results["Sensitivity","Value"]=round(TP/length(tot.0.sim_list), digits = 2)
  e.tot.0.sim_results["Positive Detection Rate", "Value"]=round(length(unique(c(tot.0.sim_called,tot.0.sim_called_wrong_ends)))/length(tot.0.sim_list), digits = 2)
  e.tot.0.sim_results["Redundancy","Value"]=round(tot.0.sim_redundancy, digits = 2)
  
  
  ### Evaluation of tot.1.sim
  #####################
  print ("Evaluation using only GENCODE simulated transcripts with >= 1 TPM")
  num_simulated = length(tot.1.sim_list)

  tot.1.sim_transcripts=as.integer(length(sqanti_data$isoform))
  tot.1.sim_called=sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list & 
                                       sqanti_data$TP==TRUE),"associated_transcript"] %>% unique()
  TP=length(tot.1.sim_called)
  RM_isoforms=sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list & 
                                  sqanti_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  tot.1.sim_transcripts_incomplete=sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list & 
                                                       sqanti_data$TP==FALSE),"isoform"]
  
  tot.1.sim_called_wrong_ends=sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list &
                                                  sqanti_data$TP==FALSE), "associated_transcript"] %>%  unique()
  PTP=length(tot.1.sim_called_wrong_ends)
  tot.1.sim_not_detected=setdiff(tot.1.sim_list,sqanti_data$associated_transcript)
  FN=length(tot.1.sim_not_detected)
  FP_tot.1.sim_detected=sqanti_data[-which(sqanti_data$associated_transcript %in% tot.1.sim_list),"isoform"]
  FP=length(FP_tot.1.sim_detected)
  
  tot.1.sim_redundancy=length(sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list),"isoform"])/length(unique(sqanti_data[which(sqanti_data$associated_transcript %in% tot.1.sim_list),"associated_transcript"]))
  
  
  # Write out results
  f.tot.1.sim_results=data.frame(row.names = c("Number of isoforms simulated", "True Positive detections (TP)", "Number of transcripts associated to TP (Reference Match)",
                                               "Partial True Positive detections (PTP)", "Number of transcripts associated to PTP",
                                               "False Negative (FN)", "False Positive (FP)", 
                                               "Sensitivity", "Positive Detection Rate",
                                               "Redundancy"))
  f.tot.1.sim_results[,"Value"]="-"
  f.tot.1.sim_results["Number of isoforms simulated","Value"]=num_simulated
  f.tot.1.sim_results["True Positive detections (TP)","Value"]=as.integer(TP)
  f.tot.1.sim_results["Number of transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  f.tot.1.sim_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  f.tot.1.sim_results["Number of transcripts associated to PTP","Value"]=as.integer(length(tot.1.sim_transcripts_incomplete))
  f.tot.1.sim_results["False Negative (FN)","Value"]=as.integer(FN)
  f.tot.1.sim_results["False Positive (FP)","Value"]=as.integer(FP)
  f.tot.1.sim_results["Sensitivity","Value"]=round(TP/length(tot.1.sim_list), digits = 2)
  f.tot.1.sim_results["Positive Detection Rate", "Value"]=round(length(unique(c(tot.1.sim_called,tot.1.sim_called_wrong_ends)))/length(tot.1.sim_list), digits = 2)
  f.tot.1.sim_results["Redundancy","Value"]=round(tot.1.sim_redundancy, digits = 2)

  ### Evaluation of tot.5.sim
  #####################
  print ("Evaluation using only GENCODE simulated transcripts with >= 5 TPM")
  num_simulated = length(tot.5.sim_list)
  
  tot.5.sim_transcripts=as.integer(length(sqanti_data$isoform))
  tot.5.sim_called=sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list & 
                                       sqanti_data$TP==TRUE),"associated_transcript"] %>% unique()
  TP=length(tot.5.sim_called)
  RM_isoforms=sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list & 
                                  sqanti_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  tot.5.sim_transcripts_incomplete=sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list & 
                                                       sqanti_data$TP==FALSE),"isoform"]
  
  tot.5.sim_called_wrong_ends=sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list &
                                                  sqanti_data$TP==FALSE), "associated_transcript"] %>%  unique()
  PTP=length(tot.5.sim_called_wrong_ends)
  tot.5.sim_not_detected=setdiff(tot.5.sim_list,sqanti_data$associated_transcript)
  FN=length(tot.5.sim_not_detected)
  FP_tot.5.sim_detected=sqanti_data[-which(sqanti_data$associated_transcript %in% tot.5.sim_list),"isoform"]
  FP=length(FP_tot.5.sim_detected)
  
  tot.5.sim_redundancy=length(sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list),"isoform"])/length(unique(sqanti_data[which(sqanti_data$associated_transcript %in% tot.5.sim_list),"associated_transcript"]))
  
  
  # Write out results
  g.tot.5.sim_results=data.frame(row.names = c("Number of isoforms simulated", "True Positive detections (TP)", "Number of transcripts associated to TP (Reference Match)",
                                               "Partial True Positive detections (PTP)", "Number of transcripts associated to PTP",
                                               "False Negative (FN)", "False Positive (FP)", 
                                               "Sensitivity", "Positive Detection Rate",
                                               "Redundancy"))
  g.tot.5.sim_results[,"Value"]="-"
  g.tot.5.sim_results["Number of isoforms simulated","Value"]=num_simulated
  g.tot.5.sim_results["True Positive detections (TP)","Value"]=as.integer(TP)
  g.tot.5.sim_results["Number of transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  g.tot.5.sim_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  g.tot.5.sim_results["Number of transcripts associated to PTP","Value"]=as.integer(length(tot.5.sim_transcripts_incomplete))
  g.tot.5.sim_results["False Negative (FN)","Value"]=as.integer(FN)
  g.tot.5.sim_results["False Positive (FP)","Value"]=as.integer(FP)
  g.tot.5.sim_results["Sensitivity","Value"]=round(TP/length(tot.5.sim_list), digits = 2)
  g.tot.5.sim_results["Positive Detection Rate", "Value"]=round(length(unique(c(tot.5.sim_called,tot.5.sim_called_wrong_ends)))/length(tot.5.sim_list), digits = 2)
  g.tot.5.sim_results["Redundancy","Value"]=round(tot.5.sim_redundancy, digits = 2)
  
  
  ### Evaluation of sim
  #####################
  print ("Evaluation using only simulated novelty")
  num_simulated = length(sim_list)
  
  sim_transcripts=as.integer(length(sqanti_data$isoform))
  sim_called=sqanti_data[which(sqanti_data$associated_transcript %in% sim_list & 
                                       sqanti_data$TP==TRUE),"associated_transcript"] %>% unique()
  TP=length(sim_called)
  RM_isoforms=sqanti_data[which(sqanti_data$associated_transcript %in% sim_list & 
                                  sqanti_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  sim_transcripts_incomplete=sqanti_data[which(sqanti_data$associated_transcript %in% sim_list & 
                                                       sqanti_data$TP==FALSE),"isoform"]
  
  sim_called_wrong_ends=sqanti_data[which(sqanti_data$associated_transcript %in% sim_list &
                                                  sqanti_data$TP==FALSE), "associated_transcript"] %>%  unique()
  PTP=length(sim_called_wrong_ends)
  sim_not_detected=setdiff(sim_list,sqanti_data$associated_transcript)
  FN=length(sim_not_detected)
  FP_sim_detected=sqanti_data[-which(sqanti_data$associated_transcript %in% sim_list),"isoform"]
  FP=length(FP_sim_detected)
  
  sim_redundancy=length(sqanti_data[which(sqanti_data$associated_transcript %in% sim_list),"isoform"])/length(unique(sqanti_data[which(sqanti_data$associated_transcript %in% sim_list),"associated_transcript"]))
  
  
  # Write out results
  h.sim_results=data.frame(row.names = c("Number of isoforms simulated", "True Positive detections (TP)", "Number of transcripts associated to TP (Reference Match)",
                                               "Partial True Positive detections (PTP)", "Number of transcripts associated to PTP",
                                               "False Negative (FN)", "False Positive (FP)", 
                                               "Sensitivity", "Positive Detection Rate",
                                               "Redundancy"))
  h.sim_results[,"Value"]="-"
  h.sim_results["Number of isoforms simulated","Value"]=num_simulated
  h.sim_results["True Positive detections (TP)","Value"]=as.integer(TP)
  h.sim_results["Number of transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  h.sim_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  h.sim_results["Number of transcripts associated to PTP","Value"]=as.integer(length(sim_transcripts_incomplete))
  h.sim_results["False Negative (FN)","Value"]=as.integer(FN)
  h.sim_results["False Positive (FP)","Value"]=as.integer(FP)
  h.sim_results["Sensitivity","Value"]=round(TP/length(sim_list), digits = 2)
  h.sim_results["Positive Detection Rate", "Value"]=round(length(unique(c(sim_called,sim_called_wrong_ends)))/length(sim_list), digits = 2)
  h.sim_results["Redundancy","Value"]=round(sim_redundancy, digits = 2)

  
  
  ### GLOBAL INFORMATION
  sqanti_data$novelGene <- "Annotated Genes"
  sqanti_data[grep("novelGene", sqanti_data$associated_gene), "novelGene"] <- "Novel Genes"
  sqanti_data$novelGene = factor(sqanti_data$novelGene,
                                 levels = c("Novel Genes","Annotated Genes"),
                                 ordered=TRUE)
  
  unique_junc=unique(sqanti_data.junc[,c("chrom", "strand", "genomic_start_coord", "genomic_end_coord", "junction_category","canonical")])
  
  genes_df=aggregate(sqanti_data$isoform, by = list("associatedGene" = sqanti_data$associated_gene, "novelGene" = sqanti_data$novelGene), length)
  
  i.global_results=data.frame(row.names = c("Number of genes detected", "Number of known genes detected",
                                            "Number of transcripts detected", "Number of transcripts associated to a known gene",
                                            "Number of unique SJ detected"))
  j.global_results=data.frame(row.names = c("Novel SJ", "Non-canonical SJ"))
  
  perc_novel=nrow(unique_junc[which(unique_junc$junction_category=="novel"),])/nrow(unique_junc)
  perc_non_can=nrow(unique_junc[which(unique_junc$canonical=="non_canonical"),])/nrow(unique_junc)
  
  i.global_results$Value=NA
  i.global_results["Number of genes detected","Value"]=as.integer(nrow(genes_df))
  i.global_results["Number of known genes detected", "Value"]=as.integer(nrow(genes_df[which(genes_df$novelGene=="Annotated Genes"),]))
  i.global_results["Number of transcripts detected","Value"]=as.integer(nrow(sqanti_data))
  i.global_results["Number of transcripts associated to a known gene","Value"]=as.integer(nrow(sqanti_data[which(sqanti_data$structural_category %in% c("FSM","ISM","NIC","NNC")),]))
  i.global_results["Number of unique SJ detected","Value"]=as.integer(nrow(unique_junc))
  
  j.global_results["Novel SJ","Absolute value"]=nrow(unique_junc[which(unique_junc$junction_category=="novel"),])
  j.global_results["Novel SJ","Relative value (%)"]=round(perc_novel, digits = 2)
  j.global_results["Non-canonical SJ","Absolute value"]=nrow(unique_junc[which(unique_junc$canonical=="non_canonical"),])
  j.global_results["Non-canonical SJ","Relative value (%)"]=round(perc_non_can, digits = 2)
  
  
  
  ####Create a list with all results and save all
  ###############################################
  
  files <- ls(pattern = "_results")
  all.results <- list()
  for ( i in 1: length(files) ) {
    all.results[[i]] <- eval(parse(text = files[i]))
  }
  setwd(out.dir)
  names(all.results) <- c("FSM", "ISM", "NIC", "NNC", "tot.sim", "tot.0.sim", "tot.1.sim", "tot.5.sim", "novel.sim", "global", "global_SJ") 
  
  save(all.results , file = paste(NAME, "_results.RData", sep = ''))
  save(sqanti_data, file=paste(NAME, "_classification.RData", sep = ''))
  save(sqanti_data.junc, file=paste(NAME, "_junctions.RData", sep = ''))
  
  save(sqanti_data_FSM, file = paste(NAME, "_FSM.RData", sep=''))
  save(sqanti_data_ISM, file = paste(NAME, "_ISM.RData", sep=''))
  save(sqanti_data_NIC, file = paste(NAME, "_NIC.RData", sep=''))
  save(sqanti_data_NNC, file = paste(NAME, "_NNC.RData", sep=''))

  
}

