### New version by Fran. Jan 2021

LRGASP_calculations <- function (NAME, class.file, junc.file, out.dir, functions.dir) {
  # Get functions and spike-ins IDs
  setwd(functions.dir)
  source("LRGASP_functions.R")
  sirv_list=read.table("SIRVs_ids.txt", header = F)$V1
  ercc_list=read.table("ERCC_ids.txt", header = F)$V1
  
  # identify files in directory
  cat("Evaluation script has being run.\nData used for ", NAME, " pipeline are \n", class.file , "\n", junc.file , "\n")
  sqanti_data=read.table(class.file , sep = "\t", as.is = T, header = T)
  sqanti_data.junc=read.table(junc.file, sep = "\t", as.is = T, header = T)
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

  ### separate spike-ins and isoforms
  sirv_data=sqanti_data[grep("SIRV",sqanti_data$chrom),]
  ercc_data=sqanti_data[grep("ERCC",sqanti_data$chrom),]
  sirv_data.junc=sqanti_data.junc[grep("SIRV",sqanti_data.junc$chrom),]
  ercc_data.junc=sqanti_data.junc[grep("ERCC",sqanti_data.junc$chrom),]

  ### remove SIRV and ERCC transcripts from sqanti data
  sqanti_data=sqanti_data[grep("SIRV|ERCC",sqanti_data$chrom, invert=T),]
  sqanti_data.junc=sqanti_data.junc[grep("SIRV|ERCC",sqanti_data.junc$chrom, invert=T),]
  
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
                                         "5' reference supported (gene)", "3' reference supported (gene)", "5' CAGE supported", "3' polyA motif supported",
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
  a.FSM_results["5' CAGE supported","Absolute value"]= FSM_TPR_5primeTP_abs
  a.FSM_results["5' CAGE supported","Relative value (%)"]= round(FSM_TPR_5primeTP, digits = 2)
  a.FSM_results["3' polyA motif supported","Absolute value"]= FSM_TPR_3primeTP_abs
  a.FSM_results["3' polyA motif supported","Relative value (%)"]= round(FSM_TPR_3primeTP, digits = 2)
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
                                         "5' CAGE supported", "3' polyA motif supported",
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
  b.ISM_results["5' CAGE supported","Absolute value"]=ISM_TPR_5primeTP_abs
  b.ISM_results["5' CAGE supported","Relative value (%)"]=round(ISM_TPR_5primeTP, digits = 2)
  b.ISM_results["3' polyA motif supported","Absolute value"]=ISM_TPR_3primeTP_abs
  b.ISM_results["3' polyA motif supported","Relative value (%)"]=round(ISM_TPR_3primeTP, digits = 2)
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
                                         "5' reference supported (gene)", "3' reference supported (gene)", 
                                         "5' CAGE supported", "3' polyA motif supported",
                                         "Supported Novel Transcript Model (SNTM)", "Intron retention incidence"))
  c.NIC_results[,"Absolute value"]="-"
  c.NIC_results[,"Relative value (%)"]="-"
  c.NIC_results["Number of isoforms","Absolute value"]=as.integer(dim(sqanti_data_NIC)[1])
  c.NIC_results["5' and 3' reference supported (gene)","Absolute value"]=NIC_TPR_TP_gene_abs
  c.NIC_results["5' and 3' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_gene, digits = 2)
  c.NIC_results["5' reference supported (gene)","Absolute value"]=NIC_TPR_TP_ref5_gene_abs
  c.NIC_results["5' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_ref5_gene, digits = 2)
  c.NIC_results["3' reference supported (gene)","Absolute value"]=NIC_TPR_TP_ref3_gene_abs
  c.NIC_results["3' reference supported (gene)","Relative value (%)"]=round(NIC_TPR_TP_ref3_gene, digits = 2)
  c.NIC_results["5' CAGE supported","Absolute value"]=NIC_TPR_5primeTP_abs
  c.NIC_results["5' CAGE supported","Relative value (%)"]=round(NIC_TPR_5primeTP, digits = 2)
  c.NIC_results["3' polyA motif supported","Absolute value"]=NIC_TPR_3primeTP_abs
  c.NIC_results["3' polyA motif supported","Relative value (%)"]=round(NIC_TPR_3primeTP, digits = 2)
  c.NIC_results["Supported Novel Transcript Model (SNTM)", "Absolute value" ]=NIC_TPR_allTP_abs
  c.NIC_results["Supported Novel Transcript Model (SNTM)", "Relative value (%)" ]=round(NIC_TPR_allTP, digits = 2)
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
                                         "5' reference supported (gene)", "3' reference supported (gene)", 
                                         "5' CAGE supported", "3' polyA motif supported",
                                         "Supported Novel Transcript Model (SNTM)", "Non-canonical SJ incidence",
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
  d.NNC_results["5' CAGE supported","Absolute value"]=NNC_TPR_5primeTP_abs
  d.NNC_results["5' CAGE supported","Relative value (%)"]=round(NNC_TPR_5primeTP, digits = 2)
  d.NNC_results["3' polyA motif supported","Absolute value"]=NNC_TPR_3primeTP_abs
  d.NNC_results["3' polyA motif supported","Relative value (%)"]=round(NNC_TPR_3primeTP, digits = 2)
  d.NNC_results["Supported Novel Transcript Model (SNTM)", "Absolute value"]=NNC_TPR_allTP_abs
  d.NNC_results["Supported Novel Transcript Model (SNTM)", "Relative value (%)"]=round(NNC_TPR_allTP, digits = 2)
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
  
  ### Evaluation of SIRVs
  ##############################################
  print ("SIRVs evaluation")
  sirv_data$TP=apply(sirv_data,1,TP_function)
  SIRVs_transcripts=as.integer(length(sirv_data$isoform))
  SIRVs_called=intersect(sirv_data[which(sirv_data$structural_category=="FSM" & 
                                           sirv_data$TP==TRUE),"associated_transcript"],
                         sirv_list)
  TP=length(SIRVs_called)
  RM_isoforms=sirv_data[which(sirv_data$structural_category=="FSM" & 
                 sirv_data$TP==TRUE),"isoform"]
  RM=length(RM_isoforms)
  
  SIRVs_transcripts_incomplete=sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM") & 
                                                 sirv_data$TP==FALSE),"isoform"]
  
  SIRVs_called_wrong_ends=setdiff(intersect(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"],
                                    sirv_list),SIRVs_called)
  PTP=length(SIRVs_called_wrong_ends)
  SIRVs_not_detected=setdiff(sirv_list,sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"])
  FN=length(SIRVs_not_detected)
  FP_sirvs_detected=sirv_data[-which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"]
  FP=length(FP_sirvs_detected)
  
  SIRVs_redundancy=length(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"isoform"])/length(unique(sirv_data[which(sirv_data$structural_category %in% c("FSM","ISM")),"associated_transcript"]))
  

  # Write out results
  e.SIRVs_results=data.frame(row.names = c("SIRV transcripts", "True Positive detections (TP)", "SIRV transcripts associated to TP (Reference Match)",
                                           "Partial True Positive detections (PTP)", "SIRV transcripts associated to PTP",
                                           "False Negative (FN)", "False Positive (FP)", 
                                           "Sensitivity", "Precision",
                                           "Non Redundant Precision","Positive Detection Rate",
                                           "False Discovery Rate","Redundancy"))
  e.SIRVs_results[,"Value"]="-"
  e.SIRVs_results["SIRV transcripts","Value"]=SIRVs_transcripts
  e.SIRVs_results["True Positive detections (TP)","Value"]=as.integer(TP)
  e.SIRVs_results["SIRV transcripts associated to TP (Reference Match)","Value"]=as.integer(RM)
  e.SIRVs_results["Partial True Positive detections (PTP)","Value"]=as.integer(PTP)
  e.SIRVs_results["SIRV transcripts associated to PTP","Value"]=as.integer(length(SIRVs_transcripts_incomplete))
  e.SIRVs_results["False Negative (FN)","Value"]=as.integer(FN)
  e.SIRVs_results["False Positive (FP)","Value"]=as.integer(FP)
  e.SIRVs_results["Sensitivity","Value"]=round(TP/length(sirv_list), digits = 2)
  e.SIRVs_results["Precision","Value"]=round(RM/SIRVs_transcripts, digits = 2)
  e.SIRVs_results["Non Redundant Precision","Value"]=round(TP/SIRVs_transcripts, digits = 2)
  e.SIRVs_results["Positive Detection Rate", "Value"]=round(length(unique(c(SIRVs_called,SIRVs_called_wrong_ends)))/length(sirv_list), digits = 2)
  e.SIRVs_results["False Discovery Rate","Value"]=round((SIRVs_transcripts - RM)/SIRVs_transcripts, digits = 2)
  e.SIRVs_results["Redundancy","Value"]=round(SIRVs_redundancy, digits = 2)
  

  ### GLOBAL INFORMATION
  sqanti_data$novelGene <- "Annotated Genes"
  sqanti_data[grep("novelGene", sqanti_data$associated_gene), "novelGene"] <- "Novel Genes"
  sqanti_data$novelGene = factor(sqanti_data$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  unique_junc=unique(sqanti_data.junc[,c("chrom", "strand", "genomic_start_coord", "genomic_end_coord", "junction_category","canonical")])
  
  genes_df=aggregate(sqanti_data$isoform, by = list("associatedGene" = sqanti_data$associated_gene, "novelGene" = sqanti_data$novelGene), length)
  
  f.global_results=data.frame(row.names = c("Number of genes detected", "Number of known genes detected",
                                            "Number of transcripts detected", "Number of transcripts associated to a known gene",
                                            "Number of unique SJ detected"))
  g.global_results=data.frame(row.names = c("Novel SJ", "Non-canonical SJ"))
  
  perc_novel=nrow(unique_junc[which(unique_junc$junction_category=="novel"),])/nrow(unique_junc)
  perc_non_can=nrow(unique_junc[which(unique_junc$canonical=="non_canonical"),])/nrow(unique_junc)
  
  f.global_results$Value=NA
  f.global_results["Number of genes detected","Value"]=as.integer(nrow(genes_df))
  f.global_results["Number of known genes detected", "Value"]=as.integer(nrow(genes_df[which(genes_df$novelGene=="Annotated Genes"),]))
  f.global_results["Number of transcripts detected","Value"]=as.integer(nrow(sqanti_data))
  f.global_results["Number of transcripts associated to a known gene","Value"]=as.integer(nrow(sqanti_data[which(sqanti_data$structural_category %in% c("FSM","ISM","NIC","NNC")),]))
  f.global_results["Number of unique SJ detected","Value"]=as.integer(nrow(unique_junc))
  
  g.global_results["Novel SJ","Absolute value"]=nrow(unique_junc[which(unique_junc$junction_category=="novel"),])
  g.global_results["Novel SJ","Relative value (%)"]=round(perc_novel, digits = 2)
  g.global_results["Non-canonical SJ","Absolute value"]=nrow(unique_junc[which(unique_junc$canonical=="non_canonical"),])
  g.global_results["Non-canonical SJ","Relative value (%)"]=round(perc_non_can, digits = 2)
  
                              
                  
  
  ####Create a list with all results and save all
  ###############################################
  
  files <- ls(pattern = "_results")
  all.results <- list()
  for ( i in 1: length(files) ) {
    all.results[[i]] <- eval(parse(text = files[i]))
  }
  setwd(out.dir)
  names(all.results) <- c("FSM", "ISM", "NIC", "NNC", "SIRV", "global", "global_SJ") 
  
  save(all.results , file = paste(NAME, "_results.RData", sep = ''))
  save(sqanti_data, file=paste(NAME, "_classification.RData", sep = ''))
  save(sqanti_data.junc, file=paste(NAME, "_junctions.RData", sep = ''))
  
  save(sqanti_data_FSM, file = paste(NAME, "_FSM.RData", sep=''))
  save(sqanti_data_ISM, file = paste(NAME, "_ISM.RData", sep=''))
  save(sqanti_data_NIC, file = paste(NAME, "_NIC.RData", sep=''))
  save(sqanti_data_NNC, file = paste(NAME, "_NNC.RData", sep=''))
  
  save(sirv_data, file=paste(NAME, "_SIRVs_class.RData", sep=''))
  save(sirv_data.junc, file=paste(NAME, "_SIRVs_junc.RData", sep=''))
  
}

