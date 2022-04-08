
#### LRGASP Functions
TP_function=function(X){
  abs(as.integer(X["diff_to_TSS"]))<=50 & abs(as.integer(X["diff_to_TTS"]))<=50
}

TP_gene_function=function(X){
  abs(as.integer(X["diff_to_gene_TSS"]))<=50 & abs(as.integer(X["diff_to_gene_TTS"]))<=50
}

ref5TP_function=function(X){
  abs(as.integer(X["diff_to_TSS"]))<=50 
}

ref5TP_gene_function=function(X){
  abs(as.integer(X["diff_to_gene_TSS"]))<=50 
}

ref3TP_function=function(X){
  abs(as.integer(X["diff_to_TTS"]))<=50 
}

ref3TP_gene_function=function(X){
  abs(as.integer(X["diff_to_gene_TTS"]))<=50 
}

fiveTP_function=function(X){
  if(isTRUE(as.logical(X["within_cage_peak"]))){TRUE}else{FALSE}
}

threeTP_function=function(X){
  if(isTRUE(as.logical(X["within_polya_site"])) || !is.na(X["polyA_motif"]) ){TRUE}else{FALSE}
}

quantseqTP_function=function(X){
  if(isTRUE(as.logical(X["within_polya_site"]))){TRUE}else{FALSE}
}

allTP_function=function(X){
  (abs(as.integer(X["diff_to_gene_TSS"]))<=50 || isTRUE(as.logical(X["within_cage_peak"]))) && 
  (abs(as.integer(X["diff_to_gene_TTS"]))<=50 || !is.na(X["polyA_motif"]))
}

allTP_norm=function(X){
  (as.logical(as.logical(X["TP_ref5"]) | as.logical(X["TP_5prime"])) & as.logical(as.logical(X["TP_ref3"]) | as.logical(X["TP_3prime"])))
}

mean_cov_novel=function(X , sqanti_data_junc){
  novel_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==as.character(X["isoform"]) & sqanti_data_junc$junction_category=="novel"),]
  mean(novel_SJ$total_coverage)
}

mean_cov_known=function(X , sqanti_data_junc){
  known_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==as.character(X["isoform"]) & sqanti_data_junc$junction_category=="known"),]
  mean(known_SJ$total_coverage)
}

mean_cov_all=function(X, sqanti_data_junc){
  all_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==X["isoform"]),]
  mean(all_SJ$total_coverage)
}


SJ_wo_cov=function(X, sqanti_data_junc ){
  all_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==X["isoform"]),]
  length(which(all_SJ$total_coverage==0))
}

SJ_w_cov_perc=function(X){
  if (as.integer(X["exons"])!=1){
  (1-(as.integer(X["SJ_wo_cov"])/(as.integer(X["exons"])-1)))*100
  }else{0}
}

novel_SJ_isof=function(X, sqanti_data_junc ){
  all_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==X["isoform"]),]
  as.integer(length(which(all_SJ$junction_category=="novel")))
}

novel_SJ_isof_perc=function(X, sqanti_data_junc ){
  all_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==X["isoform"]),]
  as.integer(length(which(all_SJ$junction_category=="novel")))*100/dim(all_SJ)[1]
}

non_canonical_SJ=function(X, sqanti_data_junc ){
  all_SJ=sqanti_data_junc[which(sqanti_data_junc$isoform==X["isoform"]),]
  as.integer(length(which(all_SJ$canonical=="non_canonical")))
}

allTP_function_novel=function(X){
  return(as.logical((abs(as.integer(X["diff_to_gene_TSS"]))<=50) | isTRUE(as.logical(X["within_cage_peak"]))) &
           (as.logical(abs(as.integer(X["diff_to_gene_TTS"]))<=50) | !is.na(X["polyA_motif"])) &
           (as.logical(as.integer(X["min_cov"])>0)))
}

distancias <- function (CLASS, category, dist) {
  cap2 <- function (list) {
    b <- list()
    for ( i in 1:11 ) {
      a <- list[[i]]
      a[a > 1000] <- 1000
      a[a < -1000] <- -1000
      b[[i]] <- a
    }
    b
  }
  myPalette = c(rep("blue",5), rep("red",6))
  diff_to_TSS_data <- sapply(CLASS, function (x) subset(x, structural_category==category)[,dist])
  type.dist <- colnames(CLASS[[1]][dist])
  diff_to_TSS_data_capped <- cap2(diff_to_TSS_data)
  #myhist <- hist(diff_to_TSS_data_capped[[9]], breaks = 100,  plot = FALSE)
  #plot( myhist$mids, myhist$counts,type = "l",
  #   main = paste(category, "\n", type.dist),
  #   xlab = "Distance",
  #   ylab = "", col = myPalette[1], 
  #   ylim = c(0,5000), xlim = c(-1100,1100)
  #   )
  #for ( j in 1 : 8) {
  #myhist <- hist(diff_to_TSS_data_capped[[j]], breaks = 100,  plot = FALSE)
  #lines(myhist$mids, myhist$counts, col = myPalette[j])
  #}
  d=hist(diff_to_TSS_data_capped[[1]], breaks=100, plot = F)
  plot(d$mids, d$counts,col = myPalette[1],
    type="l", main = paste(category, "\n", type.dist),
       xlab = "Distance",
       ylab = "",
       #ylim = c(0,10000),
    xlim = c(-1100,1100)
       )
  for (j in 2:11){
    d=hist(diff_to_TSS_data_capped[[j]], breaks=100, plot = F)
    lines(d$mids, d$counts, type="l", col = myPalette[j])
  }
  
  dens=density(diff_to_TSS_data_capped[[5]])
  PM <- round(sapply(diff_to_TSS_data, function (x) length(which(x == 0)) / length(x)) * 100,2)
  barplot(PM, names.arg = c(1:11), col =  myPalette,
          ylim = c(0,100), yaxp = c(0, 100, 5) , ylab = "Distance",
          main = paste("% Perfect reference match", "\n", category, "\n",type.dist))
}

missing_exons_function=function(X){
  as.integer(X["ref_exons"]) - as.integer(X["exons"])
}

### intrapriming if > 60% in next 20bp after TTS

is_intrapriming=function(X){
  if (as.integer(X["perc_A_downstream_TTS"])>=60){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

#### LRGASP_id

isoformTags <- function(junctions_file) {
  df <- junctions_file[, c("isoform", "chrom", "strand")] # df with isoforms in *junctions.txt
  dt <- data.table::data.table(df)
  dt <- dt[,coord:=paste0(junctions_file$genomic_start_coord, "_", junctions_file$genomic_end_coord)]
  dt <-
    dt[, list(tagcoord = paste0(coord, collapse = "_")),
       by = c("isoform", "chrom", "strand")]
  df <- as.data.frame(dt)
  df <- df[order(df$isoform),]
  df$LRGASP_id <- paste(df$chrom, df$strand, df$tagcoord, sep = "_")
  tag_df <- df[,c("isoform", "LRGASP_id")]
  return(tag_df)
}

monoexon_tag <- function(iso_classif){
  if (is.na(iso_classif["LRGASP_id"])){
    if (iso_classif["strand"]=='-'){
      init <- as.integer(iso_classif["TTS_genomic_coord"])
      end <- as.integer(iso_classif["TSS_genomic_coord"])
    } else {
      init <- as.integer(iso_classif["TSS_genomic_coord"])
      end <- as.integer(iso_classif["TTS_genomic_coord"])
    }
    init <- round(init, digits = -2)
    end <- round(end, digits = -2)
    tag <- paste0(iso_classif["chrom"], "_", iso_classif["strand"], '_', init, '_', end)
  }else{
    tag <- iso_classif["LRGASP_id"]
  }
  tag
}

addSC <- function(class_file){
  class_file$LRGASP_id <- paste0(class_file$structural_category, "_", class_file$LRGASP_id)
  return(class_file)
}
