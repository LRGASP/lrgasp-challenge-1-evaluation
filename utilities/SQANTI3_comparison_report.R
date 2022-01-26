#######################################
#                                     #
#  SQANTI3 output comparison report   #
#             generation              #
#                                     #
#######################################


# Author: Jorge Martinez Tomas & Alejandro Paniagua
# Updated and modified: Francisco J. Pardo-Palacios
# Last modified: 10/20/2021 by Francisco J. Pardo-Palacios


#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(DT))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(optparse))
suppressMessages(library(rmarkdown))
suppressMessages(library(tidyverse))
suppressMessages(library(ggridges))
suppressMessages(library(reshape))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorConesa))
suppressMessages(library(jaccard))
suppressMessages(library(corrplot))

#######################################
#                                     #
#             FUNCTIONS               #
#                                     #
#######################################

# -------------------- Read data

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}


# -------------------- Tags and basic comparison P/A


uniquetag <- function(class_file) {
  # Aggregate by UJC and calculate many metrics when coords ar given
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      associated_gene=list(unique(associated_gene)),
      FL=as.numeric(sum(FL)),
      l_FL=list(FL),
      exons=unique(exons),
      TSS_genomic_coord=list(TSS_genomic_coord),
      TTS_genomic_coord=list(TTS_genomic_coord),
      length=list(length),
      medianTSS = round(median(TSS_genomic_coord)),
      medianTTS = round(median(TTS_genomic_coord)),
      sdTSS = sd(TSS_genomic_coord),
      sdTTS = sd(TTS_genomic_coord),
      isoform = list(isoform)
    ), by = c("LRGASP_id", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$LRGASP_id, dt.out$structural_category),])
}


uniquetag_simple <- function(class_file) {
  # Aggregate by UJC and calculate some basic metrics
  dt <- data.table::data.table(class_file)
  dt.out <-
    dt[, list(
      associated_gene=list(unique(associated_gene)),
      FL=as.numeric(sum(FL)),
      l_FL=list(FL),
      exons=unique(exons),
      length=list(length),
      isoform = list(isoform)
    ), by = c("LRGASP_id", "structural_category")]
  dt.out <- as.data.frame(dt.out)
  return(dt.out[order(dt.out$LRGASP_id, dt.out$structural_category),])
}


multipleComparison <- function(l_class){
  # Presence and ausence of UJC in the samples
  a <- c(rbind(names(l_class), paste0(names(l_class), "SC")))
  
  l_class %>%
    purrr::map(~ data.frame(col = .$LRGASP_id, .$LRGASP_id,.$structural_category, stringsAsFactors = FALSE)) %>%
    purrr::reduce(full_join, by = "col") %>%
    select(-col) %>%
    setNames(a)
}


# -------------------- Isoform (UJC) analysis

SD_TSS_TTS <- function(l_class){
  # Calculate UJC SD
  a <- c("LRGASP_id", rbind(paste0(names(l_class), "TSS"), paste0(names(l_class), "TTS")))
  TSS_TTS_params <- list()
  for (i in 1:length(l_class)){
    TSS_TTS_params[[i]] <- l_class[[i]][,c("LRGASP_id", "TSS_genomic_coord", "TTS_genomic_coord")]
  }
  TSS_TTS_params <- TSS_TTS_params %>% 
    purrr::reduce(full_join, by="LRGASP_id") %>% 
    setNames(a)
  
  a <- paste0(names(l_class), "TSS")
  b <- paste0(names(l_class), "TTS")
  allTSS <- TSS_TTS_params[, a]
  allTSS[allTSS == "NULL"] <- NA
  allTTS <- TSS_TTS_params[, b]
  allTTS[allTTS == "NULL"] <- NA
  
  TSS_TTS_df <- data.frame(LRGASP_id=TSS_TTS_params$LRGASP_id)
  
  # max and min value
  sapplycolumns <- function(data, func){
    tmp <- list()
    for (i in 1:ncol(data)){
      tmp[[names(data)[i]]] <- sapply(data[,i], func)
    }
    for (i in 1:length(tmp)){
      tmp[[i]][tmp[[i]]=="NULL"] <- NA
    }
    return(as.data.frame(tmp))
  }
  
  minNA <- function(x) ifelse(length(x) > 1, min(x), NA)
  
  
  maxTSS <- sapplycolumns(allTSS, max)
  maxTSS[maxTSS=="-Inf"] <- NA
  
  maxTTS <- sapplycolumns(allTTS, max)
  maxTTS[maxTTS=="-Inf"] <- NA
  
  minTSS <- sapplycolumns(allTSS, minNA)
  minTSS[minTSS=="Inf"] <- NA
  
  minTTS <- sapplycolumns(allTTS, minNA)
  minTTS[minTTS=="Inf"] <- NA
  
  minmaxTSS <- cbind(minTSS, maxTSS)
  minmaxTTS <- cbind(minTTS, maxTTS)
  
  TSS_TTS_df$minmax.SD.TSS <- apply(minmaxTSS, 1, function(x) sd(unlist(x), na.rm = TRUE))
  
  TSS_TTS_df$minmax.SD.TTS <- apply(minmaxTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  
  # Median value
  
  medianTSS <- sapplycolumns(allTSS, median)
  medianTTS <- sapplycolumns(allTTS, median)
  
  
  TSS_TTS_df$median.SD.TSS <- apply(medianTSS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  TSS_TTS_df$mean.median.TSS <- apply(medianTSS, 1, function(x) mean(unlist(x),na.rm = TRUE))
  TSS_TTS_df$median.SD.TTS <- apply(medianTTS, 1, function(x) sd(unlist(x),na.rm = TRUE))
  TSS_TTS_df$mean.median.TTS <- apply(medianTTS, 1, function(x) mean(unlist(x),na.rm = TRUE))
  
  # Max SD
  
  TSS_TTS_df$SD.TSS <- apply(TSS_TTS_df[,c("minmax.SD.TSS", "median.SD.TSS")],1,max)
  TSS_TTS_df$SD.TTS <- apply(TSS_TTS_df[,c("minmax.SD.TTS", "median.SD.TTS")],1,max)
  
  return(TSS_TTS_df[, c("LRGASP_id", "SD.TSS", "SD.TTS", "mean.median.TSS", "mean.median.TTS")])
}

calculate_CPM <- function(class_table){
  total_counts=sum(class_table$FL)
  class_table$CPM=apply(class_table,1, function(x){(as.numeric(x["FL"]) * 10^6)/total_counts})
  return(class_table)
}

iso_analysis <- function(l_class){
  l_class <- purrr::map(l_class,calculate_CPM)
  class_bind <- bind_rows(l_class)
  class_bind <- data.table::data.table(class_bind)

  class_compact <- class_bind[, list(
    exons=unique(exons),
    length=as.numeric(median(unlist(length))),
    FL_cpm=as.numeric(median(CPM)),
    gene=unique(unlist(associated_gene)) %>% paste(collapse="-")
  ), by="LRGASP_id"]
  
  class_compact <- as.data.frame(class_compact)
  
  if (TSS_TTS_coord == TRUE) {
    TSS_TTS_df <- SD_TSS_TTS(l_class)
    df_iso <- list(TSS_TTS_df, class_compact) %>% 
      purrr::reduce(full_join, by="LRGASP_id")
  } else {df_iso <- data.frame(
    LRGASP_id = class_compact$LRGASP_id,
    SD.TSS=NA,
    SD.TTS=NA,
    mean.median.TSS=NA,
    mean.median.TTS=NA,
    exons=class_compact$exons,
    length=class_compact$length,
    iso_exp=class_compact$FL_cpm,
    gene=class_compact$gene
  )}
  
  return(df_iso)
  
}

# ---------------- Jaccard index across pipelines 

get_jaccard_matrix <- function(pa, pipe){
  only_pipelines <- pa[, pipe]
  pip <- colnames(only_pipelines)
  id <- expand.grid(pipe,pipe)
  pairwise_intersection_matrix <- matrix( colSums( only_pipelines[ , id[,1] ] == only_pipelines[ , id[,2] ] & only_pipelines[ , id[,1] ]!=0) , ncol = length(pip) )
  pairwise_union_matrix <-  matrix( colSums( (only_pipelines[ , id[,1] ]==1 | only_pipelines[ , id[,2] ]==1) ),
                                    ncol = length(pip) )
  
  jaccard_index_matrix <- pairwise_intersection_matrix/pairwise_union_matrix
  colnames(jaccard_index_matrix) <- pipe
  rownames(jaccard_index_matrix) <- pipe

  return(jaccard_index_matrix)
}

# ------------------ Sum presence of UJC across pipelines
sum_pipelines <- function(x, pip){
  x <- x[pip]
  x <- as.numeric(x)
  z <- sum(x)
  return(z)
}

# -------------------- Gene analysis

del_novel_genes <- function(df){
  df.out <- df[!grepl("novelGene|SIRV", df$associated_gene), ]
  df.out$associated_gene <- as.character(df.out$associated_gene)
  return(df.out)
}

tags_per_gene <- function(res_class){
  # Delete novel genes and SIRV
  ref_genes <- del_novel_genes(res_class)
  ref_genes <- data.table::data.table(ref_genes)
  
  ref_genes.out <- ref_genes[, list(
    N_UJC=length(unique(LRGASP_id)),
    LRGASP_id=list(LRGASP_id)
  ), by="associated_gene"]
  
  return(as.data.frame(ref_genes.out))
}

jaccard <- function(x, comb){
  x <- x[2:length(x)]
  val <- c()
  for (i in 1:ncol(comb)){
    inter <- length(intersect(x[[comb[1,i]]], x[[comb[2,i]]]))
    uni <- length(union(x[[comb[1,i]]], x[[comb[2,i]]]))
    jac <- inter/uni
    val <- c(val, jac)
  }
  return(median(val))
}

get_jaccard <- function(l_df, l_class){
  a <- c("associated_gene", paste0("LRGASP_id_", names(l_class)))
  l <- list()
  for (i in 2:length(l_df)){
    l[[i-1]] <- l_df[[i]][,c("associated_gene","LRGASP_id")]
  }
  df <- l %>% 
    purrr::reduce(full_join, by="associated_gene") %>% 
    setNames(a)
  
  n <- colnames(df[,2:ncol(df)])
  comb <- combn(n, 2)
  
  jc <- apply(df, 1, function(x) jaccard(x, comb))
  df <- data.frame(associated_gene=df$associated_gene, jaccard=jc)
  return(df)
}

gene_expr <- function(res_class){
  ref_genes <- del_novel_genes(res_class)
  ref_genes <- data.table::data.table(ref_genes)
  
  ref_genes.UJC <- ref_genes[, list(
    associated_gene=unique(associated_gene),
    FL=as.numeric(median(FL))
  ), by="LRGASP_id"]
  
  ref_genes.out <- ref_genes.UJC[, list(
    FL_cpm=sum(FL)
  ), by="associated_gene"]
  
  ref_genes.out <- as.data.frame(ref_genes.out)
  
  tot_exp <- sum(ref_genes.out$FL_cpm)
  ref_genes.out$FL_cpm <- (ref_genes.out$FL_cpm * 10^6) / tot_exp
  
  return(ref_genes.out)
}

gene_length <- function(df_gtf, df_gene){
  cond <- which(df_gtf$feature == "exon" & df_gtf$associated_gene %in% df_gene$associated_gene)
  df_gtf <- df_gtf[cond,]
  df_gtf$diff <- df_gtf$end - df_gtf$start
  df_gtf <- data.table::data.table(df_gtf)
  df_gtf.out <- df_gtf[, list(
    length=sum(diff)
  ), by="associated_gene"]
  return(as.data.frame(df_gtf.out))
}

gene_analysis <- function(l_class){
  l_res_class <- list()
  for (i in 1:length(l_class)){
    l_res_class[[i]] <- del_novel_genes(l_class[[i]])
  }
  class_bind <- bind_rows(l_res_class)
  
  l_gene_df <- list()
  l_gene_df[[1]] <- tags_per_gene(class_bind)
  for (i in 1:length(l_class)){
    l_gene_df[[i+1]] <- tags_per_gene(l_class[[i]])
  }
  
  a <- c("associated_gene", "N_UJC", paste0("N_UJC_", names(l_class)))
  df_gene <- l_gene_df %>% 
    purrr::map(~ data.frame(col = .$associated_gene, .$N_UJC, stringsAsFactors = FALSE)) %>% 
    purrr::reduce(full_join, by="col") %>% 
    setNames(a)
  
  df_jac <- get_jaccard(l_gene_df, l_class)
  
  df_exp <- gene_expr(class_bind)
  l_df_metrics <- list(df_gene, df_jac, df_exp)
  
  
  df_gene <- l_df_metrics %>% 
    purrr::reduce(full_join, by="associated_gene")
  
  return(df_gene)
}


# -------------------- Final comparison function

compareTranscriptomes <- function(l_iso){
  print("Generating unique junction chains (UJC)...")
  n <- names(l_iso)
  l_class <- list()
  for ( i in 1:length(l_iso)) {

    class.filtered <- l_iso[[i]][[1]]
    
    if (TSS_TTS_coord){
      #class.swap <- swapcoord(class.filtered) # swap - strand
      #class.uniquetags <- as.data.frame(uniquetag(class.swap)) # group by tag
      class.out <- as.data.frame(uniquetag(class.filtered)) 
    } else{
      class.out <- uniquetag_simple(class.filtered) # group by tag
    }
    
    l_class[[i]] <- class.out # add to list
  }
  
  print("Performing presence/ausence comparison...")
  names(l_class) <- n # add names
  comptags <- multipleComparison(l_class)
  comptags.SC <- cbind(structural_category =
                         do.call(dplyr::coalesce, comptags[,paste0(n,"SC")]),
                       comptags[,n]
  )
  comptags.out <- cbind( TAGS =
                           do.call(dplyr::coalesce, comptags[,n]),
                         comptags.SC
  )
  
  comptags.PA <- comptags.out
  comptags.PA[,3:ncol(comptags.PA)][!is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 1
  comptags.PA[,3:ncol(comptags.PA)][is.na(comptags.PA[,3:ncol(comptags.PA)])] <- 0
  
  print("Performing UJC analysis...")
  iso.metrics <- iso_analysis(l_class)
  
  print("Analysing associated genes...")
  gene.metrics <- gene_analysis(l_class)
  
  ### Known genes
  known_genes <- list()
  for ( i in 1:length(l_iso)) {
    class.filtered <- del_novel_genes(l_iso[[i]][[1]])
    known_genes[[names(l_iso)[i]]] <- as.character(class.filtered$associated_gene)
  }
  
  
  output <-
    list(
      classifications = l_class,
      comparison = comptags.out,
      comparisonPA = comptags.PA,
      iso_metrics = iso.metrics,
      gene_metrics = gene.metrics,
      genes_detected = known_genes
    )
  
  return(output)
}


# -------------------- Functions for the report generation

iso2url <- function(id){
  # Isoform coords to Genome Browser
  if (id != "NA") {
    id.split <- str_split(id,"_")
    chr <- id.split[[1]][2]
    start <- id.split[[1]][4]
    end <- id.split[[1]][length(id.split[[1]])]
    name <- substr(id, 1, 30)
    url <- paste0(
      "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=",
      chr, "%3A", start, "%2D", end, "&hgsid=1143169919_jAraPbUWtMCdAHfgTHk4sDHQHW7R",
      "'>", name, "...</a>")
    return(url)
  } else { return(NA) }
}



#######################################
#                                     #
#                MAIN                 #
#                                     #
#######################################

# -------------------- Argument parser

option_list <- list(
  make_option(c("-d", "--dir"), type = "character", default = NULL,
              help="directory with input files (classification and junction files)",
              metavar = "DIRIN"),
  make_option(c("-o", "--outdir"), type = "character", default = ".",
              help="Output directory for the report and CSV file [default= %default]",
              metavar = "DIROUT"),
  make_option(c("-n", "--name"), type = "character", default = "comparison_output",
              help="Output name for the HTML report and CSV file (without extension) [default= %default]",
              metavar = "OUTNAME"),
  make_option(c("--lrgasp"), action="store_true",type = "character", default = FALSE,
              help="Use lrgasp metrics",
              metavar = "LRGASP"),
  make_option(c("--code"), type="character", default = "code.txt",
              help="Code to convert pipline ID to Library Platform Data Category",
              metavar = "CODE"),
  make_option(c("--upset"), action="store_true",type = "character", default = FALSE,
              help="Generate UpSet plots if n<31",
              metavar = "UPSET")
)

opt_parser = OptionParser(
  usage = "usage: %prog [-i DIRIN] [-o DIROUT] [-n OUTNAME] [--lrgasp] [--upset] [--code]",
  option_list=option_list
)

opt = parse_args(opt_parser)

directory <- opt$dir
output_directory <- opt$outdir
output_name <- opt$name
lrgasp <- opt$lrgasp
pdf.report.file <- paste0(output_directory,"/",output_name,"_comparison_report.pdf")

if (is.null(directory)) {
  stop("\n\nAt least one argument must be supplied.\nThe -d argument is required (directory containing input files)")
}

code=read.csv(opt$code, header=T)
if (dim(code)[2]==4){
  colnames(code) <- c("pipeline", "Library", "Platform", "Data_Category")
}else if (dim(code)[2]){
  colnames(code) <- c("pipeline", "Library", "Platform", "Data_Category","Lab")
}
code=code[order(code$pipeline),]


# -------------------- Load data

print("Reading input data...")
if (dir.exists(directory)){
  dir_in <- directory
} else {
  dir_in <- paste(getwd(), directory, sep="/")
  if (!dir.exists(dir_in)){
    stop(paste0("\n\nCould not find the input directory (", directory, ").\nPlease enter a valid path"))
  }
}

class_in <-
  list.files(dir_in,
             pattern = "*_classification.txt",
             all.files = FALSE,
             full.names = TRUE)
junct_in <-
  list.files(dir_in,
             pattern = "*_junctions.txt",
             all.files = FALSE,
             full.names = TRUE)

if (length(class_in) != length(junct_in)){
  stop("ERROR: There is a different number of classification and junction files in the directory")
} else if (length(class_in) == 0){
  stop(paste0("ERROR: No classification and junction files were found in the directory: ", dir_in))
}

f_in <- list()
lrgasp.res <- list()

if (lrgasp == TRUE){
  lrgasp.files <- 
    list.files(dir_in,
               pattern = "*_results.RData",
               all.files = FALSE,
               full.names = TRUE)
  if (length(class_in) != length(lrgasp.files)){
    print("ERROR: Issue loading LRGASP files")
    print("Different number of LRGASP files than samples")
    lrgasp <- FALSE
  } 
}


for (i in 1:length(class_in)) {
  f <- class_in[[i]]
  start <- stringr::str_locate(f, dir_in)[[2]]
  end <- stringr::str_locate(f, "_classification.txt")[[1]]
  idx <- substring(f, (start+2), (end-1))
  if (idx %in% code$pipeline){
    classification <- read.table(class_in[[i]], header = T, sep = "\t")
    classification$structural_category <- gsub(classification$structural_category, pattern="Genic_", replacement = "Genic")
    classification$LRGASP_id <- gsub(classification$LRGASP_id, pattern="Genic_", replacement = "Genic")
    junctions <- read.table(junct_in[[i]], header = T, sep = "\t")
    f_in[[idx]] <- list(classification, junctions)
  }
  if (lrgasp == TRUE){
    lrgasp_f <- lrgasp.files[i]
    l_start <- stringr::str_locate(lrgasp_f, dir_in)[[2]]
    l_end <- stringr::str_locate(lrgasp_f, "_results.RData")[[1]]
    l_idx <- substring(lrgasp_f, (l_start+2), (l_end-1))
    if (l_idx %in% code$pipeline){
      lrgasp.res[[l_idx]] <- loadRData(lrgasp_f)
    }
  }
}


# ----- LRGASP input

#if (lrgasp == TRUE){
#  lrgasp.files <- 
#    list.files(dir_in,
#               pattern = "*_results.RData",
#               all.files = FALSE,
#               full.names = TRUE)
#  if (length(class_in) != length(lrgasp.files)){
#    print("ERROR: Issue loading LRGASP files")
#    print("Different number of LRGASP files than samples")
#    lrgasp <- FALSE
#  } else {
#    lrgasp.res <- list()
#    for (i in 1:length(lrgasp.files)){
#      f <- lrgasp.files[i]
#      start <- stringr::str_locate(f, dir_in)[[2]]
#      end <- stringr::str_locate(f, "_results.RData")[[1]]
#      idx <- substring(f, (start+2), (end-1))
#      if (idx %in% code$pipeline){
#          lrgasp.res[[names(f_in)[[i]]]] <- loadRData(f)
#      }
#    }
#  }
#}


# -------------------- Check for TSS and TTS genomic coords

TSS_TTS_coord <- TRUE
for ( i in 1:length(f_in)){
  if (!("TSS_genomic_coord" %in% colnames(f_in[[i]][[1]]))){
    TSS_TTS_coord <- FALSE
  }
}

# --------------------  P/A comparison and Iso & Gene analysis

res <- try({
  compareTranscriptomes(f_in)
}, silent = TRUE)

if (class(res) == "try-error"){
  print("ERROR: An error has ocurred during the comparison. A partial report will be generated. Here's the original error message:")
  print(geterrmessage())
}


# -------------------- Comparison with res classification files
print("Comparing results between samples...")

# ----- Define max number of samples in plots

if (length(f_in) <= 100){
  limit <- length(f_in)
} else {limit <- 100}


# ----- Vector of structural categories

#str_cat <- unique(f_in[[1]][[1]]$structural_category)
#str_cat <- c("full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "antisense", "fusion", "genic", "intergenic")
str_cat <- c("FSM", "ISM", "NIC", "NNC", "Antisense", "Fusion", "GenicGenomic", "GenicIntron", "Intergenic")


# ----- Generates dataframe with a summary of the SQANTI3 classification files

df_summary.1 <- data.frame(ID = names(f_in))

# Count total isoform from SQANTI3
n <- c()
for (i in f_in) {
  n <- c(n, nrow(i[[1]]))
}
df_summary.1$total <- n

# Count isoforms for each structural category
for (i in str_cat) {
  n <- c()
  for (j in f_in) {
    k <- j[[1]]
    n <- c(n, nrow(k[k$structural_category == i,]))
  }
  df_summary.1[,i] <- n
}


# ----- Generates dataframe with a summary of the unique tag comparison

df_summary.2 <- data.frame(ID = names(res[[2]][3:ncol(res[[2]])]))

# Count unique tags
n <- c()
for (i in res[[1]]) {
  n <- c(n, nrow(i))
}
df_summary.2$uniq_id <- n

# Count unique tags for each category
for (i in str_cat) {
  n <- c()
  for (j in res[[1]]) {
    n <- c(n, nrow(j[j$structural_category == i,]))
  }
  df_summary.2[,i] <- n
}


# ----- Add GenomeBrowser URL to the P/A table

df.PA <- res[[3]]
df.PA$TAGS <- lapply(df.PA$TAGS, iso2url)


# ----- Counts per gene and exon structure

countpergene <- c()
exonstructure <- c()
for (i in 1:limit){
  data <- f_in[[i]]
  data.class <- data[[1]]
  data.class <- data.class[grep("SIRV|ERCC",data.class$chrom, invert=T),]
  data.class$novelGene <- "Annotated Genes"
  data.class[grep("novelGene", data.class$associated_gene), "novelGene"] <- "Novel Genes"
  data.class$novelGene = factor(data.class$novelGene,
                                levels = c("Novel Genes","Annotated Genes"),
                                ordered=TRUE)
  
  isoPerGene = aggregate(data.class$isoform,
                         by = list("associatedGene" = data.class$associated_gene,
                                   "novelGene" = data.class$novelGene,
                                   "FSM_class" = data.class$FSM_class),
                         length)
  
  data.class[which(data.class$exons>1), "exonCat"] <- "Multi-Exon"
  data.class[which(data.class$exons==1), "exonCat"] <- "Mono-Exon"
  data.class$exonCat = factor(data.class$exonCat,
                              levels = c("Multi-Exon","Mono-Exon"),
                              ordered=TRUE)
  
  canonical.labels=c("Canonical", "Non-canonical")
  data.class$all_canonical = factor(data.class$all_canonical,
                                    labels=canonical.labels,
                                    levels = c("canonical","non_canonical"),
                                    ordered=TRUE)
  
  countpergene <- c(
    countpergene,
    sum(isoPerGene$x == 1),
    sum(isoPerGene$x == 2 | isoPerGene$x == 3),
    sum(isoPerGene$x == 4 | isoPerGene$x == 5),
    sum(isoPerGene$x >= 6)
  )
  
  exonstructure <- c(
    exonstructure,
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Novel Genes" & data.class$exonCat == "Multi-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Mono-Exon"),
    sum(data.class$novelGene == "Annotated Genes" & data.class$exonCat == "Multi-Exon")
    
  )
  
}

sample <- c(rep(names(f_in[1:limit]), each=4))
number <- rep(c("1","2-3","4-5", ">=6"), times=limit)
isoPerGene <- data.frame(sample, number, countpergene)

category <- rep(c("Novel-Monoexon", "Novel-Multiexon", "Annotated-Monoexon", "Annotated-Multiexon"), times=limit)
exonstructure <- data.frame(sample, category, exonstructure)


# ----- Summary dataframe pivoted

df_SC <- df_summary.1[1:limit,]
df_SC$total <- NULL
df_SC <- df_SC %>% 
  pivot_longer(!"ID", "SC")

# ----- Distance to TSS, TTS and CAGE peak

dist.list <- list()
dist.msr <- c("diff_to_TSS", "diff_to_TTS", "dist_to_cage_peak")
dist.SC <- c("FSM", "ISM")
contador <- 1
for (i in dist.msr){
  for (j in dist.SC){
    sample <- c()
    dist <- c()
    for (k in 1:limit){
      data.class <- f_in[[k]][[1]]
      cond <- which(data.class$structural_category == j)
      x <- data.class[cond, i]
      sample <- c(
        sample,
        rep(names(f_in)[k], times=length(x))
      )
      dist <- c(
        dist,
        x
      )
    }
    dist.list[[contador]] <- data.frame(sample, dist)
    names(dist.list)[length(dist.list)] <- paste(i,j,sep="_")
    contador <- contador + 1
  }
}

# ----- RT-switching

FSM <- c()
NIC <- c()
NNC <- c()
for (i in 1:limit){
  data.class <- f_in[[i]][[1]]
  df <- group_by(data.class, structural_category, RTS_stage) %>% dplyr::summarise(count=dplyr::n(), .groups = 'drop')
  FSM.match <- df$count[which(df$structural_category == "FSM")]
  FSM <- c(FSM, ((FSM.match[2]/(FSM.match[1]+FSM.match[2]))*100))
  
  NIC.match <- df$count[which(df$structural_category == "NIC")]
  NIC <- c(NIC, ((NIC.match[2]/(NIC.match[1]+NIC.match[2]))*100))
  
  NNC.match <- df$count[which(df$structural_category == "NNC")]
  NNC <- c(NNC, ((NNC.match[2]/(NNC.match[1]+NNC.match[2]))*100))
}

sample <- names(f_in)[1:limit]
FSM.RT <- data.frame(sample, FSM)
NIC.RT <- data.frame(sample, NIC)
NNC.RT <- data.frame(sample, NNC)

# ----- List of unique tags for each sample

l <- list()
for (i in 3:ncol(res[[2]])){
  l[[i-2]] <- res[[2]][,i] %>% na.omit()  
}
names(l) <- colnames(res[[2]])[3:ncol(res[[2]])]


# -------------------- LRGASP metrics

if (lrgasp == TRUE){
  sirv.metrics <- list()
  for (i in 1:length(lrgasp.res)){
    sirv.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["SIRV"]
  }
  sirv.metrics <- bind_cols(sirv.metrics)
  colnames(sirv.metrics) <- names(lrgasp.res)
  # FSM
  fsm.metrics <- list()
  fsm.metrics_perc <- list()
  for (i in 1:length(lrgasp.res)){
    fsm.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["FSM"] %>% as.data.frame() %>% select(FSM.Absolute.value)
    fsm.metrics_perc[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["FSM"] %>% as.data.frame() %>% select(FSM.Relative.value....)
  }
  fsm.metrics <- bind_cols(fsm.metrics)
  colnames(fsm.metrics) <- names(lrgasp.res)
  fsm.metrics_perc <- bind_cols(fsm.metrics_perc)
  colnames(fsm.metrics_perc) <- names(lrgasp.res)
  # ISM
  ism.metrics <- list()  
  ism.metrics_perc <- list()
  for (i in 1:length(lrgasp.res)){
    ism.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["ISM"] %>% as.data.frame() %>% select(ISM.Absolute.value)
    ism.metrics_perc[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["ISM"] %>% as.data.frame() %>% select(ISM.Relative.value....)
  }
  ism.metrics <- bind_cols(ism.metrics)
  colnames(ism.metrics) <- names(lrgasp.res)
  ism.metrics_perc <- bind_cols(ism.metrics_perc)
  colnames(ism.metrics_perc) <- names(lrgasp.res)

  # NIC
  nic.metrics <- list()
  nic.metrics_perc <- list()
  for (i in 1:length(lrgasp.res)){
    nic.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NIC"] %>% as.data.frame() %>% select(NIC.Absolute.value)
    nic.metrics_perc[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NIC"] %>% as.data.frame() %>% select(NIC.Relative.value....)
  }
  nic.metrics <- bind_cols(nic.metrics)
  colnames(nic.metrics) <- names(lrgasp.res)
  nic.metrics_perc <- bind_cols(nic.metrics_perc)
  colnames(nic.metrics_perc) <- names(lrgasp.res)

  # NNC
  nnc.metrics <- list()  
  nnc.metrics_perc <- list()
  for (i in 1:length(lrgasp.res)){
    nnc.metrics[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NNC"] %>% as.data.frame() %>% select(NNC.Absolute.value)
    nnc.metrics_perc[[i]] <- lrgasp.res[[names(lrgasp.res)[i]]]["NNC"] %>% as.data.frame() %>% select(NNC.Relative.value....)
  }
  nnc.metrics <- bind_cols(nnc.metrics)
  colnames(nnc.metrics) <- names(lrgasp.res)
  nnc.metrics_perc <- bind_cols(nnc.metrics_perc)
  colnames(nnc.metrics_perc) <- names(lrgasp.res)
}

#####
pa_table_sum=res$comparisonPA
pa_table_sum$found_by <- apply(pa_table_sum,1, sum_pipelines, code$pipeline)

pa_coord_merged=merge(res$iso_metrics, pa_table_sum[,c("TAGS", "structural_category", "found_by")], by.x = "LRGASP_id", by.y="TAGS")

####
length_distribution <- list()
for (i in 1:limit){
  data <- f_in[[i]]
  data.class <- data[[1]]
  data.class <- data.class[grep("SIRV|ERCC",data.class$chrom, invert=T),]
  length_distribution[[i]] <- data.class %>%  select(LRGASP_id,experiment_id, length)
}
  
length_distribution <- bind_rows(length_distribution)




#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: summary table
#    pt1.1: summary barplot
# t1.2: summary table (unique UJC)
#    pt1.2: summary barplot
# t2: presence/ausence table
# t3: SIRV metrics

# -------------------- 
# -------------------- 
# PLOT INDEX
# p1: gene characterization
#   p1.1: isoforms per gene
#   p1.2: exon structure
# p2: structural category distribution
# p3: splice junction distribution for each SC
# p4: distance to TSS
# p5: distance to TTS
# p6: distance to CAGE peak
# p7: bad quality features
# p8: good quality features
# p9: Jaccard index diagrams 
# p10: UpSet plot
# p12: Num pipelines detected UJC vs FL counts
# p13: TSS standard deviation
# p14: TTS standard deviation
# p16: iso analysis
# p17: gene analysis


print("Generating plots for the report...")

# -------------------- Global plot parameters
# COPY-PASTE FROM SQANTI3 REPORT

#myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

cat.palette = c( "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                 "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                 "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4")

cat.palette1 = c("total"="black", "FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "GenicGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "GenicIntron"="#41B6C4", "uniq_id"="black")


mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=6),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=6) ) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size=10), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15.5)) +
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "cm")) +
  theme(legend.position = "bottom")

# -------------------- 
# TABLE 1: summary table

sorted_pipelines = stringr::str_sort(code$pipeline)

#title1.1 <- grid::textGrob("Table 1.1. Summary from SQANTI3 output comparison\n", gp=gpar(fontface="italic", fontsize=10), vjust = -3.2)
#table1.1 <- tableGrob(df_summary.1[,c("ID","FSM","ISM","NIC","NNC","Antisense","Intergenic")], rows = NULL )
#t1.1 <- gTree(children=gList(table1.1, title1.1))

print("Writting Summary table for structural categories...")
write.table(df_summary.1, paste0(output_directory, "/", output_name, ".summary_table_SC.csv" ), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

melted_summary.1 <- melt(df_summary.1, id.vars = "ID")
melted_summary.1$ID <- factor(melted_summary.1$ID, levels=sorted_pipelines)
pt1.1 <- ggplot2::ggplot(melted_summary.1, aes(x=ID, y=value, group=variable)) + 
  geom_line(linetype = "dashed", aes(color=variable)) +
  geom_point(aes(color=variable)) +
  mytheme + scale_color_manual(values=cat.palette1, name="Structural Category")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Isoforms detected by pipeline",
       x="", y="Counts")



#title1.2 <- grid::textGrob("Table 1.2. Summary from SQANTI3 output comparison\n Only unique UJC", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
#table1.2 <- tableGrob(df_summary.2[,c("ID","FSM","ISM","NIC","NNC","Antisense","Intergenic")], rows = NULL)
#t1.2 <- gTree(children=gList(table1.2, title1.2))

print("Writting Summary table for structural categories. Only unique UJC...")
write.table(df_summary.2, paste0(output_directory, "/", output_name, ".summary_table_UJC_SC.csv" ), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

melted_summary.2 <- melt(df_summary.2, id.vars = "ID")
melted_summary.2$ID <- factor(melted_summary.2$ID, levels=sorted_pipelines)
pt1.2 <- ggplot2::ggplot(melted_summary.2, aes(x=ID, y=value, group=variable)) + 
  geom_line(linetype = "dashed", aes(color=variable)) +
  geom_point(aes(color=variable)) +
  mytheme + scale_color_manual(values=cat.palette1, name="Structural Category")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "UJC detected by pipeline",
       x="", y="Counts UJC")




# TABLE 2: presence/ausence isoforms

#t2 <- DT::datatable(df.PA,
#              escape = FALSE,
#              options = list(
#                pageLength = 10,
#                autoWidth = TRUE,
#                columnDefs = list(list(width = '10px', targets = "_all"))
#              ),
#              rownames = FALSE,
#              caption = "Table 3. Presence/Ausence of all the isoform models")

# -------------------- 
# TABLE 3: SIRV metrics



if (lrgasp == TRUE){
  ## table theme
  table_theme=ttheme_default(base_size = 8, base_colour = "black", base_family = "",
                 parse = FALSE, padding = unit(c(2, 2), "mm"))
  ## SIRVs
  #title3 <- grid::textGrob("Table 3. SIRV metrics comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table3 <- tableGrob(sirv.metrics, theme = table_theme)
  #t3 <- gTree(children=gList(table3, title3))
  
  print("Writting SIRVs metrics table...")
  write.table(sirv.metrics, paste0(output_directory, "/", output_name, ".SIRVS_metrics.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t3 <- melt(as.matrix(sirv.metrics), id.vars = "X1")
  colnames(melted_t3) <- c("Metric", "Pipelines", "value")
  melted_t3$value <- as.numeric(melted_t3$value)
  melted_t3$Pipelines <- factor(melted_t3$Pipelines, levels=sorted_pipelines)
  pt3 <- list()
  for (i in 1:length(rownames(sirv.metrics))) {
    metric <- rownames(sirv.metrics)[i]
    selected_variable <- melted_t3 %>% filter(Metric==metric)
    pt3[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle = 90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("SIRVs:\n\n", metric),
           x="", y=metric) +
      scale_fill_conesa(palette = "complete") + theme(legend.position = "none")
  }
  
  ## FSM
  #title4.1 <- grid::textGrob("Table 4.1 FSM metrics comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table4.1 <- tableGrob(fsm.metrics, theme = table_theme)
  #t4.1 <- gTree(children=gList(table4.1, title4.1))
  
  print("Writting FSMs metrics table...")
  write.table(fsm.metrics, paste0(output_directory, "/", output_name, ".FSM_metrics.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t4.1 <- melt(as.matrix(fsm.metrics), id.vars = c("X1","X2"))
  colnames(melted_t4.1) <- c("Metric", "Pipelines", "value")
  melted_t4.1$value <- as.numeric(melted_t4.1$value)
  melted_t4.1$Pipelines <- factor(melted_t4.1$Pipelines, levels=sorted_pipelines)
  pt4.1 <- list()
  for (i in 1:length(rownames(fsm.metrics))) {
    metric <- rownames(fsm.metrics)[i]
    selected_variable <- melted_t4.1 %>% filter(Metric==metric)
    if (metric != "Reference redundancy Level"){
      pt4.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
        mytheme + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = paste0("FSM:\n\n", metric),
             x="", y="Count") +
        scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
    }else{
      pt4.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
        mytheme + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = paste0("FSM:\n\n", metric),
             x="", y="Redundancy level") +
        scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
    }
  }
  # FSM perc
  #title4.2 <- grid::textGrob("Table 4.2 FSM metrics percentages comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table4.2 <- tableGrob(fsm.metrics_perc, theme = table_theme)
  #t4.2 <- gTree(children=gList(table4.2, title4.2))
  
  print("Writting FSMs metrics table (%)...")
  write.table(fsm.metrics_perc, paste0(output_directory, "/", output_name, ".FSM_metrics_perc.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t4.2 <- melt(as.matrix(fsm.metrics_perc), id.vars = c("X1","X2"))
  colnames(melted_t4.2) <- c("Metric", "Pipelines", "value")
  melted_t4.2$value <- as.numeric(melted_t4.2$value)
  melted_t4.2$Pipelines <- factor(melted_t4.2$Pipelines, levels=sorted_pipelines)
  pt4.2 <- list()
  for (i in 1:length(rownames(fsm.metrics_perc))) {
    metric <- rownames(fsm.metrics_perc)[i]
    selected_variable <- melted_t4.2 %>% filter(Metric==metric)
    pt4.2[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("FSM:\n\n", metric),
           x="", y="Percentage (%)") +
      scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }
  
  
  ## ISM
  #title5.1 <- grid::textGrob("Table 5.1 ISM metrics comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table5.1 <- tableGrob(ism.metrics, theme = table_theme)
  #t5.1 <- gTree(children=gList(table5.1, title5.1))
  
  print("Writting ISMs metrics table...")
  write.table(ism.metrics, paste0(output_directory, "/", output_name, ".ISM_metrics.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t5.1 <- melt(as.matrix(ism.metrics), id.vars = c("X1","X2"))
  colnames(melted_t5.1) <- c("Metric", "Pipelines", "value")
  melted_t5.1$value <- as.numeric(melted_t5.1$value)
  melted_t5.1$Pipelines <- factor(melted_t5.1$Pipelines, levels=sorted_pipelines)
  pt5.1 <- list()
  for (i in 1:length(rownames(ism.metrics))) {
    metric <- rownames(ism.metrics)[i]
    selected_variable <- melted_t5.1 %>% filter(Metric==metric)
    if (metric != "Reference redundancy Level"){
      pt5.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
        mytheme + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = paste0("ISM:\n\n", metric),
             x="", y="Count") +
        scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
    }else{
      pt4.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
        mytheme + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
        labs(title = paste0("ISM:\n\n", metric),
             x="", y="Redundancy level") +
        scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
    }
 }
  # ISM perc
  #title5.2 <- grid::textGrob("Table 5.2 ISM metrics percentages comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table5.2 <- tableGrob(ism.metrics_perc, theme = table_theme)
  #t5.2 <- gTree(children=gList(table5.2, title5.2))
  
  print("Writting ISMs metrics table (%)...")
  write.table(ism.metrics_perc, paste0(output_directory, "/", output_name, ".ISM_metrics_perc.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t5.2 <- melt(as.matrix(ism.metrics_perc), id.vars = c("X1","X2"))
  colnames(melted_t5.2) <- c("Metric", "Pipelines", "value")
  melted_t5.2$value <- as.numeric(melted_t5.2$value)
  melted_t5.2$Pipelines <- factor(melted_t5.2$Pipelines, levels=sorted_pipelines)
  pt5.2 <- list()
  for (i in 1:length(rownames(ism.metrics_perc))) {
    metric <- rownames(ism.metrics_perc)[i]
    selected_variable <- melted_t5.2 %>% filter(Metric==metric)
    pt5.2[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("ISM:\n\n", metric),
           x="", y="Percentage (%)")+
      scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }

  ## NIC
  #title6.1 <- grid::textGrob("Table 6.1 NIC metrics comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table6.1 <- tableGrob(nic.metrics, theme = table_theme)
  #t6.1 <- gTree(children=gList(table6.1, title6.1))
  
  print("Writting NICs metrics table...")
  write.table(nic.metrics, paste0(output_directory, "/", output_name, ".NIC_metrics.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t6.1 <- melt(as.matrix(nic.metrics), id.vars = c("X1","X2"))
  colnames(melted_t6.1) <- c("Metric", "Pipelines", "value")
  melted_t6.1$value <- as.numeric(melted_t6.1$value)
  melted_t6.1$Pipelines <- factor(melted_t6.1$Pipelines, levels=sorted_pipelines)
  pt6.1 <- list()
  for (i in 1:length(rownames(nic.metrics))) {
    metric <- rownames(nic.metrics)[i]
    selected_variable <- melted_t6.1 %>% filter(Metric==metric)
    pt6.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("NIC:\n\n", metric),
           x="", y="Count")+
      scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }
  # NIC perc
  #title6.2 <- grid::textGrob("Table 6.2 NIC metrics percentages comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table6.2 <- tableGrob(nic.metrics_perc, theme = table_theme)
  #t6.2 <- gTree(children=gList(table6.2, title6.2))
  
  print("Writting NICs metrics table (%)...")
  write.table(nic.metrics_perc, paste0(output_directory, "/", output_name, ".NIC_metrics_perc.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t6.2 <- melt(as.matrix(nic.metrics_perc), id.vars = c("X1","X2"))
  colnames(melted_t6.2) <- c("Metric", "Pipelines", "value")
  melted_t6.2$value <- as.numeric(melted_t6.2$value)
  melted_t6.2$Pipelines <- factor(melted_t6.2$Pipelines, levels=sorted_pipelines)
  pt6.2 <- list()
  for (i in 1:length(rownames(nic.metrics_perc))) {
    metric <- rownames(nic.metrics_perc)[i]
    selected_variable <- melted_t6.2 %>% filter(Metric==metric)
    pt6.2[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("NIC:\n\n", metric),
           x="", y="Percentage (%)")+
      scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }
  
  ## NNC
  #title7.1 <- grid::textGrob("Table 7.1 NNC metrics comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table7.1 <- tableGrob(nnc.metrics, theme = table_theme)
  #t7.1 <- gTree(children=gList(table7.1, title7.1))
  
  print("Writting NNCs metrics table...")
  write.table(nnc.metrics, paste0(output_directory, "/", output_name, ".NNC_metrics.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t7.1 <- melt(as.matrix(nnc.metrics), id.vars = c("X1","X2"))
  colnames(melted_t7.1) <- c("Metric", "Pipelines", "value")
  melted_t7.1$value <- as.numeric(melted_t7.1$value)
  melted_t7.1$Pipelines <- factor(melted_t7.1$Pipelines, levels=sorted_pipelines)
  pt7.1 <- list()
  for (i in 1:length(rownames(nnc.metrics))) {
    metric <- rownames(nnc.metrics)[i]
    selected_variable <- melted_t7.1 %>% filter(Metric==metric)
    pt7.1[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("NNC:\n\n ", metric),
           x="", y="Count")+
        scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }
  # NNC perc
  #title7.2 <- grid::textGrob("Table 7.2 NNC metrics percentages comparison\n", gp=gpar(fontface="italic", fontsize=17), vjust = -3.2)
  #table7.2 <- tableGrob(nnc.metrics_perc, theme = table_theme)
  #t7.2 <- gTree(children=gList(table7.2, title7.2))
  
  print("Writting NNC metrics table (%)...")
  write.table(nnc.metrics_perc, paste0(output_directory, "/", output_name, ".NNC_metrics_perc.csv" ), sep=",", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  melted_t7.2 <- melt(as.matrix(nnc.metrics_perc), id.vars = c("X1","X2"))
  colnames(melted_t7.2) <- c("Metric", "Pipelines", "value")
  melted_t7.2$value <- as.numeric(melted_t7.2$value)
  melted_t7.2$Pipelines <- factor(melted_t7.2$Pipelines, levels=sorted_pipelines)
  pt7.2 <- list()
  for (i in 1:length(rownames(nnc.metrics_perc))) {
    metric <- rownames(nnc.metrics_perc)[i]
    selected_variable <- melted_t7.2 %>% filter(Metric==metric)
    pt7.2[[i]] <-ggplot2::ggplot(selected_variable, aes(x=Pipelines, y=value, fill=Pipelines)) + 
      geom_bar(stat="identity") +
      geom_text(aes(label=value), hjust=1.1, size=3, angle=90)+
      mytheme + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
      labs(title = paste0("NNC:\n\n", metric),
           x="", y="Percentage (%)")+
      scale_fill_conesa(palette = "complete")+ theme(legend.position = "none")
  }
  
  
  
}


# -------------------- 
# PLOT 1: gene characterization
# PLOT 1.1: isoforms per gene
isoPerGene$sample <- factor(isoPerGene$sample, levels=sorted_pipelines)

p1.1 <- ggplot(isoPerGene, aes(fill=number, y=countpergene, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  mytheme + scale_fill_manual(values = myPalette, name="Num. isoforms per gene") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Number of isoforms per gene",
       x="", y="Percentage of genes")

# PLOT 1.2: exon structure
exonstructure$sample <- factor(exonstructure$sample, levels=sorted_pipelines)

p1.2 <- ggplot(exonstructure, aes(fill=category, y=exonstructure, x=sample)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = myPalette, name="") + mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  labs(title = "Type of isoforms", x="", y="Percentage (%)")

# PLOT 2: structural category distribution

df_SC$ID <- factor(df_SC$ID, levels=sorted_pipelines)
p2 <- ggplot(df_SC, aes(fill=SC, y=value, x=ID)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = cat.palette, name="Structural Category") + mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distribution of Structural Categories", x="", y="Percentage of isoforms")


# PLOT 4: distance to TSS
for (i in 1:length(dist.list)){
  dist.list[[i]]$sample <- factor(dist.list[[i]]$sample, levels=sorted_pipelines)
}
p4.1 <- ggplot(dist.list[[1]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distance to TSS of FSM isoforms", x="", y="Distance (bp)") +
  ylim(c(-200,200))  +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

p4.2 <- ggplot(dist.list[[2]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distance to TSS of ISM isoforms", x="", y="Distance (bp)") +
  ylim(c(-2000,2000))  +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

# PLOT5: distance to TTS
p5.1 <- ggplot(dist.list[[3]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distance to TTS of FSM isoforms", x="", y="Distance (bp)") +
  ylim(c(-200,200)) +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none") 

p5.2 <- ggplot(dist.list[[4]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distance to TTS of ISM isoforms", x="", y="Distance (bp)") +
  ylim(c(-2000,2000)) +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

# PLOT 6: distance to CAGE peak
p6.1 <- ggplot(dist.list[[5]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  ylim(c(-1000,200))+
  labs(title = "Distance to closest CAGE peak (FSM only)", y="Distance (bp)", x="") +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

p6.2 <- ggplot(dist.list[[6]], aes(x=sample, y=(dist*-1))) + 
  geom_violin(alpha=0.7, aes(fill=sample)) +
  geom_boxplot(width=0.05, outlier.shape = NA) +
  mytheme +
  ylim(c(1000,-200))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Distance to closest CAGE peak (ISM only)", y="Distance (bp)", x="") +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")


# PLOT 7: bad quality features
# PLOT 7.1: RT-switching
FSM.RT$sample <- factor(FSM.RT$sample, levels=sorted_pipelines)
p7.1.1 <- ggplot(FSM.RT, aes(x=sample, y=FSM, fill=sample)) + geom_bar(stat="identity") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "RT-Switching incidence on FSM isoforms", x="", y="Percentage of isoforms") +
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

NIC.RT$sample <- factor(NIC.RT$sample, levels=sorted_pipelines)
p7.1.2 <- ggplot(NIC.RT, aes(x=sample, y=NIC, fill=sample)) + geom_bar(stat="identity") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "RT-Switching incidence on NIC isoforms", x="", y="Percentage of isoforms")+
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

NNC.RT$sample <- factor(NNC.RT$sample, levels=sorted_pipelines)
p7.1.3 <- ggplot(NNC.RT, aes(x=sample, y=NNC, fill=sample)) + geom_bar(stat="identity") +
  mytheme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(title = "RT-Switching incidence on NNC isoforms", x="", y="Percentage of isoforms")+ 
  scale_fill_conesa(palette = "complete") + theme(legend.position = "none")

# PLOT 8: Length Distribution
if (max(length_distribution$length)>30000){
  p8 <- ggplot(length_distribution, aes(x=experiment_id, y=length))+
    geom_violin(alpha=0.7 , aes(fill=experiment_id)) +
    geom_boxplot(width=0.05, outlier.shape = NA) +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    ylim(c(0,30000)) +
    labs(title = "Transcript length distribution", x="", y="Length (bp)") + 
    scale_fill_conesa(palette = "complete") + theme(legend.position = "none") 
  
}else{
  p8 <- ggplot(length_distribution, aes(x=experiment_id, y=length))+
    geom_violin(alpha=0.7 , aes(fill=experiment_id)) +
    geom_boxplot(width=0.05, outlier.shape = NA) +
    mytheme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(title = "Transcript length distribution", x="", y="Length (bp)") + 
    scale_fill_conesa(palette = "complete") + theme(legend.position = "none") 
}


# PLOT 9: Jaccard Index

jac_index_mat <- get_jaccard_matrix(pa_table_sum, as.character(code$pipeline))

#mean_jac <- jac_index_mat[lower.tri(jac_index_mat, diag = F)] %>% median()


# PLOT 10: UpSet plot

if (opt$upset){
  upset_data = ComplexHeatmap::list_to_matrix(l)
  comb_mat = ComplexHeatmap::make_comb_mat(l, mode = "distinct")
  code = code[match(set_name(comb_mat), code$pipeline),]
  sorted_code =  stringr::str_sort(code$pipeline)
  print("Upset all...")
  p10 <- 
    UpSet(comb_mat, 
          pt_size=unit(1, "mm") , lwd = 0.5,
          comb_order = order(comb_size(comb_mat), decreasing = T),
          set_order = sorted_code,
          row_title = "ALL",
          top_annotation = HeatmapAnnotation(
            degree = as.character(comb_degree(comb_mat)),
            "Intersection\nsize" = anno_barplot(comb_size(comb_mat), 
                                                border = FALSE, 
                                                gp = gpar(fill = "black"),
                                                height = unit(10, "cm")),
            annotation_name_side = "left"),
          right_annotation = rowAnnotation(
            "Set size" = anno_barplot(set_size(comb_mat), 
                                      border = FALSE, 
                                      gp = gpar(fill = "black")),
            group=code$Library 
          ),
          row_names_max_width = unit(1.5, "cm"),
          row_names_gp = gpar(fontsize = 10)
    )
  
  p10.1 <- 
    UpSet(comb_mat, 
          pt_size=unit(1, "mm") , lwd = 0.5,
          comb_order = order(comb_size(comb_mat), decreasing = T),
          set_order = sorted_code,
          row_title = "ALL",
          top_annotation = HeatmapAnnotation(
            degree = as.character(comb_degree(comb_mat)),
            "Intersection\nsize" = anno_barplot(comb_size(comb_mat), 
                                                border = FALSE, 
                                                gp = gpar(fill = "black"),
                                                height = unit(10, "cm")),
            annotation_name_side = "left"),
          right_annotation = rowAnnotation(
            "Set size" = anno_barplot(set_size(comb_mat), 
                                      border = FALSE, 
                                      gp = gpar(fill = "black")),
            group=code$Platform 
          ),
          row_names_max_width = unit(1.5, "cm"),
          row_names_gp = gpar(fontsize = 10)
    )
  ## Upset for genes
  print("Upset Genes...")
  gd = res$genes_detected
  upset_data = ComplexHeatmap::list_to_matrix(gd)
  comb_mat = ComplexHeatmap::make_comb_mat(gd, mode = "distinct")
  sorted_code =  stringr::str_sort(code$pipeline)
  
  p11 <-
    UpSet(comb_mat, 
          pt_size=unit(1, "mm") , lwd = 0.5,
          comb_order = order(comb_size(comb_mat), decreasing = T),
          set_order = sorted_code,
          row_title = "Known genes detected",
          top_annotation = HeatmapAnnotation(
            degree = as.character(comb_degree(comb_mat)),
            "Intersection\nsize" = anno_barplot(comb_size(comb_mat), 
                                                border = FALSE, 
                                                gp = gpar(fill = "black"),
                                                height = unit(10, "cm")),
            annotation_name_side = "left"),
          right_annotation = rowAnnotation(
            "Set size" = anno_barplot(set_size(comb_mat), 
                                      border = FALSE, 
                                      gp = gpar(fill = "black")),
            group=code$Library 
          ),
          row_names_max_width = unit(1.5, "cm"),
          row_names_gp = gpar(fontsize = 10)
    )
}

# PLOT 12: Num pipleines that found a UJC vsd mean FL counts
pa_coord_merged_wo_SIRV=pa_coord_merged[!grepl("SIRV", pa_coord_merged$LRGASP_id),]
p12 <- ggplot(pa_coord_merged_wo_SIRV, aes(x=found_by)) +
  geom_histogram(aes(fill=structural_category), binwidth = 1) + 
  mytheme +
  labs(x="Num. Pipelines ", y="Count", 
       title="Number of pipelines that found a certain UJC",
       subtitle = "Coloured by Structural Categories")+
  scale_fill_manual(values=cat.palette, name="Structural Category")  +
  stat_bin(aes(y=..count.. , label=..count..), geom="text", size=2, vjust=-0.2, angle=0, binwidth = 1) 

p12.1 <- list()
sc_levels <- levels(as.factor(pa_coord_merged_wo_SIRV$structural_category))
for (i in 1:length(sc_levels)) {
    metric <- sc_levels[i]
    selected_variable <- pa_coord_merged_wo_SIRV %>% filter(structural_category==metric)
    p12.1[[i]] <- ggplot(selected_variable, aes(x=factor(found_by), y=log10(FL_cpm)))+
      geom_violin(alpha=0.7, size=0.5, aes(fill=structural_category)) +
      geom_boxplot(width=0.05, outlier.shape = NA) +
      mytheme +
      labs(x="Num. Pipelines", y="log10(median CPM)", 
      title="Number of pipelines that found \n\na certain UJC vs median CPM",
      subtitle=metric) +
      scale_fill_manual(values=cat.palette[metric])
}


p12.2 <- ggplot(pa_coord_merged_wo_SIRV, aes(x=factor(found_by), y=log10(FL_cpm)))+
  geom_violin(alpha=0.7, aes(fill=structural_category)) +
  geom_boxplot(width=0.05, outlier.shape = NA) + 
  mytheme +
  labs(x="Num. Pipelines ", y="log10(median CPM)", 
       title="Number of pipelines that found \n\na certain UJC vs median CPM") +
  scale_fill_manual(values=cat.palette) +
  facet_wrap(~ structural_category, nrow=3)


if (TSS_TTS_coord == TRUE) {
  
  # PLOT 13: TSS standard deviation per pipeline
  
  a <- bind_rows(res$classifications, .id = "experiment_id")
  a$experiment_id <- factor(a$experiment_id, levels=sorted_pipelines)
  p13 <- ggplot(a, aes(y=experiment_id, x=(sdTSS))) +
    geom_violin(alpha=0.7, aes(fill=experiment_id)) +
    geom_boxplot(width=0.05, outlier.shape = NA) +
    mytheme +
    labs(title = "Standard deviation of genomic TSS coordinates", x="SD", y="") +
    scale_fill_conesa(palette = "complete") + theme(legend.position = "none") 
  
  
  # PLOT 14: TTS standard deviation per pipeline
  p14 <- ggplot(a, aes(y=experiment_id, x=(sdTTS))) +
    geom_violin(alpha=0.7, aes(fill=experiment_id)) +
    geom_boxplot(width=0.05, outlier.shape = NA) +
    mytheme +
    labs(title = "Standard deviation of genomic TTS coordinates", x="SD", y="") +
    scale_fill_conesa(palette = "complete") + theme(legend.position = "none") 
  
}

# -------------------- Output report

print("Generating the PDF report...")
script.path <- getwd()
source(paste(script.path, "generatePDFreport.R", sep = "/"))

invisible(generatePDFreport())

#Sys.setenv(RSTUDIO_PANDOC = "/usr/lib/rstudio/bin/pandoc")
#rmarkdown::render(
#  input = paste(getwd(), "SQANTI3_comparison_report.Rmd", sep = "/"),
#  output_dir = output_directory,
#  output_file = paste0(output_name, ".html")
#)

# -------------------- Output csv

print("Writting CSV file with P/A table of UJC...")
write.table(res$comparisonPA, paste0(output_directory, "/", output_name, ".pa.csv" ), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

print("Writting pre-BED...")
write.table(res$iso_metrics, paste0(output_directory, "/", output_name, ".UJC_info.csv" ), sep=",", quote = FALSE, row.names = FALSE, col.names = TRUE)

print("Saving the environment...")
save.image(file=paste0(output_directory,"/",output_name,".environment.RData"))

print("DONE\nExecution finished")

