# libraries
library(tidyverse)
library(stringr)
library(Rmpfr)

# data arguments
args <- commandArgs(trailingOnly = TRUE)
pa_file <- args[1]
coord_file <- args[2]
code_file <- args[3]
name_bed <- args[4]
# read data
pa_table=read_csv(pa_file)
coord_table=read_csv(coord_file , col_types = "cccccccc") 
code=read_csv(code_file)

# read code list
# platforms (freestyle discarded)
pb_pipelines <- code %>% filter(Platform %in% c("PB","PacBio", "PacBio+Illumina", "PB+Illumina")) %>% select(pipelineCode) 
pb_pipelines <- pb_pipelines$pipelineCode %>% as.character()
ont_pipelines <- code %>% filter(Platform %in% c("ONT", "ONT+Illumina")) %>% select(pipelineCode)
ont_pipelines <- ont_pipelines$pipelineCode %>% as.character()

# libraries
cdna_pipelines <- code %>% filter(Library_Preps %in% c("cDNA")) %>% select(pipelineCode) 
cdna_pipelines <- cdna_pipelines$pipelineCode %>% as.character()
drna_pipelines <- code %>% filter(Library_Preps %in% c("dRNA")) %>% select(pipelineCode) 
drna_pipelines <- drna_pipelines$pipelineCode %>% as.character()
r2c2_pipelines <- code %>% filter(Library_Preps %in% c("R2C2")) %>% select(pipelineCode) 
r2c2_pipelines <- r2c2_pipelines$pipelineCode %>% as.character()
captrap_pipelines <- code %>% filter(Library_Preps %in% c("CapTrap")) %>% select(pipelineCode) 
captrap_pipelines <- captrap_pipelines$pipelineCode %>% as.character()
# data categories
lo_pipelines <- code %>% filter(Data_Category %in% c("Long only", "LO")) %>% select(pipelineCode) 
lo_pipelines <- lo_pipelines$pipelineCode %>% as.character()
ls_pipelines <- code %>% filter(Data_Category %in% c("Long short", "LS")) %>% select(pipelineCode) 
ls_pipelines <- ls_pipelines$pipelineCode %>% as.character()
fs_pipelines <- code %>% filter(Data_Category %in% c("Freestyle","FS")) %>% select(pipelineCode)
fs_pipelines <- fs_pipelines$pipelineCode %>% as.character()
# platform and lib
ont_cdna_pipelines <- code %>% filter(Library_Preps %in% c("cDNA")) %>% filter(Platform %in% c("ONT", "ONT+Illumina")) %>% select(pipelineCode) 
ont_cdna_pipelines <- ont_cdna_pipelines$pipelineCode %>% as.character()
ont_drna_pipelines <- code %>% filter(Library_Preps %in% c("dRNA")) %>% filter(Platform %in% c("ONT", "ONT+Illumina")) %>% select(pipelineCode) 
ont_drna_pipelines <- ont_drna_pipelines$pipelineCode %>% as.character()
ont_r2c2_pipelines <- code %>% filter(Library_Preps %in% c("R2C2")) %>% filter(Platform %in% c("ONT", "ONT+Illumina")) %>% select(pipelineCode) 
ont_r2c2_pipelines <- ont_r2c2_pipelines$pipelineCode %>% as.character()
ont_captrap_pipelines <- code %>% filter(Library_Preps %in% c("CapTrap")) %>% filter(Platform %in% c("ONT", "ONT+Illumina")) %>% select(pipelineCode) 
ont_captrap_pipelines <- ont_captrap_pipelines$pipelineCode %>% as.character()

pb_cdna_pipelines <- code %>% filter(Library_Preps %in% c("cDNA")) %>% filter(Platform %in% c("PB","PacBio", "PacBio+Illumina", "PB+Illumina")) %>% select(pipelineCode) 
pb_cdna_pipelines <- pb_cdna_pipelines$pipelineCode %>% as.character()
pb_captrap_pipelines <- code %>% filter(Library_Preps %in% c("CapTrap")) %>% filter(Platform %in% c("PB","PacBio", "PacBio+Illumina", "PB+Illumina")) %>% select(pipelineCode) 
pb_captrap_pipelines <- pb_captrap_pipelines$pipelineCode %>% as.character()
## functions

get_bed <- function(coord){
  UJC_split <- unlist(strsplit(as.character(coord[1]), '_'))
  chrom_type="normal"
  if (as.character(UJC_split[3])!= "+" & as.character(UJC_split[3])!= "-"){
    if (startsWith(as.character(UJC_split[4]), "random")){
      chromosome <- paste(c(as.character(UJC_split[2]),as.character(UJC_split[3]),as.character(UJC_split[4])),
                          collapse = "_")
      strand <- as.character(UJC_split[5])
      chrom_type="random"
    }else{
      chromosome <- paste(c(as.character(UJC_split[2]),as.character(UJC_split[3])), collapse = "_")
      strand <- as.character(UJC_split[4])
      chrom_type="unknown"
    }
  }else{
    chromosome <- as.character(UJC_split[2])
    strand <- as.character(UJC_split[3])
  }
  coord["exons"]=as.integer(coord["exons"])
  tss=as.bigz(gsub("\\..*","",coord['mean.median.TSS']))
  tts=as.bigz(gsub("\\..*","",coord['mean.median.TTS']))
  start=min(tss,tts)
  end=max(tss,tts)
  if (coord["exons"]==1){
    l=end-start
    all <- c(
      chromosome,
      as.character(start),
      as.character(end),
      as.character(coord['LRGASP_id']),
      '40',
      strand, 
      as.character(start),
      as.character(start),
      "255,0,0",
      as.character(coord["exons"]),
      as.character(l),
      "0",
      as.character(coord["gene"])
    )
  }else{
    if (chrom_type == "normal"){
      blocks=as.bigz(c(start,(UJC_split[seq(5,length(UJC_split),by = 2)])))
      block_ends=c(as.bigz(UJC_split[seq(4,length(UJC_split),by = 2)]), end)
    }else if (chrom_type == "unknown"){
      blocks=as.bigz(c(start,(UJC_split[seq(6,length(UJC_split),by = 2)])))
      block_ends=c(as.bigz(UJC_split[seq(5,length(UJC_split),by = 2)]), end)
    } else {
      blocks=as.bigz(c(start,(UJC_split[seq(7,length(UJC_split),by = 2)])))
      block_ends=c(as.bigz(UJC_split[seq(6,length(UJC_split),by = 2)]), end)
    }
    block_starts=blocks - (start)
    block_sizes=block_ends - (blocks+1)
    if (block_sizes[length(block_sizes)]==0){
      block_sizes[-1]=1
      end=end+1
    }
    block_starts=as.character(block_starts)
    block_sizes=as.character(block_sizes)
    block_sizes=paste(block_sizes, collapse = ",")
    block_starts=paste(block_starts, collapse=",")
    all <- c(
      chromosome,
      as.character(start),
      as.character(end-1),
      as.character(coord['LRGASP_id']),
      '40',
      strand,
      as.character(start),
      as.character(start),
      "255,0,0",
      as.character(coord["exons"]),
      block_sizes,
      block_starts,
      as.character(coord["gene"])
    )
  }
  all <- paste(all,collapse = ";") 
#  colnames(all) <- c("chrom", "start", "end", "id", "Q", "strand", "thickStart", "thickEnd", "rgb", "blockCount", "blockSizes", "blockStarts")
  return(all)
}

get_barcode <- function(pa){
  n_pb=pa[pb_pipelines] %>% as.numeric() %>% sum()
  n_ont=pa[ont_pipelines] %>% as.numeric() %>% sum()
  n_fs=pa[fs_pipelines] %>% as.numeric() %>% sum()
  
  n_cdna=pa[cdna_pipelines] %>% as.numeric() %>% sum()
  n_drna=pa[drna_pipelines] %>% as.numeric() %>% sum()
  n_captrap=pa[captrap_pipelines] %>% as.numeric() %>% sum()
  n_r2c2=pa[r2c2_pipelines] %>% as.numeric() %>% sum()
  
  n_lo=pa[lo_pipelines] %>% as.numeric() %>% sum()
  n_ls=pa[ls_pipelines] %>% as.numeric() %>% sum()
  
  n_ont_cdna=pa[ont_cdna_pipelines] %>% as.numeric() %>% sum()
  n_ont_drna=pa[ont_drna_pipelines] %>% as.numeric() %>% sum()
  n_ont_captrap=pa[ont_captrap_pipelines] %>% as.numeric() %>% sum()
  n_ont_r2c2=pa[ont_r2c2_pipelines] %>% as.numeric() %>% sum()
  n_pb_cdna=pa[pb_cdna_pipelines] %>% as.numeric() %>% sum()
  n_pb_captrap=pa[pb_captrap_pipelines] %>% as.numeric() %>% sum()
  
  plat=paste(c(n_pb,n_ont,n_fs), collapse=",")
  lib=paste(c(n_cdna,n_drna,n_captrap,n_r2c2,n_fs), collapse = ",")
  dat=paste(c(n_lo,n_ls,n_fs), collapse = ",")
  plat_lib=paste(c(n_pb_cdna,n_pb_captrap,n_ont_cdna,n_ont_drna,n_ont_captrap,n_ont_r2c2), collapse = ",")
  bc=str_c(plat,lib,dat,plat_lib,
           n_pb, n_ont, n_fs,
           n_cdna, n_drna, n_captrap, n_r2c2, n_fs,
           n_lo, n_ls, n_fs,
           sep = ";")
  return(bc)
}


### actual script
#coord_table <- as.data.frame(coord_table)
bed_file <- apply(coord_table,1,get_bed)
bed_file <- str_split(bed_file, pattern = ";", simplify = T) %>% as.data.frame()
colnames(bed_file) <- c("chrom", "start", "end", "id", "Q", "strand", "thickStart", "thickEnd", "rgb", "blockCount", "blockSizes", "blockStarts", "associated_gene")
rownames(bed_file) <- bed_file$id
bed_file$num = seq_along(bed_file$id)
bed_file$associated_gene = apply(bed_file,1, function(x){
  g <- x["associated_gene"]
  if (startsWith(g, "novelGene")){
    g_id <- paste("novelGene",  trimws(as.character(x["num"])), sep = "-")
  }else{
    if (nchar(g)>250){
      tmp_genes <- strsplit(g,split="-") %>% unlist()
      g <- paste(tmp_genes[1:3], collapse="-")
    }
    g_id <- g
  }
  return(g_id)
})

### barcode
pa_table <- pa_table %>% as.data.frame()
rownames(pa_table) <- pa_table$TAGS
pa_table$barcode <- apply(pa_table,1,get_barcode)
split_bc <- str_split(pa_table$barcode, pattern = ";", simplify = T) %>% as.data.frame()
rownames(split_bc) <- pa_table$TAGS
colnames(split_bc) <- c("Platform", "Library", "Data", "Platform_Library",
                        "PB_count", "ONT_count", "PlatFS_count",
                        "cDNA_count", "dRNA_count", "CapTrap_count", "R2C2_count", "LibFS_count",
                        "LO_count", "LS_count", "DataCatFS_count")

final_bed=merge(bed_file, split_bc, by=0)
coord_table$num = seq_along(coord_table$LRGASP_id)
coord_table <- coord_table %>%  as.data.frame()
coord_table$small_id = apply(coord_table,1, function(x){
  UJC_split <- strsplit(as.character(x[1]), '_') %>%  unlist()
  paste(c(UJC_split[1], trimws(as.character(x["num"]))), collapse = "-")
  })
coord_table[which(is.na(coord_table$SD.TSS)),"SD.TSS"] <- "NaN"
coord_table[which(is.na(coord_table$SD.TTS)),"SD.TTS"] <- "NaN"

final_bed=merge(final_bed, coord_table[,c("LRGASP_id", "SD.TSS", "SD.TTS", "exons", "length", "FL_cpm", "small_id")], by.x="id", by.y="LRGASP_id")

# Change RGB color so they will represent the relative abundance of the isoform
RGB_SCALE = c('0,0,0', '26,0,0', '51,0,0', '77,0,0', '102,0,0',
             '128,0,0', '153,0,0', '179,0,0', '204,0,0', '230,0,0',
             '255,0,0', '255,26,26', '255,51,51', '255,77,77', '255,102,102',
             '255,128,128', '255,153,153', '255,179,179', '255,204,204', '255,230,230')
RGB_SCALE = rev(RGB_SCALE)
NUM_RGB = length(RGB_SCALE)

final_bed <- final_bed %>% 
  group_by(associated_gene) %>%
  mutate(new_rgb=cut(log10(as.numeric(FL_cpm)), breaks=NUM_RGB, labels=RGB_SCALE)) %>% 
  as.data.frame()

final_bed=final_bed[,c("chrom", "start", "end", "small_id", "Q", "strand", "thickStart", "thickEnd", "new_rgb", "blockCount", "blockSizes", "blockStarts",
                       "Platform", "Library", "Data","Platform_Library", "SD.TSS", "SD.TTS", "exons", "length", "FL_cpm", "associated_gene", 
                       "PB_count", "ONT_count", "PlatFS_count",
                       "cDNA_count", "dRNA_count", "CapTrap_count", "R2C2_count", "LibFS_count",
                       "LO_count", "LS_count", "DataCatFS_count")]

bed_file <- paste(c(name_bed, ".bed"), collapse="")
write.table(final_bed, file=bed_file,sep='\t', quote = FALSE, col.names = FALSE, row.names = FALSE)
