# IMETAGENE HELPER.R
library(Rsamtools)
library(GenomicRanges)

ismetagene <- function(loaded_mg) {
  ret = 0
  
  if(!"metagene" %in% class(loaded_mg)){
    ret = 1
  } else {
    if(is.null(loaded_mg$flip_regions)) {  # old version of metagene
      ret = 2
    }
  }
  
  
  return(ret)
}

getBEDlevels <- function(beds) {
  
  regions <- c()
  
  for(bed in beds){
    df <- read.table(file = bed, header = FALSE)
    regions <- c(regions, levels(df[,1]))
  }
  
  unique(regions)
}

getBAMlevels <- function(bam) {
  names(scanBamHeader(bam)[[1]]$targets)
}

getBEDregionSize <- function(beds) {
  min_size <- NULL
  max_size <- NULL
  
  for(bed in beds){
    df <- read.table(file = bed, header = FALSE)
    tmp_max <- max((df[,3] - df[,2]))
    tmp_min <- min((df[,3] - df[,2]))
    
    
    if(is.null(min_size)) {
      min_size <- tmp_min
    } else {
      if(min_size > tmp_min) {
        min_size <- tmp_min
      }
    }
    
    if(is.null(max_size)) {
      max_size <- tmp_max
    } else {
      if(max_size < tmp_max) {
        max_size <- tmp_max
      }
    }
    
  }
  
  sizes <- list(min = min_size,max = max_size)
  return(sizes)
}


resizeRegions <- function(bed_files, new_size) {
  regions <- c()
  
  for(bed in bed_files){
    df <- read.table(file = bed, header = FALSE)
    regions <- c(regions, 
                 GRanges(seqnames = df[,1], ranges = IRanges(df[,2],df[,3])))
  }
  
  for(i in seq(1,length(regions))) {
    regions[[i]] <- resize(regions[[i]], width = new_size, fix = "center")
  }
  
  regions
}


extract_file_name <- function(x){
  fn <- strsplit(x = x , split = .Platform$file.sep)[[1]]
  fn[length(fn)]
}






