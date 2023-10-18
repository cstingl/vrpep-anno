library(tidyverse)

read_tsv_multi <- function(files){
  
  for(file in files){
    
    curr <- read_tsv(file)
    
    if(file == files[1]){
      tbl_return <- curr
    } else {
      tbl_return <- bind_rows(tbl_return, curr)
    }
  }
  
  return(tbl_return)
}

read_blast_multi <- function(files){
  
  for(file in files){
    curr <- read_tsv(file, col_types = "ccdiiiiiiiddcc", na = c("N/A", "", "NA"),
                     col_names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseq", "qseq"))
    
    if(file == files[1]){
      tbl_return <- curr
    } else {
      tbl_return <- bind_rows(tbl_return, curr)
    }
  }
  
  return(tbl_return)
}


read_csv_multi <- function(files){
  
  for(file in files){
    
    curr <- read_csv(file, show_col_types = FALSE)
    
    if(file == files[1]){
      tbl_return <- curr
    } else {
      tbl_return <- bind_rows(tbl_return, curr)
    }
  }
  
  return(tbl_return)
}


read_denovo_multi <- function(files){
  
  for(file in files){
    
    curr <- read_csv(file, col_types = "icccciddididddddcccc", na = c("","-"))
    
    if(file == files[1]){
      tbl_return <- curr
    } else {
      tbl_return <- bind_rows(tbl_return, curr)
    }
  }
  
  return(tbl_return)
}


read_tsv_multi_coltypes <- function(files, coltypes, append_file_name = FALSE){
  
  for(file in files){
    
    if(append_file_name == FALSE){
      curr <- read_tsv(file, col_types = c(coltypes))
    } else {
      file_name = basename(file)
      curr <- read_tsv(file, col_types = c(coltypes)) %>% mutate(source = file_name)
    }
    
    if(file == files[1]){
      tbl_return <- curr
    } else {
      tbl_return <- bind_rows(tbl_return, curr)
    }
  }
  
  return(tbl_return)
}

vDate <- function(x) { 
  return = str_c("v",str_sub(gsub("-","",Sys.Date()),-6)) 
  return(return)
}

transformINQ <- function(seq){
  seq.INQ = gsub("I", "L", seq)
  seq.INQ = gsub("n", "D", seq.INQ)
  seq.INQ = gsub("q", "E", seq.INQ)
  return(str_to_upper(seq.INQ))
}

transformBOJ <- function(seq, ver = 2){
  seq.BOJ = gsub("[IL]", "J", seq)
  if(is.na(ver)){ver = 2}
  if( ver == 1){
    seq.BOJ = gsub("[Dn]", "B", seq.BOJ)
    seq.BOJ = gsub("[Eq]", "O", seq.BOJ)
  } else {
    seq.BOJ = gsub("[DN]", "B", str_to_upper(seq.BOJ))
    seq.BOJ = gsub("[EQ]", "O", seq.BOJ)
    
  }
  return(str_to_upper(seq.BOJ))
}


seqsimil <- function(seq1, seq2){
  l = max(str_length(seq1), str_length(seq2))
  d = adist(seq1, seq2) 
  return(round((l-d[1])/l*100,1))
}

center_cut_value <- function(x){
  S <- tibble(x = x) %>%
    mutate(min    = if_else(is.na(x), as.numeric(NA), as.numeric(gsub(".(.*)\\,.*", "\\1", x, perl = TRUE))),
           max    = if_else(is.na(x), as.numeric(NA), as.numeric(gsub(".*\\,(.*).", "\\1", x, perl = TRUE))),
           center = (min+max)/2)
  
  return(S %>% pull(center))
}

max_cut_value <- function(x){
  S <- tibble(x = x) %>%
    mutate(min    = if_else(is.na(x), as.numeric(NA), as.numeric(gsub(".(.*)\\,.*", "\\1", x, perl = TRUE))),
           max    = if_else(is.na(x), as.numeric(NA), as.numeric(gsub(".*\\,(.*).", "\\1", x, perl = TRUE))),
           center = (min+max)/2)
  
  return(S %>% pull(max))
}
