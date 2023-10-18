{# DATA LOADING  -----
  if(ANALYSE_FAIMS == 1){ # Loading scan info (FAIMS CVs) ----
      
      scan <- read_tsv( scan_tbl_file, col_types = "ccccc") %>%
        filter(msLevel == "ms2" & rt != "rt") %>%
        mutate(scan = as.numeric(scan),
               rt = as.numeric(rt), 
               cv = as.numeric(cv)) %>%
        rename(file = filename)
  
  }
  
  { # .. denovo results (ver: local conf)----
    denovo.CACHE   = gsub("#TBL", "denovo-lc", gsub("#SET", data_set_name, temptbl.CACHE))
    
    if(file.exists(denovo.CACHE)){
      denovo <- read_tsv(denovo.CACHE, col_types = "cdccdddddddddcccccccdccdddd") # Col_types until 28. Marchl "cdccdddddddddccccclldccdddd", thus PW and POS where log, now char
                                                     
    } else {
      denovo.files <- str_c(peaks_export, "~de-novo-peptides.csv")
      tag.files    <- str_c(peaks_export, "~de-novo-peptides.INQTAGS80.index")
      
      lctags.0A <- read_tsv_multi_coltypes(tag.files, coltypes = "cdcccd", append_file_name = TRUE) %>%
        mutate( localconf = as.numeric(gsub(".*de-novo-peptides\\.INQTAGS(..)\\.index", "\\1", source, perl = TRUE))) %>%
        select(-source) %>%
        group_by(file, scan, modpepseq, pepseq, INQ, localconf, tagseq)

      lcatags.0B <- lctags.0A %>%
        mutate(localconf = str_c("LC", localconf)) %>%
        summarize(.groups = "drop_last") %>%
        mutate(taglen = str_length(tagseq)) %>%
        summarize(taglen = sum(taglen), .groups = "drop") %>%
        spread(localconf, taglen)
        lctags.0 <- lctags.0A %>%
          filter(localconf == localconf.thresh) %>%
          left_join(lcatags.0B, by = c("file", "scan", "modpepseq", "pepseq", "INQ")) %>%
          ungroup()
      

      denovo <- read_denovo_multi(denovo.files) %>%
        mutate( file = gsub("\\.(.*)$", "", `Source File`, perl = TRUE),
                modpepseq = gsub("M(+15.99)", "m", Peptide, fixed = TRUE),
                modpepseq = gsub("C(+45.99)", "c", modpepseq, fixed = TRUE),
                modpepseq = gsub("N(+.98)", "n", modpepseq, fixed = TRUE),
                modpepseq = gsub("Q(+.98)", "q", modpepseq, fixed = TRUE),
                pepseq    = str_to_upper(modpepseq),
                Scan      = if_else(grepl("^F", Scan), Scan, str_c("F0:",Scan)) ) %>%
        separate(Scan, c("Fx", "scan"), sep = ":") %>%
        select(  file, scan, pepseq, modpepseq, 
               tag_length   = `Tag Length`,
               denovo_score = `Denovo Score`,
               ALC          = `ALC (%)`, 
               length,
               mz           = `m/z`,
               z, RT, Area, ppm,
               localconf    = `local confidence (%)`) %>%
        left_join(index, by = "file") %>%
        mutate(scan = as.numeric(scan)) %>%
        left_join(lctags.0 %>% select(-localconf), by = c("file", "scan", "modpepseq", "pepseq"))

      if(ANALYSE_FAIMS == 1){
        denovo <- denovo %>%
          left_join(scan %>% select(file, scan, cv), by = c("file","scan")) %>%
          mutate(frac = cv) %>% select(-cv)
      } 
      
      denovo %>% write_tsv(denovo.CACHE)  
      
      remove(tag.files, lctags.0, lctags.0A, denovo.files)
    }
    remove(denovo.CACHE) 

  }
   
  { # .. peaks search results ----
    peaksdb.CACHE  = gsub("#TBL", "peaksdb", gsub("#SET", data_set_name, temptbl.CACHE))
    
    if(file.exists(peaksdb.CACHE)){
      peaksdb <- read_tsv(peaksdb.CACHE, col_types = "cdccdcddddddccccccdc")
    } else {
      peaksdb.files <- str_c(peaks_export, "~DB-search-psm.csv")
      peaksdb.0 <- read_csv_multi(peaksdb.files) 
      peaksdb.1 <- peaksdb.0 %>%
         mutate( file = gsub("\\.(.*)$", "", `Source File`, perl = TRUE),
                modpepseq = gsub("M(+15.99)", "m", Peptide, fixed = TRUE),
                modpepseq = gsub("C(+45.99)", "c", modpepseq, fixed = TRUE),
                modpepseq = gsub("N(+.98)", "n", modpepseq, fixed = TRUE),
                modpepseq = gsub("Q(+.98)", "q", modpepseq, fixed = TRUE),
                pepseq    = str_to_upper(modpepseq),
                pep_len   = str_length(pepseq), 
                Scan      = if_else(grepl("^F", Scan), Scan, str_c("F0:",Scan)) ) %>%
        separate(Scan, c("Fx", "scan"), sep = ":") %>%
        select(  file, scan, pepseq, modpepseq, pep_len,
                 accnr        = `Accession`,
                 peaks_score  = `-10lgP`,
                 mz           = `m/z`,
                 z = Z,
                 RT, Area, ppm) %>%
        filter(!is.na(accnr)) %>% # e.g. non-human background from crap db (trypson)
        left_join(index, by = "file")
      
      if(ANALYSE_FAIMS == 1){
        peaksdb <- peaksdb.1 %>%
          mutate(scan = as.numeric(scan)) %>%
          left_join(scan %>% select(file, scan, cv), by = c("file","scan")) %>%
          #filter(is.na(cv)) # -> check here if NO LINES 
          mutate(frac = cv) %>% select(-cv)
      } else {
        peaksdb <- peaksdb.1
      }
      
      peaksdb %>% write_tsv(peaksdb.CACHE)
      remove(peaksdb.0, peaksdb.1)
    }
    remove(peaksdb.CACHE)
  } 
  
  { # .. BLAST UPSP & mprot results (same as db search) ----
    BLAST.CACHE = gsub("#TBL", "BLAST", gsub("#SET", data_set_name, temptbl.CACHE))
    if(file.exists(BLAST.CACHE)){
      blast <- read_tsv(BLAST.CACHE, col_types = "ccddddddddddcccddcc")
                                                   
    } else {
      blast.files <- str_c(peaks_export, "~de_novo_peptides_INQ.blast~CS-2139_upsp-and-mprot_v211005_I2L.txt")
      blast.qsid <- read_tsv_multi_coltypes(str_c(peaks_export, "~de_novo_peptides_INQ.txt"), coltypes = "ccdd")
      blast.0 <- read_blast_multi(blast.files) 

      blast <- blast.0 %>%
        group_by(qseqid) %>%
        mutate(rank_bitscr = rank(-bitscore, ties.method = "min")) %>%
        filter(rank_bitscr == 1) %>%
        ungroup() %>% select(-rank_bitscr) %>%
        separate(sseqid, c("undef","accnr", "ID"), sep="\\|", fill = "right") %>%
        filter(!grepl("^S0", accnr)) %>%
        select(-undef, -ID) %>%
        left_join(blast.qsid %>% rename(qseqid = accnr), by = "qseqid") %>%
        mutate(patcode = str_sub(accnr, 1,8)) %>% 
        left_join(patcodes, by = "patcode") 

      
      blast %>% write_tsv(BLAST.CACHE)
      remove(blast.0)
      
    }
    remove(BLAST.CACHE)
  }
  
  
  { # .. IMGT VDJC BLAST INQ & anno results using TAGs (ver 25. May 2022) ----

    BLAST.CACHE = gsub("#TBL", str_c("IMGT-VDJC-BLASTINQ", localconf.thresh), gsub("#SET", data_set_name, temptbl.CACHE))

    if(file.exists(BLAST.CACHE)){
      blastINQ.lc <- read_tsv(BLAST.CACHE, col_types = "cccccddddddddddccccdcccdccddccdddddddddcdc")
    } else {
      blast.files <- str_c(peaks_export, "~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.txt")
      blastINQ.qsid <- read_tsv_multi_coltypes(str_c(peaks_export, "~de_novo_peptides_INQ.txt"), coltypes = "ccdd")
      blast.0 <- read_blast_multi(blast.files) 

      blastINQ.anno <- read_tsv_multi_coltypes(str_c(peaks_export, "~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.anno80.txt"), coltypes = "ccddddddddddccccdcccdccdd") %>%
        mutate(sseqid = gsub("^..\\|","", sseqid, perl = TRUE)) %>%
        mutate(nsep = str_count(sseqid, "\\|")) %>%
        mutate(sseqid = if_else(nsep == 1, str_c(sseqid,"|---"), sseqid)) %>%
        separate(sseqid, sep = "\\|", into = c("accnr", "IMGT_gene", "IMGT_domain"), remove = FALSE) %>%
        select(-nsep) %>%
        left_join(blastINQ.qsid %>% select(qseqid = accnr, pepINQ = pepseq), by = "qseqid")
      
      blastINQ.vreg <- read_tsv_multi_coltypes(str_c(peaks_export, "~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.vreg80.txt"), coltypes = "cccdddddddddccdc") %>%
        rename(tag_status_reg = tag_status)

      blastINQ.lc <- blastINQ.anno %>%
        group_by(qseqid) %>%
        mutate(rank_bitscr = rank(-bitscore, ties.method = "min")) %>%
        filter(rank_bitscr == 1) %>%
        ungroup() %>% select(-rank_bitscr) %>%
        left_join(blastINQ.vreg, by = c("qseqid", "sseqid", "pepseq"))
      
      blastINQ.lc %>% write_tsv(BLAST.CACHE)
      remove(blastINQ.anno, blastINQ.vreg, blastINQ.qsid)
    }
    remove(BLAST.CACHE)
  }
}
