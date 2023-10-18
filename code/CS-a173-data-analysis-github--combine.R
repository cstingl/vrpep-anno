# run first: main script until data loading
combined.CACHE   = gsub("#TBL", "combined", gsub("#SET", data_set_name, temptbl.CACHE))

if(file.exists(combined.CACHE)){
  combined <- read_tsv(combined.CACHE, col_types = "ccddccddccdcccdddccdcccddccdcdddd")
} else { # Combine PEAKSDB~DENOVO~BLAST -----

  { # HAP/MAP annotation of peaksdb
    HAPMAP.nonIg <- HAPMAP2 %>% filter(!grepl("^Immuno", Description, perl = TRUE)) 

    peaksdb.temp0 <- peaksdb %>% 
      mutate(db_INQseq = transformINQ(modpepseq)) %>%
      select(sample, file, scan, frac, z, mz, RT, db_pepseq = pepseq, db_INQseq, db_score = peaks_score, db_accnr = accnr) %>% 
      separate_rows(db_accnr, sep=":") %>%
      left_join(HAPMAP.nonIg %>% select(db_HAPMAP = gene, db_accnr = accnr), by = "db_accnr") %>%
      mutate(is_db_HAPMAP = if_else(is.na(db_HAPMAP),0,1))
    
    peaksdb_HAPMAP <- peaksdb.temp0 %>%
      group_by(sample, file, scan, frac, z, mz, RT, db_pepseq, db_INQseq, db_score) %>%
      summarize(n_db_HAPMAP     = sum(is_db_HAPMAP), .groups      = "drop") %>%
      left_join(peaksdb.temp0 %>% filter(!is.na(db_accnr)) %>%
                  group_by(sample, file, scan, frac, z, mz, RT, db_pepseq, db_INQseq, db_score) %>%
                  summarize(db_accnrs    = paste0(unique(db_accnr) , collapse = ","), .groups      = "drop"),
                by = c("sample", "file", "scan", "frac", "z", "mz", "RT", "db_pepseq", "db_INQseq", "db_score")) %>%
      left_join(peaksdb.temp0 %>% filter(!is.na(db_HAPMAP)) %>%
                  group_by(sample, file, scan, frac, z, mz, RT, db_pepseq, db_INQseq, db_score) %>%
                  summarize(db_HAPMAPs   = paste0(unique(db_HAPMAP), collapse = ","), .groups      = "drop"),
                by = c("sample", "file", "scan", "frac", "z", "mz", "RT", "db_pepseq", "db_INQseq", "db_score")) 
  }

  { # denovo+blast(UPSP)+HAPMAP
    denovo.tempA <- denovo %>% 
      mutate(dn_INQseq = transformINQ(modpepseq)) %>%
      left_join(blast %>% 
                  filter(evalue <= filter.BLAST_evalue & length >= filter.BLAST_length) %>%
                  mutate(pepseq = gsub("(.*)\\..*", "\\1", qseqid, perl = TRUE)) %>% select(pepseq, accnr, pident, length_A = length, evalue, bitscore),
                by = "pepseq") %>%
      select(sample, file, scan, frac, z, mz, RT, dn_length = length, dn_pepseq = pepseq, dn_INQseq, dn_score = denovo_score, ALC, dn_accnr = accnr, pident, l_A = length_A, evalue, bitscore) %>%
      left_join(HAPMAP.nonIg %>% select(dn_HAPMAP = gene, dn_accnr = accnr), by = "dn_accnr") %>%
      mutate(is_dn_HAPMAP = if_else(is.na(dn_HAPMAP),0,1))
    
    denovo_blast.A <- denovo.tempA %>%
      group_by(sample, file, scan, frac, z, mz, RT, dn_pepseq, dn_INQseq, dn_score, ALC) %>%
      summarize(n_dn_HAPMAP     = sum(is_dn_HAPMAP), .groups      = "drop") %>%
      left_join(denovo.tempA %>% filter(!is.na(dn_accnr)) %>%
                  group_by(sample, file, scan, frac, z, mz, RT, dn_pepseq, dn_INQseq, dn_score) %>%
                  summarize(dn_accnrs    = paste0(unique(dn_accnr) , collapse = ","), .groups      = "drop"),
                by = c("sample", "file", "scan", "frac", "z", "mz", "RT", "dn_pepseq", "dn_INQseq", "dn_score")) %>%
      left_join(denovo.tempA %>% filter(!is.na(dn_HAPMAP)) %>%
                  group_by(sample, file, scan, frac, z, mz, RT, dn_pepseq, dn_INQseq, dn_score) %>%
                  summarize(dn_HAPMAPs   = paste0(unique(dn_HAPMAP), collapse = ","), .groups      = "drop"),
                by = c("sample", "file", "scan", "frac", "z", "mz", "RT", "dn_pepseq", "dn_INQseq", "dn_score")) 
    
  }
    
  { # denovo+blast(IMGT)
    denovo.tempB <- denovo %>% 
      mutate(dn_INQseq = transformINQ(modpepseq)) %>%
      left_join(blastINQ.lc %>% 
                  filter(evalue <= filter.BLAST_evalue & length >= filter.BLAST_length) %>%
                  mutate(pepseq = gsub("(.*)\\..*", "\\1", qseqid, perl = TRUE)) %>% 
                  select(pepseq, accnr, pident, length, evalue, bitscore, IMGT_gene, IMGT_domain, vreg, blast_length = length, tagseq, pepseqCC, aa_reg_covered, aa_reg_aligned),
                # REMARK: join using pepseq and tagseq !!
                by = c("pepseq", "tagseq") ) %>%
      select(sample, file, scan, frac, IMGT_gene, IMGT_domain, vreg, z, mz, RT, dn_pepseq = pepseq, dn_INQseq, dn_score = denovo_score, ALC, dn_accnr = accnr, pident, length, evalue, bitscore, blast_length, tagseq, pepseqCC, aa_reg_covered, aa_reg_aligned) 

    
    denovo_blast.B <- denovo.tempB %>%
      filter(!is.na(IMGT_gene)) %>%
      mutate(vreg       = factor(vreg, levels = c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "CDR3j") ),
             vdom       = if_else(IMGT_domain %in% c("V-REGION", "D-REGION", "J-REGION"), 1, 0),
             igclass4   = str_sub(IMGT_gene, 1, 4),
             igclass    = gsub("(.*)\\*.*", "\\1", IMGT_gene, perl = TRUE)) %>%
      select(sample, file, scan, frac, IMGT_gene, igclass4, igclass, IMGT_domain, vdom, vreg, dn_pepseq, dn_INQseq, dn_score, ALC) %>%
      arrange(sample, file, scan, frac, dn_pepseq, dn_INQseq, vreg) %>%
      group_by(sample, file, scan, frac, dn_pepseq, dn_INQseq, dn_score, ALC) %>%
      summarize(n_igclass4    = n_distinct(igclass4, na.rm = TRUE),
                igclasses4    = paste0(unique(igclass4), collapse = ","),
                n_igclass     = n_distinct(igclass, na.rm = TRUE), 
                igclasses     = paste0(unique(igclass), collapse = ","),
                n_IMGT_domain = n_distinct(IMGT_domain, na.rm = TRUE),
                IMGT_domains  = paste0(unique(IMGT_domain), collapse = ","),
                n_vreg        = n_distinct(vreg, na.rm = TRUE),
                vregs         = paste0(unique(vreg), collapse = "-"),
                .groups       = "drop") %>%
      mutate(vregs = if_else(vregs == "NA", as.character(NA), vregs),
             vregs = gsub("-*NA-*", "", vregs, perl = TRUE))
  }
    
  { # combine all (final() -----
    # annotate Ig related hits
    peaksdb_HAPMAP.Ig <- peaksdb_HAPMAP %>%
      left_join(  peaksdb_HAPMAP %>%
                    mutate(accnr = db_accnrs) %>% select(db_accnrs, accnr) %>% separate_rows(accnr) %>%
                    left_join(accnr.immuno, by = "accnr") %>% filter(!is.na(igene)) %>% arrange(db_accnrs, igene) %>%
                    group_by(db_accnrs) %>% 
                    summarize(n_igenes = n_distinct(igene, na.rm = TRUE), igenes = paste0(unique(igene), collapse = ","), .groups = "drop"),
                  by = "db_accnrs") %>%
      mutate(n_igenes = if_else(is.na(n_igenes),as.integer(0),n_igenes))
    
    denovo_blast.A.Ig <- denovo_blast.A %>%
      left_join(  denovo_blast.A %>%
                    mutate(accnr = dn_accnrs) %>% select(dn_accnrs, accnr) %>% separate_rows(accnr) %>%
                    left_join(accnr.immuno, by = "accnr") %>% filter(!is.na(igene)) %>% arrange(dn_accnrs, igene) %>%
                    group_by(dn_accnrs) %>% 
                    summarize(n_igenes = n_distinct(igene, na.rm = TRUE), igenes = paste0(unique(igene), collapse = ","), .groups = "drop"),
                  by = "dn_accnrs") %>%
        mutate(n_igenes = if_else(is.na(n_igenes),as.integer(0),n_igenes))
    
    combined.0 <- peaksdb_HAPMAP.Ig %>%
      mutate(scan = as.numeric(scan)) %>%
      select(sample, file, scan, frac, db_pepseq, db_INQseq, db_score, n_db_HAPMAP, db_accnrs, db_HAPMAPs, db_n_igenes = n_igenes, db_igenes = igenes) %>%
      full_join(denovo_blast.A.Ig %>%
                  select(sample, file, scan, frac, dnA_pepseq = dn_pepseq, dnA_INQseq = dn_INQseq, dnA_score = dn_score, dnA_ALC = ALC, n_dn_HAPMAP, dn_accnrs, dn_HAPMAPs, dnA_n_igenes = n_igenes, dnA_igenes = igenes),
                by = c("sample", "file", "scan", "frac")) %>%
      full_join(denovo_blast.B %>%
                  select(sample, file, scan, frac, dnB_pepseq = dn_pepseq, dnB_INQseq = dn_INQseq, dnB_score = dn_score, dnB_ALC = ALC, igclasses4, IMGT_domains, n_vreg, vregs),
                by = c("sample", "file", "scan", "frac"))
    

    # calculate the similarity of denovo sequences to sequences from database searches results
    SIMIL.A <- combined.0 %>% 
      group_by(dnA_INQseq, db_INQseq) %>%
      summarize(n_scans = length(scan), .groups = "drop") %>%
      rowwise() %>%
      mutate(match_INQseq = if_else(dnA_INQseq == db_INQseq,1,0)) %>%
      rowwise() %>%
      mutate(INQ_simil = if_else((is.na(dnA_INQseq) | is.na(db_INQseq)),0,
                                 if_else(dnA_INQseq == db_INQseq, 100,
                                seqsimil(dnA_INQseq, db_INQseq))))

    match_accnr <- combined.0 %>%
      filter(!is.na(dn_accnrs) & !is.na(db_accnrs)) %>%
      select(sample, file, scan, frac, db_accnrs, dn_accnrs) %>%
      separate_rows(dn_accnrs) %>%
      separate_rows(db_accnrs) %>%
      mutate(match = if_else(dn_accnrs == db_accnrs, 1, 0)) %>%
      group_by(sample, file, scan, frac) %>%
      summarize(match_accnr = max(match), .groups = "drop") 
  }
    
  combined <- combined.0 %>%
    left_join(SIMIL.A, by = c("dnA_INQseq", "db_INQseq") ) %>%
    left_join(match_accnr, by = c("sample", "file", "scan", "frac"))
    
  combined %>% 
    write_tsv(combined.CACHE)
    #1D:  121,104 rows
    remove(HAPMAP.nonIg)
    remove(peaksdb.temp0)
    remove(peaksdb_HAPMAP)
    remove(denovo.tempA)
    remove(denovo_blast.A)
    remove(denovo.tempB)
    remove(denovo_blast.B)
    remove(peaksdb_HAPMAP.Ig)
    remove(denovo_blast.A.Ig)
    remove(combined.0)
    remove(SIMIL.A)
    remove(match_accnr)
}
