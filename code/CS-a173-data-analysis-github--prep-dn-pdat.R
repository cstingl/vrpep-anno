# Development version, see `CS-2170-dev-plotting-of-dnB-results.v220824.R`
# this script...
#  - blastINQ.lc_vregcov: linking combined & coverage of variable region 
#  - blastINQ.lc_vregcov.long: linking combined & coverage of variable region (exceptions)
#  - blastINQ.lc_vregcov.min_on_lost: linking combined & coverage of variable region (exceptions 2) 
#  - denovo.tagseq_index: linking combined & coverage of variable region (exception 3)
#  - combined.B.vrec_cov: returns coverage of variable region (by peptide and region)
#  - pdat.IMGT_regions:

denovo_cover_anno.CACHE =  gsub("#TBL", "denovo-cover-anno", temptbl.CACHE)

if(refresh_denovo_cover == 0 & file.exists(denovo_cover_anno.CACHE) ){
  pdat.dnB_vreg_cov <- read_tsv(denovo_cover_anno.CACHE, col_types = "ccddccddccdcdc") %>%
    mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j")))

} else {
    { # 3.1 blast resuls: keys = pepseq, tagseqs and vreg
        blastINQ.lc_vregcov <- blastINQ.lc %>% 
          filter(evalue <= filter.BLAST_evalue & length >= filter.BLAST_length) %>%
          select(qseqid,  tagseq, vreg, aa_reg_covered) %>%
          group_by(qseqid,  tagseq, vreg) %>%
          summarize(aa_reg_covered = median(aa_reg_covered), .groups = "drop") %>%
          mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j")),
                 dnB_pepseq  = gsub("(.*)\\..*", "\\1", qseqid, perl = TRUE),
                 tagseq = transformBOJ(tagseq, ver=2)) %>%
          select(-qseqid) %>%
          arrange(dnB_pepseq, tagseq, vreg) %>%
          group_by(dnB_pepseq, vreg) %>%
          summarize(tagseqs = paste0(unique(tagseq), collapse = "+"),
                    aa_reg_covered = median(aa_reg_covered), .groups = "drop")
    }
    { # 3.1B - fixing nested tagseqs
    
      blastINQ.lc_vregcov.long <- blastINQ.lc %>% 
        filter(evalue <= filter.BLAST_evalue & length >= filter.BLAST_length) %>%
        select(qseqid,  tagseq, vreg, aa_reg_covered) %>%
        mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j")),
               dnB_pepseq  = gsub("(.*)\\..*", "\\1", qseqid, perl = TRUE),
               tagseqs = transformBOJ(tagseq, ver=2)) %>%
        select(-qseqid) %>%
        group_by(dnB_pepseq,  tagseqs, vreg) %>%
        summarize(aa_reg_covered_B = median(aa_reg_covered), .groups = "drop")
    }
    { # 3.2C - fixing cases of (falsely) nested AND concatenated tagseqs (by offering the shortest coverage though join on pepseq+region only)
      blastINQ.lc_vregcov.min_on_lost <- blastINQ.lc %>% 
        filter(evalue <= filter.BLAST_evalue & length >= filter.BLAST_length) %>%
        select(qseqid,  tagseq, vreg, aa_reg_covered) %>%
        mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j")),
               dnB_pepseq  = gsub("(.*)\\..*", "\\1", qseqid, perl = TRUE),
               tagseqs = transformBOJ(tagseq, ver=2)) %>%
        select(-qseqid) %>%
        group_by(dnB_pepseq,  vreg) %>%
        summarize(aa_reg_covered_C = min(aa_reg_covered), .groups = "drop")
    }
    { # 3.3 make tagseqs index (from denovo table)
      denovo.tagseq_index <- denovo %>% group_by(file, scan, frac, pepseq) %>% 
        mutate(tagseq = transformBOJ(tagseq, ver=2) ) %>%
        arrange(file, scan, frac, pepseq, tagseq) %>%
        summarize(tagseqs = paste0(tagseq, collapse = "+"), .groups = "drop" )
    }
    { # 3.4 Join tagseqs index an coverage table on combined.B (long version)
      combined.B.vrec_cov <- combined.B %>% #15,380 rows
        select(-frac) %>%
        left_join(denovo.tagseq_index %>% rename(dnB_pepseq = pepseq) %>% mutate(frac = if_else(is.na(frac), 0, frac)), 
                  by = c("file", "scan", "dnB_pepseq")) %>% # 15,380 -> right
        rename(vreg = vregs) %>%
        separate_rows(vreg, sep = "-") %>% # 20,033 rows
        left_join(blastINQ.lc_vregcov, by = c("dnB_pepseq", "tagseqs", "vreg"))  %>% # 20,033
        left_join(blastINQ.lc_vregcov.long, by = c("dnB_pepseq", "tagseqs", "vreg"))  %>% # 20,033
        mutate(aa_reg_covered = if_else(is.na(aa_reg_covered) & !is.na(aa_reg_covered_B), aa_reg_covered_B, as.numeric(aa_reg_covered) )) %>%
        left_join(blastINQ.lc_vregcov.min_on_lost, by = c("dnB_pepseq", "vreg")) %>%
        mutate(aa_reg_covered = if_else(is.na(aa_reg_covered) & !is.na(aa_reg_covered_C), aa_reg_covered_C, as.numeric(aa_reg_covered) )) %>%
        select(sample, file, scan, frac, dnB_pepseq, dnB_INQseq, dnB_score, dnB_ALC, igclasses4, IMGT_domains, n_vreg, vreg, aa_reg_covered) %>%
        mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j")))
  
      pdat.dnB_vreg_cov <- combined.B.vrec_cov %>%
        mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j"))) %>%
        filter(!is.na(vreg)) 
      
      pdat.dnB_vreg_cov %>% 
        mutate(dataset = !!dataset) %>%
        write_tsv(denovo_cover_anno.CACHE)
    }
    remove(blastINQ.lc_vregcov)
    remove(blastINQ.lc_vregcov.long)
    remove(blastINQ.lc_vregcov.min_on_lost)
    remove(denovo.tagseq_index)
    remove(combined.B.vrec_cov)

}


{ # - Distribution of IMGT regions (theoretical)

  pdat.IMGT_regions <- read_tsv("external-data/IMGT-V-seq_v210614.igblast7.tidy.txt", col_types = "ccccdddddc") %>%
    separate(query, c("IMGTID", "IMGT_gene"), sep = "\\|", extra = "drop") %>%
    mutate(vreg = gsub(" \\(germline\\)", "", gsub("-IMGT", "", gsub("FR", "FWR", alignment)))) %>%
    filter(vreg != "Total") %>%
    select(IMGT_gene, vreg, from, to, length) %>%
    mutate(vreg = factor(vreg, levels = c( "FWR1", "CDR1", "FWR2", "CDR2",  "FWR3", "CDR3", "CDR3j"))) 

}


remove(denovo_cover_anno.CACHE)
