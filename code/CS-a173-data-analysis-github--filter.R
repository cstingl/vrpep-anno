opt_denovo.CACHE   = gsub("#TBL", "opt-denovo", gsub("#SET", data_set_name, temptbl.CACHE))

if(refresh_opt_denovo == 0 & file.exists(opt_denovo.CACHE)){
  OPT_denovo_score.0 <- read_tsv(opt_denovo.CACHE, col_types = "cdddddddd")
  OPT_denovo_score <- OPT_denovo_score.0 %>%
    mutate(peplen = factor(peplen, levels = c(OPT_denovo_score.0 %>% select(peplen, Q_len) %>% unique() %>% arrange(Q_len) %>%pull(peplen)) ) )
  
} else {
  combined.A <- combined %>%
    select(sample, file, scan, db_pepseq, db_score, dnA_pepseq, dnA_score, dnA_ALC, db_INQseq, dnA_INQseq, match_INQseq, INQ_simil, n_db_HAPMAP, db_HAPMAPs, db_n_igenes, dn_accnrs, n_dn_HAPMAP, dnA_n_igenes, igclasses4, match_accnr) %>%
    filter((n_db_HAPMAP > 0)   &!is.na(dnA_pepseq) & !is.na(dn_accnrs) ) %>%
    mutate(match_pepseq = if_else(db_pepseq == dnA_pepseq, 1,0),
           dnA_peplen   = str_length(dnA_pepseq)) %>%
    group_by(sample, file, scan, dnA_pepseq, dnA_peplen) %>%
    summarize( db_score = max(db_score),
               dnA_score = max(dnA_score),
               dnA_ALC       = max(dnA_ALC),
               match_pepseq = max(match_pepseq),
               match_INQseq = max(match_INQseq),
               match_accnr = max(match_accnr),
               INQ_simil    = max(INQ_simil),
               .groups      = "drop") %>%
    filter(!is.na(db_score)) %>%
    left_join(denovo %>% 
                group_by(sample, file, scan, frac, dnA_pepseq = pepseq) %>% 
                summarize(LC80 = max(LC80, na.rm = TRUE), .groups = "drop"), 
              by = c("sample", "file", "scan", "dnA_pepseq"))
  
    
  OPT_denovo_score.0 <- combined.A %>%
    left_join(denovo %>% select(sample, file, scan, frac, z) %>% unique(), by = c("sample", "file", "scan", "frac")) %>%
    filter( z > 1) %>%
    mutate(peplen  = cut(dnA_peplen, quantile(dnA_peplen, probs = seq(0, 1, 0.25))),
           Q_len   = cut(dnA_peplen, quantile(dnA_peplen, probs = seq(0, 1, 0.25)), labels = FALSE) ) %>%
    filter(!is.na(peplen))

  dn_score_range = OPT_denovo_score.0 %>% pull(dnA_score) %>% range()
  
  OPT_denovo_score <- tibble( peplen            = factor(as.character(), levels = levels(OPT_denovo_score.0$peplen)),
                              Q_len             = as.integer(),
                              #z                 = as.integer(),
                              n_PSM             = as.integer(),
                              mean_SIMIL        = as.numeric(),
                              mean_INQseq_match = as.numeric(),
                              mean_match_accnr  = as.numeric(),
                              Q_n_total         = as.integer(),
                              recovery          = as.integer(),
                              denovo_score      = as.integer())
  
  for (curr_len_Q in 1:4){
    curr_n_total = OPT_denovo_score.0 %>% filter(Q_len == curr_len_Q ) %>% pull(scan) %>% length()
    #curr_n_total = OPT_denovo_score.0 %>% filter(Q_len == curr_len_Q & z == curr_z) %>% pull(scan) %>% length()
    for (curr_score in seq(dn_score_range[1], dn_score_range[2],  1)){
      OPT_denovo_score.curr <- OPT_denovo_score.0 %>% 
        filter(dnA_score >= curr_score & Q_len == curr_len_Q) %>%
        group_by(peplen, Q_len) %>%
        summarize(n_PSM = length(scan), 
                  mean_SIMIL = mean(INQ_simil), 
                  mean_INQseq_match = mean(match_INQseq),
                  mean_match_accnr  = mean(match_accnr),
                  .groups = "drop") %>%
        mutate(Q_n_total = curr_n_total,
               recovery  = if_else(is.na(n_PSM), 0, n_PSM/curr_n_total),
               denovo_score = curr_score)
        
    OPT_denovo_score <- bind_rows(OPT_denovo_score, OPT_denovo_score.curr)
    }
  }

  remove(OPT_denovo_score.curr)
  
  OPT_denovo_score %>%
    write_tsv(opt_denovo.CACHE)
}

{ # Plot denovo score threshold optimization
  OPT_denovo_score.thresh <- OPT_denovo_score %>%
    filter(mean_SIMIL >= !!filter.seq_simil) %>%
    group_by(peplen) %>%
    mutate(rank_score = rank(denovo_score, ties.method = "min")) %>%
    filter(rank_score == 1)
  
  pdat.OPT_denovo_score <- OPT_denovo_score %>%
    mutate( mean_match_accnr  =  mean_match_accnr*100,
            mean_INQseq_match =  mean_INQseq_match*100,
            recovery          = recovery*100) %>%
    select(peplen, denovo_score, mean_SIMIL, mean_INQseq_match, recovery, mean_match_accnr) %>%
    gather(-peplen, -denovo_score, key = "measure", value = "proportion")
  
  plot.parameter_name <- tibble(measure   = c("mean_SIMIL", "mean_INQseq_match", "recovery", "mean_match_accnr"),
                                parameter = c("sequence (INQ) similarity", "exact sequence (INQ) match", "PSM recovery", "protein (accnr) match (BLAST)") ) %>%
    mutate(parameter = factor(parameter, levels = c("exact sequence (INQ) match", "sequence (INQ) similarity", "protein (accnr) match (BLAST)",  "PSM recovery") ) )
  
  plot.OPT_denovo_score <-  pdat.OPT_denovo_score %>%
    left_join(plot.parameter_name, by = "measure") %>%
    ggplot(aes(x = denovo_score, y = proportion, color = parameter)) + 
    geom_line(size = 1) + 
    geom_hline(yintercept = 90, linetype = "dotted", color = "darkred"  )+
    facet_grid(.~peplen) + 
    scale_y_continuous(name = "percentage", breaks = seq(0,100,by=10)) +
    xlab("denovo score (peaks)") + 
    theme(legend.position = "bottom") + 
    scale_color_manual(values = c("#0a7b83", "#2aa876", "#ffd265", "#e8554e")) + 
    ggtitle("Agreement of MS/MS database and denovo search, by denovo score and peptide length")
  
  plot.OPT_denovo_score %>%
    ggsave( filename = gsub("#PLOT", "opt-denovo-score", plot_base_file),
            width = 9, height = 6)
    
  pdat.OPT_denovo_score %>%
    mutate(dataset = !!dataset) %>%
    write_tsv(gsub("#PLOT", "opt-denovo-score", pdat_base_file))
    
  remove(plot.parameter_name)
    
}

{ # APPLYING DENOVO FILTERS
  
  OPT_denovo_score.upper <- OPT_denovo_score.thresh %>%
    mutate(`peptide length (AA)` = center_cut_value(peplen),
           upper_AA_limit = max_cut_value(peplen)) %>%
    ungroup() %>%
    select(upper_AA_limit, denovo_score)
  
  upper <- c(OPT_denovo_score.upper %>% pull(upper_AA_limit))
  score <- c(OPT_denovo_score.upper %>% pull(denovo_score))
  upper.l = 6
  
  denovo_thresh <- tibble(peplen = as.integer(), score_thresh= as.numeric())

  for(i in 1:4) {
    denovo_thresh <- bind_rows(denovo_thresh, 
                               tibble(peplen      = upper.l:upper[i],
                                      score_thresh = score[i]) )
    upper.l = upper[i] + 1
  }
  remove(upper.l, upper, score, OPT_denovo_score.upper)
  
  combined.A <- combined%>%
    filter(!is.na(dnA_pepseq)) %>%
    mutate(peplen = str_length(dnA_pepseq)) %>%
    left_join(denovo_thresh, by = "peplen") %>%
    filter(dnA_score >= score_thresh | dnA_score >= score_thresh )
  
  combined.B <- combined%>%
    filter(!is.na(dnB_pepseq)) %>%
    mutate(peplen = str_length(dnB_pepseq)) %>%
    left_join(denovo_thresh, by = "peplen") %>%
    filter(dnB_score >= score_thresh | dnB_score >= score_thresh )
  
}

{ # tidy up ...
  remove(OPT_denovo_score)
  remove(OPT_denovo_score.0)
  remove(OPT_denovo_score.thresh)
  remove(pdat.OPT_denovo_score)
  remove(plot.OPT_denovo_score)
  
}
  