setwd("D:/local/DATA/a/a173/github/")

{ # libraries and functions ----
  library(tidyverse)
  library(tictoc)
  source("code/CS-a173-data-analysis-github--fun.R")
}


{ # parameters & sources -----
  { # dataset specific parameters

    # variable dataset: used for filtering and anotations IN plots
    dataset.ALL            = c("1D", "FAIMS", "2D")
    
    # dataset ID: string will be part (prefix) of each output file (plot and data)
    data_set_name.ALL <- c("CS-a173-1D", "CS-a173-FAIMS" , "CS-a173-2D" )
    peaks_export.ALL <- c("data/1D-all", "data/FAIMS-all", "data/2D-pat1;data/2D-pat2;data/2D-pat3;data/2D-pat4")    

    # CV fraction from scan table: just required for FAIMS dataset
    scan_tbl_file.ALL       <- c("",
                                 "data/FAIMS-all~220204OLc1_CS-2170-td210712-Xy-AP20x--DDAft90-FAIMS.scantbl.txt",
                                 "")
    ANALYSE_FAIMS.ALL      = c(0,1,0)
  
    DO_PLOTTING = 1 # set to 1 for  plot generation
    plot.outdir = "results/"
    pdat.outdir = "results/"
    
  }
  { # filters and thresholds
    localconf.thresh = 80
    
    filter.BLAST_evalue = 0.1
    filter.BLAST_length = 6
    filter.seq_simil    = 90
    
    filter.min_vreg_cover = 3
  }
  
  { # load global parameters & settings
    source("code/CS-a173-data-analysis-github--global-param.R")
  }
}

#  _                      __        _ _      _      _                _
# | |___  ___ _ __   ___ / _|  __ _| | |  __| |__ _| |_ __ _ ___ ___| |_ ___
# | / _ \/ _ \ '_ \ / _ \  _| / _` | | | / _` / _` |  _/ _` (_-</ -_)  _(_-<
# |_\___/\___/ .__/ \___/_|   \__,_|_|_| \__,_\__,_|\__\__,_/__/\___|\__/__/
#            |_|
# START LOOP THROUGH 3 DATA SETS -----
i.dataset  = 1
refresh_combined   = 1
refresh_opt_denovo = 1
refresh_denovo_cover = 1

curr_vDate = vDate()

for(i.dataset in 1:3){ # starting main loop -----
  {
    peaks_export = unlist(str_split(peaks_export.ALL[i.dataset], ";"))
    data_set_name      = data_set_name.ALL[i.dataset]
    dataset            = dataset.ALL[i.dataset]
    scan_tbl_file      = scan_tbl_file.ALL[i.dataset]
    ANALYSE_FAIMS      = ANALYSE_FAIMS.ALL[i.dataset]
  
    plot_base_file = gsub("#DATE", curr_vDate, str_c(plot.outdir, data_set_name, "--#PLOT.#DATE.png") )
    pdat_base_file = gsub("#DATE", curr_vDate, str_c(pdat.outdir, data_set_name, "--#PLOT.#DATE.txt") )
    temptbl.CACHE = gsub("#SET", data_set_name, "data/#SET.#TBL.CACHE")
    
    cat(str_c("# STARTING DATA SET ", dataset, "\n"))
  }  

  { ## DATA LOADING ----
    
    source("code/CS-a173-data-analysis-github--loading.R")

    # select primary hits only
    denovo  <- denovo  %>% group_by(scan, file) %>% mutate(score_rank = rank(-denovo_score, ties.method = "min")) %>% filter(score_rank == 1) %>% select(-score_rank)
    peaksdb <- peaksdb %>% group_by(scan, file) %>% mutate(score_rank = rank(-peaks_score, ties.method = "min")) %>% filter(score_rank == 1) %>% select(-score_rank)
    
  }
  
  { ## COMBINING DATA of db+dn+blast ----
    source ("code/CS-a173-data-analysis-github--combine.R")
  }
  
  { # FILTER DATA
    source("code/CS-a173-data-analysis-github--filter.R")
  }
  
  #       _     _   _   _
  #  _ __| |___| |_| |_(_)_ _  __ _
  # | '_ \ / _ \  _|  _| | ' \/ _` |
  # | .__/_\___/\__|\__|_|_||_\__, |
  # |_|                       |___/
  ## PLOTTING AND ANALYZING --------

  if(DO_PLOTTING == 1){  

    { ### PEAKSDB  ----
      
      { #### Fig2: distributions of peptide length and charge ----
        
        curr.palette <- c("#5dd39e", "#348aa7", "#525174", "#513b56", "#ffb997", "#f67e7d", "#843b62")

        peplen_breaks <- c(0,11,18,100)
        peplen_labels.0 <- c("short (<12)", "middle (12-18)", "long (>18)")
        peplen_labels <- tibble(peplen.bins      = factor(c("(0,11]", "(11,18]", "(18,100]")),
                                `Peptide length` = factor(peplen_labels.0, levels = peplen_labels.0) ) 
          
        dat.peaksdb_pepids.peplen <- peaksdb %>% 
          mutate(peplen.bins = cut(pep_len, breaks = peplen_breaks)) %>%
          left_join(peplen_labels, by = "peplen.bins") %>% select(-peplen.bins)
            
        dat.peaksdb_pepids.charge <- peaksdb %>% 
          filter(z > 1) %>%
          mutate(`Peptide charge` = if_else(z > 4, str_c("\u2265","5+"), str_c(z,"+")),
                 `Peptide charge` = factor(`Peptide charge`, levels = c("2+", "3+", "4+", str_c("\u2265","5+")))) 
        
        pdat.peaksdb_pepids.all <- bind_rows( 
          dat.peaksdb_pepids.charge %>%
            mutate(scanid = str_c(file, ".", scan),
                   parameter = "Peptide charge") %>%
            group_by(parameter, sample, pat, tp, frac, bin = `Peptide charge`) %>%
            summarize(n_PSM = n_distinct(scanid), n_peptides = n_distinct(pepseq), .groups = "drop"),
          dat.peaksdb_pepids.charge %>%
            mutate(scanid = str_c(file, ".", scan),
                   parameter = "Peptide charge",
                   sample = "ALL", pat = "ALL", tp = "ALL") %>%
            group_by(parameter, sample, pat, tp, frac, bin = `Peptide charge`) %>%
            summarize(n_PSM = n_distinct(scanid), n_peptides = n_distinct(pepseq), .groups = "drop"),
          # peptide lengths
          dat.peaksdb_pepids.peplen %>%
            mutate(scanid = str_c(file, ".", scan),
                   parameter = "Peptide length") %>%
            group_by(parameter, sample, pat, tp, frac, bin = `Peptide length`) %>%
            summarize(n_PSM = n_distinct(scanid), n_peptides = n_distinct(pepseq), .groups = "drop"),
          dat.peaksdb_pepids.peplen %>%
            mutate(scanid = str_c(file, ".", scan),
                   parameter = "Peptide length",
                   sample = "ALL", pat = "ALL", tp = "ALL") %>%
            group_by(parameter, sample, pat, tp, frac, bin = `Peptide length`) %>%
            summarize(n_PSM = n_distinct(scanid), n_peptides = n_distinct(pepseq), .groups = "drop")
        ) %>% mutate( frac = if_else(is.na(frac), 0, frac),
                      frac = as.factor(frac))
        
        counts = "n_PSM"
        values = "absolute"
          
        plot_name = str_c(values," ",counts)
        file_tag  = str_c("peaks-", str_sub(counts,3,5),"-",str_sub(values,1,3), "-cnts")

        nfrac = pdat.peaksdb_pepids.all %>% pull(frac) %>% unique() %>% length()
                
        plot.peaksdb_pepids.all <- pdat.peaksdb_pepids.all %>%
          mutate(use = if_else(nfrac > 5 & sample != "ALL", 0, 1)) %>%
          filter(use == 1) %>%
          gather(n_PSM, n_peptides, key = "count", value = "n") %>%
          filter(count == !!counts) %>%
          ggplot(aes(x = sample, y = n, fill = bin))
        
        if(values == "relative"){
          plot.peaksdb_pepids.all <- plot.peaksdb_pepids.all +
            geom_bar(stat="identity", position = "fill") + 
            scale_y_continuous(labels = scales::percent_format()) 
        } else {
          plot.peaksdb_pepids.all <- plot.peaksdb_pepids.all +
            geom_bar(stat="identity") 
        }
                
        plot.peaksdb_pepids.all <- plot.peaksdb_pepids.all +
          facet_grid(parameter~frac, scales = "free") +
          scale_fill_manual(values = curr.palette) + 
          ylab(counts) 
        
        if(nfrac > 5){
          plot.width = 10
        } else {
          plot.width = 1.5 + (nfrac)*1.75  
        }
        
        plot.peaksdb_pepids.all %>%
          ggsave( filename = gsub("#PLOT", file_tag, plot_base_file),
                  width = plot.width, height = 6, scale = 1.25)
        
        pdat.peaksdb_pepids.all %>%
          mutate(dataset = !!dataset) %>%
          write_tsv(gsub("#PLOT", "peaks-pep-n-PSM-stats", pdat_base_file))

        remove(dat.peaksdb_pepids.peplen)
        remove(dat.peaksdb_pepids.charge)
        remove(pdat.peaksdb_pepids.all)
        remove(plot.peaksdb_pepids.all)
      }

      { #### Fig2: box plots of n unique peptides by precursor type ----
        peakdsdb_HAPMAP_anno <- peaksdb %>%
          left_join(
            peaksdb %>%
              ungroup() %>%
              select(pepseq, accnr) %>%
              unique() %>%
              mutate(accnrs = accnr) %>%
              separate_rows(accnr, sep = ":") %>%
              inner_join(HAPMAP2, by = "accnr") %>%
              group_by(accnrs, pepseq) %>%
              summarize(gene = paste0(unique(gene), collapse = ","),
                        type = paste0(unique(type), collapse = ","),
                        .groups = "drop") %>% 
              rename(accnr = accnrs),
            by = c("accnr", "pepseq") ) %>% 
          mutate(type = if_else(is.na(type), "other", type),
                 type = gsub("(by name)", "", type, fixed = TRUE) )
        
        data.peaksdb_pepids3 <- bind_rows( peakdsdb_HAPMAP_anno,
                                           peaksdb %>% mutate(type = "All") ) %>%
          group_by(type, sample) %>% summarize(n_PSM = length(scan),
                                               n_pep = n_distinct(pepseq),
                                               .groups = "drop") %>%
          mutate(method = dataset)
        
        plot.peaksdb_pepids3 <- data.peaksdb_pepids3 %>%
          rename(precursor = type) %>%
          ggplot(aes(x = precursor, y = n_pep, fill = precursor)) + 
          geom_boxplot(outlier.shape=NA) + 
          geom_jitter(aes(fill = precursor), width = 0.1, shape = 21, size = 2) +
          scale_fill_manual(values = c("#50514f", "#f25f5c", "#ffe066", "#247ba0") ) +  # #70c1b3
          facet_grid(.~method) +
          theme_light() + 
          xlab("protein precursor class") + 
          ylab("number of unique peptide sequences")
        
        plot.peaksdb_pepids3 %>%
          ggsave( filename = gsub("#PLOT", "peaks-pepids-byprec", plot_base_file),
                  width = 4, height = 3)
      }
      
    }
    
    { ### DENOVO ----
      
      source("code/CS-a173-data-analysis-github--prep-dn-pdat.R")
      
      { #### Fig3: Distribution of de novo n PSMs assigned to Ab regions (gene) and variable region ----
        
        combined.dn_psm <- bind_rows(combined.A %>% group_by(file, scan, sample, frac) %>% summarize(.groups = "drop"),
                                     combined.B %>% group_by(file, scan, sample, frac) %>% summarize(.groups = "drop") ) %>%
          group_by(file, scan, sample, frac) %>% summarize(.groups = "drop")
        
        pdat.dnB_n_PSM_by_igclass_frac.0 <- combined.dn_psm %>%
          left_join(combined.B %>% 
                      select(sample, file, scan, frac, dnB_pepseq, igclasses4, n_vreg, vregs),
                    by = c("file", "scan", "sample", "frac") ) %>%
          mutate(IG = if_else(is.na(igclasses4),0,1),
                 VR = if_else(is.na(vregs),0,1),
                 PSM = 1) %>%
          group_by(file, scan, sample, frac) %>%
          summarize(dn_PSM = max(PSM), IG = max(IG), VR = max(VR), .groups = "drop") %>%
          group_by(frac) %>%
          summarize(dn_hit = sum(dn_PSM) , ig_hit = sum(IG) , vreg_hit = sum(VR), .groups = "drop") %>%
          mutate(`non Ig` = dn_hit - ig_hit,
                 `Ig CR` = ig_hit - vreg_hit) %>%
          rename(`Ig VR` = vreg_hit) %>%
          select(frac, `non Ig`, `Ig CR`, `Ig VR`) %>%
          gather(-frac, key = "type", value = "n_PSM") %>%
          mutate(type = factor(type, levels = c("non Ig", "Ig CR", "Ig VR")))
        
        pdat.dnB_n_PSM_by_igclass_frac.all <- pdat.dnB_n_PSM_by_igclass_frac.0  %>%
          group_by(frac) %>% summarize(n_all = sum(n_PSM), .groups = "drop")
        
        pdat.dnB_n_PSM_by_igclass_frac_dnB <- bind_rows(
          pdat.dnB_n_PSM_by_igclass_frac.0 %>% mutate(count = "counts"),
          pdat.dnB_n_PSM_by_igclass_frac.0 %>% 
            left_join(pdat.dnB_n_PSM_by_igclass_frac.all, by = "frac") %>%
            mutate(count = "relative", n_PSM = n_PSM*100/n_all) %>% select(-n_all) ) %>%
          mutate(type = factor(type, levels = c("non Ig", "Ig CR", "Ig VR")),
                 frac = factor(frac),
                 dataset = !!dataset)
        
        plot.dnB_n_PSM_by_igclass_frac_dnB <- pdat.dnB_n_PSM_by_igclass_frac_dnB %>%
          ggplot(aes(x = frac, y = n_PSM, fill = type)) + 
          geom_bar(stat = "identity", color = "black", size = 0.5, alpha = 1) + 
          scale_fill_manual(values = rev(c("yellow", "#ff8c61", "#68b0ab", "#c8d5b9"))  ) +
          facet_grid(count~., scales = "free") +
          theme_light() + 
          xlab("protein precursor class") + 
          ylab("number/proportion of peptide spectra matches (PSM)")
      
        plot.dnB_n_PSM_by_igclass_frac_dnB %>%
          ggsave( filename = gsub("#PLOT", "denovo-PSMs-byprec-and-frac-dnB", plot_base_file),
                  width = 10, height = 6)
        
        pdat.dnB_n_PSM_by_igclass_frac_dnB %>%
          mutate(dataset = !!dataset) %>%
          write_tsv(gsub("#PLOT", "denovo-PSMs-byprec-and-frac-dnB", pdat_base_file))
      }
      

      { #### Fig4: box-plot of n(PSM) by variable region ----

        pdat.dnB_n_PSM_by_vreg <- pdat.dnB_vreg_cov %>%
          filter(!is.na(vreg)) %>%
          filter(aa_reg_covered >= !!filter.min_vreg_cover ) %>%
          # group_by(sample, file, vreg) %>%
          group_by(sample, vreg) %>%
          summarize(n_PSM = length(scan), .groups = "drop")
        
        data.denovo_vreg_pepseq.0 <- pdat.dnB_vreg_cov  %>%
          filter(!is.na(vreg)) %>%
          filter(aa_reg_covered >= !!filter.min_vreg_cover ) %>%
          group_by(sample, dnB_pepseq, n_vreg) %>%
          summarize( vregs     = paste0(vreg, collapse = "-"),
                     vregs_len = paste0(aa_reg_covered, collapse = "+"),
                     dnB_score = max(dnB_score),
                     dnB_ALC   = max(dnB_ALC),
                     igclasses = paste0(unique(igclasses4), collapse = ","),
                     IMGT_dom  = paste0(unique(IMGT_domains), collapse = ","),
                     .groups = "drop") %>%
          rowwise() 

        plot.dnB_n_PSM_by_vreg <- pdat.dnB_n_PSM_by_vreg %>%
          ggplot(aes(x = vreg, y = n_PSM, fill = vreg)) + 
          geom_boxplot(outlier.shape=NA) + 
          geom_jitter(aes(fill = vreg), width = 0.1, shape = 21, size = 2) +
          scale_fill_manual(values = c("grey70", "#6b5c7b","grey70", "#c06c84", "grey70", "#f2727f", "#f2727f") ) + 
          xlab(str_c("Variable region covered by PSM (≥ ",filter.min_vreg_cover," AAs)" ) ) + 
          ylab("Number of PSMs (per sample)")
        
        plot.dnB_n_PSM_by_vreg %>%
          ggsave( filename = gsub("#PLOT", "dnB-n-PSM-by-vreg", plot_base_file),
                  width = 4, height = 5, scale = 1.25)
        
        pdat.dnB_n_PSM_by_vreg %>%
          mutate(dataset = !!dataset) %>%
          write_tsv(gsub("#PLOT", "dnB-n-PSM-by-vreg", pdat_base_file)) 
      } 
        
      { #### Fig4: box-plots of n(PSM) by Ig class and chain ----
        
        pdat.dnB_n_PSM_by_igclass <- combined.B %>%
          select(sample, file, scan, igclasses4) %>%
          separate_rows(igclasses4, sep = ",") %>%
          mutate(igchain = str_sub(igclasses4, 3, 3)) %>%
          group_by(sample, igclasses4, igchain) %>%
          summarize(n_PSM = length(scan), .groups = "drop") 
        
        plot.dnB_n_PSM_by_igclass <- pdat.dnB_n_PSM_by_igclass %>%
          ggplot(aes(x = igclasses4, y = n_PSM, fill = igchain)) + 
          geom_boxplot(outlier.shape=NA) + 
          geom_jitter(aes(fill = igchain), width = 0.1, shape = 21, size = 2) +
          coord_flip() +
          xlab("Ig class and chain") + ylab("Number of PSM") + 
          scale_fill_manual(values = c("#E69D45", "#308695", "#D45769"))
        
        plot.dnB_n_PSM_by_igclass %>%
          ggsave( filename = gsub("#PLOT", "dnB-n-PSM-by-igclass", plot_base_file),
                  width = 4, height = 5, scale = 1.0)
        
        pdat.dnB_n_PSM_by_igclass %>%
          mutate(dataset = !!dataset) %>%
          write_tsv(gsub("#PLOT", "dnB-n-PSM-by-igclass", pdat_base_file)) 
        
      } 

      { #### Fig4: n PSM ~ VR peptides + frac ----
        
        pdat.dnB_n_PSM_by_vreg_frac <- pdat.dnB_vreg_cov %>%
          filter(!is.na(vreg)) %>%
          filter(aa_reg_covered >= !!filter.min_vreg_cover ) %>%
          # group_by(sample, file, vreg) %>%
          group_by(sample, vreg, frac) %>%
          summarize(n_PSM = length(scan), .groups = "drop")
        
        plot.dnB_n_PSM_by_vreg_frac <- pdat.dnB_n_PSM_by_vreg_frac %>%
          group_by(vreg, frac) %>%
          summarize(sumPSM = sum(n_PSM), .groups = "drop") %>%
          mutate(frac = factor(frac)) %>%
          ggplot(aes(x = frac, y = sumPSM, fill = vreg)) + 
          geom_bar(stat = "identity") + 
          scale_fill_manual(values = c("grey70", "#6b5c7b","grey70", "#c06c84", "grey70", "#f2727f", "#f2727f") ) + 
          xlab(str_c("Variable region covered by PSM (≥ ",filter.min_vreg_cover," AAs)" ) ) + 
          ylab("Number of PSMs (per sample)") + 
          facet_grid(vreg~.)
          
        plot.dnB_n_PSM_by_vreg_frac %>%
          ggsave( filename = gsub("#PLOT", "dnB-n-PSM-by-vreg-frac", plot_base_file),
                  width = 4 + 1/2*pdat.dnB_n_PSM_by_vreg_frac %>% pull(frac) %>% unique() %>% length(), height = 5, scale = 1.25)
          
        pdat.dnB_n_PSM_by_vreg_frac %>%
          mutate(dataset = !!dataset) %>%
          write_tsv(gsub("#PLOT", "dnB-n-PSM-by-vreg-frac", pdat_base_file)) 
      } 
        
    }
  }
}
  