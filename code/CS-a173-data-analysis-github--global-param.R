{ # global  parameters (for all datasets) ----
  # short patent codes
  patcodes <- tibble(pat     = c("D", "E", "G", "H"),
                     patcode = c("001004IG", "015005IG", "036001IG", "044012IG") )
  
  # data file meta info about WP, POS, WPPOS, frac, pat_tp, sample, method, filename, and runid
  index_file     = "data/CS-2170-prep-swdrv-of-1D-2D-runs-211231OLc1.index.txt"
  index <- read_tsv(index_file, col_types = "cccdccccc") %>%
    mutate(method = str_sub(sample, -10),
           file   = gsub(".raw", "", filename, fixed=TRUE ),
           pat    = str_sub(pat_tp, 1,1),
           tp     = str_sub(pat_tp, -1)) %>%
    left_join(patcodes, by = "pat") %>%
    select(-filename, -runid, -sample, -WPPOS) %>%
    rename(sample = pat_tp) %>%
    relocate(file, sample, pat, patcode, tp )

  # types of HAPMAPS, with rank
  plasma_protein_anno_tbl.rank <- tibble(type = c("H/MAPP", "IG (by name)"), rank = c(1,2) )
  plasma_protein_anno_tbl.file = "data/plasma-protein-anno-table--HMAP-IG.v220121.txt"
  HAPMAP2 <- read_tsv(plasma_protein_anno_tbl.file, col_types = "cccc")
  
  # Accnr numbers of Immunoglobulins (for annotation only)
  accnr.immuno <- read_tsv("data/Ig-anno-table.txt", col_types = "cc")
  
  remove(index, plasma_protein_anno_tbl.file)
}
