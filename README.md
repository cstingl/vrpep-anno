# Variable region peptide annotation (vrpep-anno)

Toolbox and workflow for annotating the variable antibody region of denovo sequenced peptides

The project, data acquisition and data processing are scribed in "Improved detection of tryptic immunoglobulin variable region peptides by chromatographic and gas-phase fractionation techniques." ([Stingl et al. 2023](); submitted manuscript). Acquired RAW data, processed meta and generated result data are publicly available on the proteomexchange/PRIDE[^1] repository (identified PXD046072) after manuscript is accepted (credentials for reviewers send with the manuscript submission). This document describes the key data analysis to annotated acquired spectra to antibody complementary determining (CDR) and framework regions (FWR).

## Input data from proteomexchange/PRIDE repository

The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the PRIDE [1] partner repository with the dataset identifier PXD046072 and 10.6019/PXD046072. ([PRIDE proteomics identification database](https://www.ebi.ac.uk/pride/))

## Data processing

1. download `1D-all-de-novo-tables.7z` from proteomexchange/repositor PXD046072 and extract content in a `data` subfolder

**Prepare external (other sources)**

2. download germline sequences from IMGT for all V-, D-, J- and C-Regions, (e.g. <http://www.imgt.org/genedb/GENElect?query=7.6+IGHV&species=Homo+sapiens> for IGHV) and assemble fasta files.
3. generate version of IMGT fasta with Isoleucine (Ile) substituted to the isobaric Leucine (Leu); e.g. by running `code\fasta-I2L.pl external-data\IMGT-VDJC_v210614.fasta`
  * conduct IgBLAST using the IMGT fasta file as sequence and query database. Use format option 7 (`-outfmt 7`) and species *H. sapiens* (`-organism human`). Reformat output to tab-separated with following columns: query, alignment, from, to, length, matches, .mismatches, gaps, percent_identity, subject_id. This table contains germline AA position of CDRs and FWRs. (e.g. run `code\parse-igblast7.pl external-data\IMGT-VDJC_v210614.igblast7.txt`)
  
**Prepare and run BLAST searches**

4. Prepare BLAST query files (fasta file from PEAKS `de novo peptides.csv`):

```
code\csv-to-fasta.pl data\1D-all~de-novo-peptides.csv
```

**Remark:** As a result of this step a second fasta file with the post-fix `_INQ.fasta` is generated that contains peptide sequences with following three substitutions: Ile to Leu, deamidated GLn to Glu and deamidated Asn to Asp. Furthermore an index file (`_INQ.txt`) is generated that contains a reference table (index) of the accession numbers-referring to the initial peptide sequences-and the peptide sequence versions.

5. run BLAST search using `de_novo_peptides_INQ` as query file and `IMGT-VDJC_v210614_I2L` as sequence file (database):

```
blastp -query de_novo_peptides_INQ.fasta -out data\de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.txt -db external-data\IMGT-VDJC_v210614_I2L -comp_based_stats 0 -matrix PAM30 -num_alignments 100 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq qseq"
```

The corresponding output file is also available for download from the proteomexchange/PRIDE repository [PXD046072](https://www.ebi.ac.uk/pride/).

6. Extract TAGs with local confidence > 80.

```
code\peaks-anno-localconf.pl data\1D-all~de-novo-peptides.csv
code\peaks-query-localconf.pl data\1D-all~de-novo-peptides.TAGS.txt 80
```

7. Assign the CDR and FWR regions to the *de novo* search sequences

```
code\blast-anno-igreg.pl data\1D-all~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.txt localconf=80

```

## Analysis and plotting of data

For further analysis and plotting of results we used the statistical programming language *R* (version 4.2.1)[^2] and the *Tidyverse* libraries[^3]). The R file `CS-s173-DataAnalysis-github.R` contains the corresponding code and conducts following processing steps:

### Data loading

Following tables are loaded through `CS-a173-data-analysis-github--loading.R`:

+ Just for the FAIMS dataset, a scan table that contains the FAIMS compensation voltages (CV) for each scan: `data/FAIMS-all~220204OLc1_CS-2170-td210712-Xy-AP20x--DDAft90-FAIMS.scantbl.txt`
+ denovo tag files `data/1D-all~de-novo-peptides.INQTAGS80.index` (sequence TAGs with local confidence > 80)
+ denovo `data/1D-all~de-novo-peptides.csv`
+ peaks `data/1D-all~DB-search-psm.csv`
+ BLAST query sequence index (accnr, pepseq) `data/1D-all~de_novo_peptides_INQ.txt`
+ blast.files (UP/SP) `data/1D-all~de_novo_peptides_INQ.blast~CS-2139_upsp-and-mprot_v211005_I2L.txt`
+ blast.files (IMGT) `data/1D-all~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.txt`
+ IMGT BLAST variable regions `data/1D-all~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.vreg80.txt`
+ IMGT BLAST annotations `data/1D-all~de_novo_peptides_INQ.blast~IMGT-VDJC_v210614_I2L.anno80.txt` 

Loading and processed data is stored in temporary cache files to speed-up plot follow-up data analysis and plot generation. The following cache files are stored:

+ denovo.CACHE `data/CS-a173-1D.denovo-lc.CACHE`
+ peaks.CACHE `data/CS-a173-1D.peaksdb.CACHE`
+ BLAST.CACHE (UP/SP BLAST) `data/CS-a173-1D.BLAST.CACHE`
+ BLAST.CACHE (IMGT BLAST) `data/CS-a173-1D.IMGT-VDJC-BLASTINQ80.CACHE`

### Combining DB, *de novo* and BLAST search data

The individual data tables are joined in the following combinations by `code/CS-a173-data-analysis-github--combine.R`:

+ peaksdb results+HAP/MAP annotation
+ denovo+blast(UPSP)+HAPMAP (A)
+ denovo+blast(IMGT) (B)
+ combine all (A+B)
+ calculate similarity between peaksdb and denovo/BLAST(UPSP) sequences
+ match accnr of peakDB search and denovo/BLAST(UPSP)

### Filtering *de novo* results

The goodness of the *de novo* hit is assessed on the basis of the sequence similarity to the corresponding database search hit. Comparison is conducted spectra-wise and database search hits with and FDR < 1% are assumed to be correct. A sequence similarity of at least 90% was defined as a valid *de novo* hit. Filtering of *de novo* hits through `code/CS-a173-data-analysis-github--filter.R` is conducted in the following steps:

+ Determine de novo score thresholds
+ Generates de novo score vs similarity plots
+ Applying de novo filters:
  * denovo/BLAST(UPSP) hits > threshold (combined.A)
  * denovo/BLAST(IMGT) hits > threshold (combined.B)

### Plotting

* Bar charts of number of DB search hits by fraction, peptide length and charge state (`results/*--peaks-PSM-abs-cnts.v*.png`)
* Boxplot of number of DB search peptide IDs by precursor type (`results/*--peaks-pepids-byprec.v231017.png`)
* Barchart of number of *de novo* PSMs grouped by alignment to variable region, constant region or non-Ig protein (`results/*--denovo-PSMs-byprec-and-frac-dnB.v*.png`)
* Boxplot of number of *de novo* PSMs by CDR and FWR (`results/*--dnB-n-PSM-by-vreg.v*.png`)
* Boxplot of number of *de novo* PSMsby Ig chain an type (`results/*--dnB-n-PSM-by-igclass.v*.png`)
* Barchart of number of *de novo* PSMs by fraction and variable region (`results/*--dnB-n-PSM-by-vreg-frac.v*.png`)

[^1]: Vizcaíno, J. A. et al. 2016 update of the PRIDE database and its related tools. Nucleic Acids Res. 44, D447–D456 (2016).
[^2]: R Core Team. R: A Language and Environment for Statistical Computing. (2020). <https://www.r-project.org/>
[^3]: Wickham, H. et al. Welcome to the Tidyverse. J. Open Source Softw. 4, 1686 (2019). <https://joss.theoj.org/papers/10.21105/joss.01686#>
