# MicrobiomeHamiltonianMedicine_repo
 
Theory model of a microbes effect on its host as a function of relatedness and transmission ecology.
Testing predictions with microbiome metagenomic data (case-control data spanning a range of diseases).
 
 
 # data

## case-control metagenomic data sources, accession, etc

- **studies_list.xslx**: liste of metagenomic study, accession number, and number of case/control samples available

- **SRAtabs_and_metadata**: directory with excel and txt tables with lists of samples and metadata for each metagenomic study

- **assembled_SRA_tables**: assembled SRA tables (by script 2_assemble_SRA_tables), ready for processing (by script 3_sample_processing)

- **Samples_accession_list.xlsx**:  complete list of each sample accession number. For processing pipeline and records.


## data to compute relatedness and transmission pattern in the microbiome

- **species_info.txt**: complete genome info about species from MIDAS database (derived from PATRIC info files). Gives exact genome name, characteristics, isolation source, etc.

- **relatedness_sporulation.txt**: relatedness and sporulation data for 101 human gut associated microbes, computed by Simonet & McNally, 2021. Publicly available, and methods descrived in Simonet & McNally 2021 (PNAS)

- **SPORULATION_SCORES_all_midas.txt**: sporulation scores computed for all organisms in the  MIDAS database (~5900 microbial organisms), following method described in brown et al 2017

## data for supplementary text and model showin negative correlation between vertical transmission and sporulation scores

- **Midas_transmission_data.xls**:  data from MIDAS paper, using MIDAS to track strain between mother-infant pairs, used to derived rates of vertical transmission. Used in supplementary material to show correlation between sporulations cores and rate of vertical transmission (phylogenetic model).

- **midas_tree_renamed.newick**: phylogeny taken from MIDAS database, with tree tips renamed to be full species name. Used for comparative analyses


 # output
 
 - **figures**: figures and table for manuscript
 
 - **metagenomes_processing**: output of script 3_sample_processing.sh and 3.2_output_sortout.sh **assembled_species_profiles.txt.bz2** is zipped file with all the metagenomes relative abundances used in the analysis. Other files are various sanity check and pipeline progress check files.
 
  - **stats_runs**: RData files. **data_preps_23062021.RData** is produced by xxxx and directly loaded in the actual statistical analyses script (4.2_statistical_analyses.R) because it takes some times to format all the data. Easier to work it out this way. AUC_crossValidation_runClean is AUC analysis. VerticalTransmission_phylogenetic_model is the supplementary text model showing negative correlation between sporulation scores and vertical transmission. MMMC_HM2_run_July2021 is not included in manuscript. Another modelling approach of the data, using multi-membership modelling with MCMCglmm.
  

# scripts

 - **Scripts 1 to 3**: data processing: assemble and format lists of sample access. Then those fed into '3_sample_processing' to download the data, process read via MOCAT, run MIDAS to estimate relative abundances. Scrip 3.2 sort out the data.
 
  - **Scripts 4.1 and 4.2**: 4.1 preps the data for statistical analyses, environment is saved in then loaded into script 4.2 for statistical analyses, and figures/tables productions.
  
   - **Scripts 5**: makes supplementary material: formats model summaries to get tex formatted tables, + supplementary model analysing vertical transmission and sporulations cores.
   
Other scripts sourced for package loading and figures formatting themes (ggplot)
 
