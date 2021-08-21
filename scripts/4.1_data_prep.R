local_project_dir="/Users/s1687811/Documents/PhD/Research/MicrobiomeHamiltonianMedicine/"
setwd(local_project_dir)
#setwd('./HamiltonianMedicineMicrobiome/')
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(vegan)

# (1) DATA PREP ----
# __ Sanity check ----

# Sanity check on reads processing
logs<- read_excel('./output/metagenome_processing/pipeline_processing.xlsx', sheet = 'hosts_with_SNPS_output') %>%
  rename(host = host_with_SNPS_output)
logs2<- logs %>% gather('step', 'nb_reads', 2:3)
head(logs2)
ggplot(logs2, aes(x = step, y = nb_reads))+geom_boxplot()
summary(logs$nb_reads_processed)
# All samples have > 1million read after processing. Keep all.

cohorts<- do.call('rbind', strsplit(logs$host, '_'))[,1]
table(cohorts)
# KAZAK and USCRC samples did not have time to run
# Keep USCRC as it is for now, has about equal amount of Healthy and Non-Healthy
# can use KAZAK samples from last year's data

logs$cohort<- cohorts
logs[logs$cohort == 'KAZAK',] # only one ran this time, delete it and replace by old ones


# __ Metdata prep & cleanup ----

md<- read.table('./data/assembled_SRA_tables.txt', header=TRUE, sep = '\t',
                na.strings = c('NA', 'N/A'),
                colClasses = c(rep('character', 6), 'character', 'character', rep('character', 3))) %>%
  select(-run, -BioSample, -BioProject) %>%
  unique()

is.numeric(md$age)
is.numeric(md$bmi)

# Recode gender
md$gender<- gsub('Female', 'female', md$gender)
md$gender<- gsub('Male', 'male', md$gender)
md$gender<- gsub('F', 'female', md$gender)
md$gender<- gsub('M', 'male', md$gender)
md$gender<- gsub('N/A', NA, md$gender, fixed = TRUE)
unique(md$gender)


# Recode country
unique(md$country)
md$country<- gsub('ESP', 'SPA', md$country)
md$country<- gsub('spanish', 'SPA', md$country)
md$country<- gsub('Kazakhstan', 'KZA', md$country)
md$country<- gsub('Chinese', 'CHN', md$country)
md$country<- gsub('China', 'CHN', md$country)
md$country<- gsub('china', 'CHN', md$country)
md$country<- gsub('Gothenburg', 'SWE', md$country)
md$country<- gsub('Germany', 'GER', md$country)
md$country<- gsub('France', 'FRA', md$country)
unique(md$country)


# Record BMI and Age and set back to numeric
unique(md$bmi)
md$bmi<- gsub('N/A', 'NA', md$bmi)
md$bmi<- gsub('>30', '31', md$bmi)
md$bmi<- gsub('<23', '22', md$bmi)
md$bmi<- as.numeric(md$bmi)

md$age<- gsub('Adult', 'NA', md$age)
md$age<- gsub('N/A', 'NA', md$age)
md$age<- as.numeric(md$age)

# Record phenotypes
unique(md$phenotype)
md$phenotype<- gsub('Ulcerative colitis', 'UC', md$phenotype)
md$phenotype<- gsub('Crohns disease', 'CD', md$phenotype)
md$phenotype<- gsub('Healthy relative', 'Healthy', md$phenotype)
md$phenotype<- gsub('hypertension', 'HYPER', md$phenotype)
md$phenotype<- gsub('pre-hypertension', 'preHYPER', md$phenotype)
md$phenotype<- gsub('cirrhosis', 'LC', md$phenotype)
md$phenotype<- gsub('Obesity', 'OBE', md$phenotype)
md$phenotype<- gsub('rheumatoid_arthritis', 'RA', md$phenotype)

# MH richness only healthy, classify BMI > 30 as obese
md[which(md$cohort == 'MHRICH' & md$bmi > 30), 'phenotype'] <- 'OBE'

# Make binary caseControl variable
md$caseControl<- ifelse(md$phenotype == 'Healthy', 1, 0)

# Filter out the healthy which are actually obese
md<- md[-which(md$phenotype == 'Healthy' & md$bmi > 30), ] # this assumes that host for which BMI info NA are not obese
length(unique(md$host))


# __ Relatedness & sporulation ----
r.df<- read.table('./data/relatedness_sporulation.txt', header = TRUE, sep = '\t')
nrow(r.df) # 101 species


# __ Abundances ----
# Load abundance table
prob.df<- read.table('./output/metagenome_processing/assembled_species_profiles.txt', header=FALSE, sep = '\t', colClasses = c('character', 'numeric', 'numeric', 'numeric', 'character'))
colnames(prob.df)<- c('species_id', 'count_reads', 'coverage', 'within_host_relative_abundance', 'host')
prob.df<-  prob.df[,c('species_id','within_host_relative_abundance', 'host')]


# Adding MetS samples from old batch
old<- read.table('./zzz_archive/CASECONTROL/data_CC_04Feb/OLD_assembled_species_profiles.txt', sep = '\t')
colnames(old)<- c('species_id', 'count_reads', 'coverage', 'within_host_relative_abundance', 'host')
old.kza<- old[which(substr(old$host, 1, 4) %in% c('713A', '713B')), ]
old.kza$host<- paste0('KAZAK_', old.kza$host, '1100')
old.kza<- old.kza %>% select(species_id, within_host_relative_abundance,  host)
head(old.kza)
head(prob.df)

# remove the one kazak sample from mew prob.df, then append old kazaks to prob.df
length(unique(prob.df$host)) # 3074 for which I had an output
prob.df<- prob.df[prob.df$host != 'KAZAK_713A0331100',]
length(unique(prob.df$host)) # 3073 now
prob.df<- rbind(prob.df, old.kza)

# Check that table is complete
length(unique(prob.df$host))  # 3074-1+86 = 3159 host for which we have an SNPS output
length(unique(prob.df$host)) * length(unique(prob.df$species_id)) # nbrows ok


# Assemble relative abundance and metadata
d<- inner_join(prob.df, md, by = 'host') # Give 3103 hosts total for which had SNP output (or had it from previous batch) AND kepts on the basis of metadata

# __ Diversity indices ----
# Compute diversity indices on WHOLE metagenome (not just focal 101 species)
d.tmp<- d %>%
  select(host, species_id, within_host_relative_abundance) %>%
  spread(species_id, within_host_relative_abundance)
rownames(d.tmp)<- d.tmp[,1]
d.tmp<- d.tmp[,-1]


dim(d.tmp) # 3103 hosts, 5952 species (So: 5952 + 1 columns)
rowSums(d.tmp[,2:5952]) # check: the relative abundances sum to 1

# Diversity indices to compare to NR

Shannon<- as.data.frame(diversity(d.tmp, index = "shannon")) %>% mutate(host = rownames(.)) %>% rename(Shannon = `diversity(d.tmp, index = "shannon")`)# Shannon (entropy)

invSimpson<- as.data.frame(diversity(d.tmp, index = "invsimpson")) %>% mutate(host = rownames(.)) %>% rename(invSimpson = `diversity(d.tmp, index = "invsimpson")`)# Simpson (eveness)

Richness<- as.data.frame(specnumber(d.tmp)) %>% mutate(host = rownames(.)) %>% rename(Richness = `specnumber(d.tmp)`)# Richness


divs.index<- left_join(Shannon, invSimpson) %>% 
  left_join(Richness) %>%
  select(host, Shannon, invSimpson, Richness)

# verified: if trim matrix to remove species present in 0 host, diversity indices are the same

head(divs.index)



# __ Select 101 ----
d.101<- d %>% filter(species_id %in% r.df$species_id)
length(unique(d.101$species_id)) # 101 species
length(unique(d.101$host))       # 3159 hosts , actually 3103 when filtering out the Healthy obese


# Check number of species found in hosts
check<- d.101 %>%
  select(host, cohort, species_id, within_host_relative_abundance, phenotype) %>%
  spread(species_id, within_host_relative_abundance)
rowSums(check[,4:104])
check$nb_sp<- rowSums(check[,4:104] > 0)
check<- check[,c(1:3, 105)]
check<- arrange(check, nb_sp)
check<- left_join(check, logs, by = 'host')
ggplot(check, aes(x = nb_reads_processed, y = nb_sp))+
  geom_point()

ggplot(check, aes(x = phenotype, y = nb_sp))+
  geom_boxplot()

# Just a couple of hosts with a few species
# CD have less species than other diseases, kind of makes sense


# METADATA OVERVIEW
md.select<- unique(d.101[,3:11]) # coding: '1' is healthy, '0' is sick
table(md.select$caseControl)          # 1250 healthy, 1853 non-healthy
table(md.select$phenotype, md.select$cohort) # Some studies have 2 disease phenotypes (e.g. T1D & T2D) or various stages (e.g. IGT & T2D, or Adenoma and CRC)



# __ Recode & filter out hosts ----



# table(md.select[md.select$cohort == 'MHMGS','phenotype'])
# table(md.select[md.select$cohort == 'iHMPIBD','phenotype'])
# table(md.select[md.select$cohort == 'CHINAIBD','phenotype'])
# unique(md.select[md.select$phenotype %in% c('T1D', 'T2D', 'IGT'),'cohort'])
# table(md.select[md.select$cohort == 'CHINAQINJ','phenotype'])
# table(md.select[md.select$cohort == 'MHIGC','phenotype'])
# table(md.select[md.select$cohort == 'SWEDISH','phenotype'])
# 
# 
# # filter out adenoma and prehyper
# d.101.trim<- d.101 %>% filter(!phenotype %in% c('adenoma', 'pre-HYPER'))
# # filter out CD from MHMGS
# d.101.trim<- d.101.trim[-which(d.101.trim$phenotype == 'CD' & d.101.trim$cohort == 'MHMGS'),]
# # filter out UC from iHMPIBD
# d.101.trim<- d.101.trim[-which(d.101.trim$phenotype == 'UC' & d.101.trim$cohort == 'iHMPIBD'),]
# # filter out T1D from MHIGC
# #d.101.trim<- d.101.trim[-which(d.101.trim$phenotype == 'T2D' & d.101.trim$cohort == 'MHIGC'),]
# # UPDATE: in fact, filter out MHIGC entirely. Its controls are all from spain and cases from denmark
# # case/control confounded with country/cohort. Unclear if it was sequenced in same batch or not.
# # drop it.
# d.101.trim<- d.101.trim[-which(d.101.trim$cohort == 'MHIGC'),]
# # filter out T2D from SWEDISH
# d.101.trim<- d.101.trim[-which(d.101.trim$phenotype == 'T2D' & d.101.trim$cohort == 'SWEDISH'),]
# 
# 
# 
# # Set disease variable
# d.101.trim$disease = NA
# 
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAJIEZ', 'SWEDACVD')]<- 'ACVD'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAYEZ')]<- 'BD'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINACRCYU', 'ITCRC', 'USCRC', 'AUTCRCFENG', 'FGCRC')]<- 'CRC'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAHYPER')]<- 'HYPER'
# d.101.trim$disease[d.101.trim$cohort %in% c('MHMGS')]<- 'UC'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAIBD', 'iHMPIBD')]<- 'CD'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINALC')]<- 'LC'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAOB', 'MHRICH')]<- 'OBE'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINARA')]<- 'RA'
# #d.101.trim$disease[d.101.trim$cohort %in% c('MHIGC')]<- 'T1D'
# d.101.trim$disease[d.101.trim$cohort %in% c('CHINAQINJ')]<- 'T2D'
# d.101.trim$disease[d.101.trim$cohort %in% c('SWEDISH')]<- 'IGT'
# d.101.trim$disease[d.101.trim$cohort %in% c('KAZAK')]<- 'MetS'
# 
# 
# md.select.trim<- unique(d.101.trim[,3:12]) # coding: '1' is healthy, '0' is sick
# nrow(md.select.trim) # 2672, i.e. : 3103 - 10 (CD from MHMGS) - 35 (UC from iHMPIBD) - 169 (ALL from MHIGC) - 51 (T2D from SWEDISH) - 110 (adenoma) - 56 (preHYPER) = 2672
# table(md.select.trim$caseControl)          # 1184 healthy, 1488 non-healthy
# table(md.select.trim$phenotype, md.select.trim$cohort) # Now all studies have a single disease
# 
# 
# # min observed number for SWEDACVD which has 10 control and 10 cases
# table(md.select.trim$disease, md.select.trim$cohort) # Now all studies have a single disease
# rowSums(table(md.select.trim$disease, md.select.trim$cohort)) # Now all studies have a single disease
# sum(table(md.select.trim$disease, md.select.trim$cohort))
# 
# 
# rowSums(table(md.select.trim$cohort, md.select.trim$disease))
# table(md.select.trim$cohort, md.select.trim$caseControl)



d.101<- left_join(d.101, divs.index)

foo<- d.101 %>%
  select(host, phenotype, cohort, caseControl, Shannon, invSimpson, Richness) %>%
  unique()

ggplot(foo, aes(x = cohort, fill = phenotype, y = Shannon))+
  geom_boxplot()




# (2) GUT RELATEDNESS & Co ----

# (2.1) Matrices ----
# Make abundance matrix: one column per species, one row per host
prob<- d.101 %>%
  select(host, species_id, within_host_relative_abundance) %>%
  spread(species_id, within_host_relative_abundance) %>%
  select(r.df$species_id) %>% # Order columns in same order as relatedness dataframe
  as.matrix()


# Cumulative relative abundance covered by these 101 species
dim(prob)
cum.101<- rowSums(prob)
hist(cum.101, xlab = 'Cumulative abundance of our 101 species') # some hosts microbiomes is poorly covered by this set of 101 species 
summary(cum.101)


# Rescale the relative abundance over that  set of 101 species
prob.rescaled<- prob/cum.101


# Multiply the rescaled relative abundance matrix by  the vectors of Relatedness, Relatedness*Sporulation
probr.rescaled<- t(t(prob.rescaled)*r.df$mean_relatedness)    # mean R
probrspo.rescaled<- t(t(prob.rescaled)*(r.df$mean_relatedness*r.df$sporulation_score)) # mean (R*SPO)


# Final data assembly
dat.rescaled<- d.101 %>%
  select(host, cohort, phenotype, caseControl, species_id, within_host_relative_abundance) %>%
  spread(species_id, within_host_relative_abundance) %>%
  select(host, cohort, phenotype, caseControl) %>%
  left_join(divs.index)


# replace cohort name by publication compatible name

dat.rescaled$cohortNew<- NA
dat.rescaled[dat.rescaled$cohort == 'CHINAJIEZ', 'cohortNew']<- 'JieZ2017_ACVD'
dat.rescaled[dat.rescaled$cohort == 'SWEDACVD', 'cohortNew']<- 'Karlsson2012_ACVD'
dat.rescaled[dat.rescaled$cohort == 'CHINAYEZ', 'cohortNew']<- 'YeZ2018_BD'
dat.rescaled[dat.rescaled$cohort == 'CHINAIBD', 'cohortNew']<- 'He2017_IBD' # CD only indeed
dat.rescaled[dat.rescaled$cohort == 'iHMPIBD', 'cohortNew']<- 'iHMP_IBD' # both CD & UC
dat.rescaled[dat.rescaled$cohort == 'AUTCRCFENG', 'cohortNew']<- 'FengQ2015_CRC'
dat.rescaled[dat.rescaled$cohort == 'ITCRC', 'cohortNew']<- 'ThomasAM2018_CRC'
dat.rescaled[dat.rescaled$cohort == 'USCRC', 'cohortNew']<- 'VogtmannE2016_CRC'
dat.rescaled[dat.rescaled$cohort == 'CHINACRCYU', 'cohortNew']<- 'YuJ2015_CRC'
dat.rescaled[dat.rescaled$cohort == 'FGCRC', 'cohortNew']<- 'ZellerG2014_CRC'
dat.rescaled[dat.rescaled$cohort == 'CHINAHYPER', 'cohortNew']<- 'LiJ2017_HYPER'
dat.rescaled[dat.rescaled$cohort == 'SWEDISH', 'cohortNew']<- 'Karlsson2013_T2D&IGT'
dat.rescaled[dat.rescaled$cohort == 'CHINALC', 'cohortNew']<- 'QinN2014_LC'
dat.rescaled[dat.rescaled$cohort == 'KAZAK', 'cohortNew']<- 'Costea2017_MetS'
dat.rescaled[dat.rescaled$cohort == 'MHRICH', 'cohortNew']<- 'LeChatelierE2013_OBE'
dat.rescaled[dat.rescaled$cohort == 'CHINAOB', 'cohortNew']<- 'Liu2017_OBE'
dat.rescaled[dat.rescaled$cohort == 'CHINARA', 'cohortNew']<- 'Zhang2015_RA'
dat.rescaled[dat.rescaled$cohort == 'MHIGC', 'cohortNew']<- 'Li2014_IGC'
dat.rescaled[dat.rescaled$cohort == 'CHINAQINJ', 'cohortNew']<- 'QinJ2012_T2D'
dat.rescaled[dat.rescaled$cohort == 'MHMGS', 'cohortNew']<- 'NielsenH2014_IBD' # 10 CD and 50 UC
unique(dat.rescaled$cohortNew)



dat.rescaled$disease_type<- NA
dat.rescaled[dat.rescaled$cohort == 'CHINAJIEZ', 'disease_type']<- 'ACVD'
dat.rescaled[dat.rescaled$cohort == 'SWEDACVD', 'disease_type']<- 'ACVD'
dat.rescaled[dat.rescaled$cohort == 'CHINAYEZ', 'disease_type']<- 'BD'
dat.rescaled[dat.rescaled$cohort == 'CHINAIBD', 'disease_type']<- 'IBD' # CD only indeed
dat.rescaled[dat.rescaled$cohort == 'iHMPIBD', 'disease_type']<- 'IBD' # both CD & UC
dat.rescaled[dat.rescaled$cohort == 'AUTCRCFENG', 'disease_type']<- 'CRC'
dat.rescaled[dat.rescaled$cohort == 'ITCRC', 'disease_type']<- 'CRC'
dat.rescaled[dat.rescaled$cohort == 'USCRC', 'disease_type']<- 'CRC'
dat.rescaled[dat.rescaled$cohort == 'CHINACRCYU', 'disease_type']<- 'CRC'
dat.rescaled[dat.rescaled$cohort == 'FGCRC', 'disease_type']<- 'CRC'
dat.rescaled[dat.rescaled$cohort == 'CHINAHYPER', 'disease_type']<- 'HYPER'
dat.rescaled[dat.rescaled$cohort == 'SWEDISH', 'disease_type']<- 'T2D_IGC' # both T2D and IGC
dat.rescaled[dat.rescaled$cohort == 'CHINALC', 'disease_type']<- 'LC'
dat.rescaled[dat.rescaled$cohort == 'KAZAK', 'disease_type']<- 'MetS'
dat.rescaled[dat.rescaled$cohort == 'MHRICH', 'disease_type']<- 'OBE'
dat.rescaled[dat.rescaled$cohort == 'CHINAOB', 'disease_type']<- 'OBE'
dat.rescaled[dat.rescaled$cohort == 'CHINARA', 'disease_type']<- 'RA'
dat.rescaled[dat.rescaled$cohort == 'MHIGC', 'disease_type']<- 'T2D_IGC' # IGC only
dat.rescaled[dat.rescaled$cohort == 'CHINAQINJ', 'disease_type']<- 'T2D_IGC' # T2D only
dat.rescaled[dat.rescaled$cohort == 'MHMGS', 'disease_type']<- 'IBD' # 10 CD and 50 UC
unique(dat.rescaled$cohortNew)




dat.rescaled<- dat.rescaled %>%
  select(host, cohortNew, cohort, phenotype, disease_type, caseControl, Shannon, invSimpson, Richness) %>%
  rename(internalName = cohort)

# Append the mean gut relatedness predictor variables
dat.rescaled$cum.101 = cum.101
dat.rescaled$nb_sp.101 = rowSums(prob.rescaled != 0)
dat.rescaled$NR = rowSums(probr.rescaled)     # mean R
dat.rescaled$NRSPO = rowSums(probrspo.rescaled) # mean (R*SPO)


# Append the relative abundance matrix (prob) and relative abundance matrix  * R (probr) to dataframe
dat.rescaled$prob = prob.rescaled
dat.rescaled$probr =  probr.rescaled
dat.rescaled$probrspo =  probrspo.rescaled



# Set data factor levels for plotting
#dat.rescaled$disease<- factor(dat.rescaled$disease,
#                              levels = c("OBE","MetS","CD","UC","CRC","IGT","T2D","LC","ACVD","RA","HYPER","BD"))

#dat.rescaled$cohort<- factor(dat.rescaled$cohort,
#                             levels = rev(c('CHINAOB', 'MHRICH', 'KAZAK','CHINAIBD','iHMPIBD',
#                                            'MHMGS','AUTCRCFENG','CHINACRCYU','FGCRC','ITCRC','USCRC',
#                                            'SWEDISH','CHINAQINJ','CHINALC','CHINAJIEZ','SWEDACVD',
#                                            'CHINARA','CHINAHYPER','CHINAYEZ')))



#dat.rescaled$caseControl<- factor(dat.rescaled$caseControl,
#                                  levels = c('0', '1'))


dim(dat.rescaled) # 2672 hosts in final dataset
table(dat.rescaled$caseControl) # 1184 healthy (caseControl == 1), 1488 non-healthy (caseControl == 0)


# Remove study MHIGC (li2014) entirely because its controls are all from spain and cases from denmark
# case/control confounded with country/cohort. Unclear if it was sequenced in same batch or not.
# drop it.

nrow(dat.rescaled)
nrow(prob.rescaled)
nrow(prob)

prob<- prob[which(dat.rescaled$internalName != 'MHIGC'),]
prob.rescaled<- prob.rescaled[which(dat.rescaled$internalName != 'MHIGC'),]
dat.rescaled<- dat.rescaled[which(dat.rescaled$internalName != 'MHIGC'),]


table(dat.rescaled$phenotype, dat.rescaled$internalName)



# (2.2) Cleanup and save ----

#save.image('./output/stats_runs/data_preps_25062021_FULL_ENVIRONMENT.RData')

sort( sapply(ls(),function(x){object.size(get(x))})) 

rm(prob.df, d, old, d.tmp, d.101, d.101.trim, old.kza, probspo.rescaled, probrspo.rescaled,
   probr.rescaled, cum.101, invSimpson, Shannon, Richness, logs2, logs, divs.index, cohorts, md, md.select,
   check, foo, logs2, logs, local_project_dir)

# keep prob.rescaled because using that matrix in plotting script
# Also keep prob for the new AUC analysis

sort( sapply(ls(),function(x){object.size(get(x))})) 

#save.image('./output/stats_runs/data_preps_23062021.RData')
#load('./output/stats_runs/data_preps_01062021.RData')
#load('./output/stats_runs/data_preps_23062021.RData')

