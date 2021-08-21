library("DBI")
library("RSQLite")
library("SRAdb")
sqlfile <- '~/Documents/PhD/Research/HamiltonianMedicine/HamiltonianMedicine_work/data/caseControl_studies/SRAdb/SRAmetadb.sqlite'
#sqlfile<- 'data/caseControl_studies/SRAdb/SRAmetadb.sqlite'
sra_con <- dbConnect(SQLite(),sqlfile) # Create a connection for later queries



library(curatedMetagenomicData)


mh<- read.csv('~/Desktop/data_CC_04Feb/MetaHIT_all.txt', header=TRUE, sep = ',') 

mh2<- mh %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title)


mh.igc<- read.csv('~/Desktop/data_CC_04Feb/PRJEB5224_MHIGC_LiJ2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title)

mh.mgs<- read.csv('~/Desktop/data_CC_04Feb/PRJEB1220_MHMGS_NielsenB2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title)

mh.rich<- read.csv('~/Desktop/data_CC_04Feb/PRJEB4336_MHRICHNESS_LechatelierE2013_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = do.call('rbind', strsplit(Alias, '-'))[,2]) %>%
  select(-Alias) %>%
  rename(Alias = host)



mh.all<- rbind(mh.igc, mh.mgs, mh.rich)



# IGC, 109 unique hosts from denmark
unique(mh.igc[mh.igc$geographic_location_country_sea_region == 'Denmark',2])
# RICHNESS, 292 unique hosts, all from denmark
unique(do.call('rbind', strsplit(mh.rich[,2], '-'))[,2])
# MGS, 177 unique hosts from denmark
unique(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Denmark',2])


unique(mh.igc[mh.igc$geographic_location_country_sea_region == 'Denmark',2]) %in% unique(do.call('rbind', strsplit(mh.rich[,2], '-'))[,2])

unique(mh.igc[mh.igc$geographic_location_country_sea_region == 'Denmark',2]) %in% unique(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Denmark',2])

unique(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Denmark',2]) %in% unique(do.call('rbind', strsplit(mh.rich[,2], '-'))[,2])

mh.igc

length(unique(mh.all[mh.all$geographic_location_country_sea_region == 'Denmark', 2]))
292+109


# IGC, 118 unique hosts from Spain
igc<- unique(do.call('rbind', strsplit(mh.igc[mh.igc$geographic_location_country_sea_region == 'Spain',2], '-'))[,1])


# MGS, 219 unique hosts from Spain
mgs<- unique(do.call('rbind', strsplit(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Spain',2], '-'))[,1])


# IGC in MGS
table(unique(do.call('rbind', strsplit(mh.igc[mh.igc$geographic_location_country_sea_region == 'Spain',2], '-'))[,1]) %in%
        unique(do.call('rbind', strsplit(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Spain',2], '-'))[,1]))

# MGS in IGC
table(unique(do.call('rbind', strsplit(mh.mgs[mh.mgs$geographic_location_country_sea_region == 'Spain',2], '-'))[,1]) %in% 
        unique(do.call('rbind', strsplit(mh.igc[mh.igc$geographic_location_country_sea_region == 'Spain',2], '-'))[,1]) )


library(VennDiagram)
venn.diagram(x = list(igc, mgs), 
             category.names = c('igc', 'mgs'),
             filename = '~/Desktop/vd.spain.jpg')


d1<- data.frame(host = igc, dataset = 'IGC')
d2<- data.frame(host = mgs, dataset = 'MGS')
d3<- rbind(d1, d2)

nrow(d3)
unique(d3$host)

table(d3$dataset)
length(unique(igc))


hosts.rich<- mh.rich$Alias

hosts.mh.all<- c(hosts.rich, unique(d3$host))

length(unique(hosts.rich))
length(unique(mh.igc$host))



# METAHIT - RICHNESS ----

mh.rich<- read.csv('~/Desktop/data_CC_04Feb/PRJEB4336_MHRICHNESS_LechatelierE2013_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = do.call('rbind', strsplit(Alias, '-'))[,2])


MD_mh.rich<- combined_metadata[combined_metadata$dataset_name == 'LeChatelierE_2013',] %>%
  select(subjectID, country, gender, age, BMI, study_condition, disease) %>%
  mutate(phenotype = 'Healthy') %>%
  mutate(cohort = 'MHRICH') %>%
  rename(host = subjectID, bmi = BMI) %>%
  select(host, country, gender, age, bmi, phenotype, cohort)
  
  
mh.rich_md<- left_join(mh.rich, MD_mh.rich, by = 'host') # 595 rows, 292 hosts




# METAHIT MGS ----

mh.mgs<- read.csv('~/Desktop/data_CC_04Feb/PRJEB1220_MHMGS_NielsenB2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = do.call('rbind', strsplit(Alias, '-'))[,1]) %>%
  filter(!host %in% mh.rich$host) %>%
  arrange(Alias)

first_alias<- mh.mgs %>% # Keep first time point sample only, all runs available for them
  select(host, Alias) %>%
  mutate(first = !duplicated(host)) %>%
  filter(first == TRUE)
  
mh.mgs<- mh.mgs %>% filter(Alias %in% first_alias$Alias)

# Add metadata
MD_mh.mgs<- read_excel('~/Desktop/data_CC_04Feb/metadata_nielsen2014.xlsx',
                       col_types = c('text', 'text', 'text', 'numeric', 'text', 'numeric', 'numeric', 'text', 'text', 'numeric', 'numeric')) %>%
  select(host, country, gender, age, bmi, phenotype) %>%
  mutate(cohort = 'MHMGS') %>%
  arrange(host) %>%
  mutate(first = !duplicated(host)) %>%
  filter(first == TRUE) %>%
  select(-first)

mh.mgs_md<- left_join(mh.mgs, MD_mh.mgs, by = 'host') # 246 rows, 141 unique hosts, all from spain



# METAHIT IGC ----

mh.igc<- read.csv('~/Desktop/data_CC_04Feb/PRJEB5224_MHIGC_LiJ2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED',
         geographic_location_country_sea_region != 'China') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = do.call('rbind', strsplit(Alias, '-'))[,1]) %>%
  filter(!host %in% mh.mgs$host)

first_alias<- mh.igc %>% # Keep first time point sample only, all runs available for them
  select(host, Alias) %>%
  mutate(first = !duplicated(host)) %>%
  filter(first == TRUE)

mh.igc<- mh.igc %>% filter(Alias %in% first_alias$Alias)


# Add metadata
treat<- as.data.frame(read_excel('~/Desktop/data_CC_04Feb/metadata_forslund2015.xlsx'))
treat.bmi<- as.data.frame(read_excel('~/Desktop/data_CC_04Feb/metadata_forslund2015.xlsx', sheet = 2))[,1:2]
treat<- left_join(treat, treat.bmi, by = 'host')
sum(unique(mh.igc$host) %in% treat$host) # 106/109 are in forslund
t<- treat[treat$host %in% unique(mh.igc$host),] %>% select(-country)


MD_mh.igc<- combined_metadata[combined_metadata$dataset_name == 'LiJ_2014',] %>%
  select(subjectID, country, gender, age, BMI, study_condition, disease) %>%
  filter(country %in% c('DNK', 'ESP')) %>%
  mutate(host = do.call('rbind', strsplit(subjectID, '-'))[,1]) %>%
  left_join(t, by = 'host') %>%
  select(host, country, gender, age, bmi, disease, status) %>%
  mutate(cohort = 'MHIGC') %>%
  arrange(host) %>%
  mutate(first = !duplicated(host)) %>%
  filter(first == TRUE) %>%
  select(-first)


mh.igc_md<- left_join(mh.igc, MD_mh.igc, by = 'host') # 200 rows, 175 unique hosts
mh.igc_md$phenotype<- gsub(' metformin+', '', mh.igc_md$status, fixed = TRUE)
mh.igc_md$phenotype<- gsub(' metformin-', '', mh.igc_md$phenotype, fixed = TRUE)
mh.igc_md$phenotype<- ifelse(is.na(mh.igc_md$phenotype), 'Healthy', mh.igc_md$phenotype)
mh.igc_md$status<- mh.igc_md$phenotype
mh.igc_md<- mh.igc_md %>% select(-phenotype) %>% rename(phenotype = status)
#mh.igc_md<- mh.igc_md %>% rename(phenotype = status)
mh.igc_md<- mh.igc_md %>% select(-disease)


# COSTEA 2017 ----

costea.kazak<- read.csv('~/Desktop/data_CC_04Feb/PRJEB17632_KAZAK_Costea2017_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED',
         geographic_location_country_sea_region == 'Kazakhstan') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(time_point = do.call('rbind', strsplit(Alias, '-'))[,3]) %>%
  filter(time_point == 0) %>%
  select(-time_point) %>%
  mutate(host = Alias)


MD_costea.kazak<- read_excel('~/Desktop/data_CC_04Feb/metadata_costea2017.xlsx', sheet = 2,
                         col_types = c(rep('text', 4), 'numeric', 'numeric', 'text', 'text', rep('numeric', 3))) %>%
  filter(country == 'Kazakhstan') %>%
  select(host, country, gender, age, bmi, `Diabetes status`) %>%
  mutate(cohort = 'KAZAK') %>%
  mutate(time_point = do.call('rbind', strsplit(host, '-'))[,3]) %>%
  filter(time_point == 0) %>%
  select(-time_point) %>%
  mutate(`Diabetes status` = ifelse(`Diabetes status` == 'CTR', 'Healthy', 'MetS')) %>%
  rename(phenotype = `Diabetes status`)


costea.kazak_md<- left_join(costea.kazak, MD_costea.kazak, by = 'host') # 341 runs, 86 hosts


# Karlsson 2013 ----

karlsson.swedish<- read.csv('~/Desktop/data_CC_04Feb/PRJEB1786_SWEDISH_Karlsson2013_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, Description, Title) %>%
  mutate(host = paste0('S',Title))


treat<- as.data.frame(read_excel('~/Desktop/data_CC_04Feb/metadata_forslund2015.xlsx'))
treat.bmi<- as.data.frame(read_excel('~/Desktop/data_CC_04Feb/metadata_forslund2015.xlsx', sheet = 2))[,1:2]
treat<- left_join(treat, treat.bmi, by = 'host')
treat.swe<- treat[treat$country == 'SWE', c('host', 'bmi')] %>%
  mutate(host = gsub('NG-5425_', 'S', fixed = TRUE, host)) %>%
  mutate(host = gsub('NG-5636_', 'S', fixed = TRUE, host))


MD_karlsson.swedish<- combined_metadata[combined_metadata$dataset_name == 'KarlssonFH_2013',] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  mutate(study_condition = ifelse(study_condition == 'control', 'Healthy', study_condition)) %>%
  mutate(phenotype = study_condition) %>%
  rename(host = subjectID) %>%
  left_join(treat.swe, by = 'host') %>%
  mutate(BMI = bmi) %>%
  select(-bmi) %>%
  rename(bmi = BMI) %>%
  select(-study_condition) %>%
  mutate(cohort = 'SWEDISH')

karlsson.swedish_md<- left_join(karlsson.swedish, MD_karlsson.swedish, by = 'host') # 145 rows/runs, 145 hosts





# QinJ 2012 - T2D ---

QinJ.chinat2d<- read.csv('~/Desktop/data_CC_04Feb/PRJNA422434_CHINAQINJ_QinJ2012_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, SampleName, BioSample, BioProject) %>%
  mutate(host = gsub('bgi-', '', fixed = TRUE, SampleName))


MD_QinJ.chinat2d<- read_excel('~/Desktop/data_CC_04Feb/metadata_costea2017.xlsx', sheet = 'Sheet1') %>%
  filter(country == 'Chinese') %>%
  mutate(cohort = 'CHINAQINJ',
         bmi = as.numeric(bmi)) %>%
  select(host, country, gender, age, bmi, `Diabetes status`, cohort) %>%
  rename(phenotype = `Diabetes status`) %>%
  mutate(phenotype = ifelse(phenotype == 'CTR', 'Healthy', phenotype)) %>%
  filter(!is.na(phenotype))


QinJ.chinat2d_md<- left_join(QinJ.chinat2d, MD_QinJ.chinat2d, by = 'host') %>%
  filter(!is.na(phenotype)) # 368 rows.runs, 368 samples, (2/370 original drop because no metadata)


# JieZ 2017 - ACVD ----

jiez.china<- read.csv('~/Desktop/data_CC_04Feb/PRJEB21528_CHINAJIEZ_JieZ2017_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = Alias)


MD_jiez.china<- read_excel('~/Desktop/data_CC_04Feb/metadata_JieZ2017.xlsx', sheet = 2) %>%
  mutate(country = 'china') %>%
  select(host, country, gender, age, bmi, phenotype) %>%
  mutate(cohort = 'CHINAJIEZ') %>%
  mutate(phenotype = ifelse(phenotype == 0, 'Healthy', 'ACVD'))


jiez.china_md<- left_join(jiez.china, MD_jiez.china, by = 'host') # 405 rows/runs, 405 samples
nrow(jiez.china_md)
length(unique(jiez.china_md$host))



# Karlsson 2012


swed.acvd<- read.csv('~/Desktop/data_CC_04Feb/PRJNA177201_SWEDACVD_Karlsson2012_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, BioSample) %>%
  mutate(host = BioSample)

MD_swed.acvd<- read_excel('~/Desktop/data_CC_04Feb/GMHI_supp1.xlsx', sheet = 'reformat') %>%
  filter(BioSample %in% swed.acvd$BioSample) %>%
  mutate(country = 'Gothenburg', cohort = 'SWEDACVD') %>%
  rename(host = BioSample) %>%
  select(host, country, gender, age, bmi, phenotype, cohort) %>%
  mutate(phenotype = ifelse(phenotype == 'Overweight', 'Healthy', phenotype)) %>% # this was a re-classification from GMHI paper
  mutate(phenotype = ifelse(phenotype == 'Symptomatic atherosclerosis', 'ACVD', phenotype))


swed.acvd_md<- left_join(swed.acvd, MD_swed.acvd, by = 'host') %>%
  filter(!is.na(phenotype)) # 25 rows/runs, 25 hosts (2/27 had to phenotype data)



# YeZ - BD ----

YeZ.china<- read.csv('~/Desktop/data_CC_04Feb/PRJNA431482_CHINAYEZ_YeZ2018_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, BioSample, SampleName) %>%
  mutate(host = SampleName)

nrow(YeZ.china)
length(unique(YeZ.china$host))


MD_YeZ.china<- combined_metadata[combined_metadata$dataset_name == 'YeZ_2018',] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition) %>%
  mutate(cohort = 'CHINAYEZ',
         host = gsub('YEZ_', '', fixed = TRUE, subjectID)) %>%
  select(host, country, gender, age, bmi, phenotype, cohort) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype))

YeZ.china_md<- left_join(YeZ.china, MD_YeZ.china, by = 'host') %>%
  filter(!is.na(phenotype)) # 65 rows/runs, 65 hosts, (11/76 original had no phenotype data)


# CRC US - Vogtman ----

crc.vogtman<- read.csv('~/Desktop/data_CC_04Feb/PRJEB12449_USCRC_VogtmannE2016_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = Alias)


MD_crc.vogtman<- combined_metadata[combined_metadata$dataset_name == 'VogtmannE_2016',] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'ASCRC') %>%
  select(host, country, gender, age, bmi, phenotype, cohort)

crc.vogtman_md<- left_join(crc.vogtman, MD_crc.vogtman, by = 'host') %>% # 417 rows/runs, 104 hosts
  filter(!is.na(phenotype))



# THOMAS - ITALIAN CRC

crc.thomas<- read.csv('~/Desktop/data_CC_04Feb/SRP136711_ITCRC_ThomasAM2018_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, BioSample, SampleName, BioProject) %>%
  mutate(host = SampleName)


MD_crc.thomas<- combined_metadata[combined_metadata$dataset_name %in% c('ThomasAM_2018a', 'ThomasAM_2018b'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'ITCRC') 


crc.thomas_md<- left_join(crc.thomas, MD_crc.thomas, by = 'host')

nrow(crc.thomas_md)
length(unique(crc.thomas_md$host))
length(unique(crc.thomas_md$BioSample))
table(crc.thomas_md$phenotype)


# YuJ 2015 , CRC china ----

crc.yuJ<- read.csv('~/Desktop/data_CC_04Feb/PRJEB10878_CHINACRCYU_YuJ2015_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = Alias)

nrow(crc.yuJ)

MD_crc.yuJ<- combined_metadata[combined_metadata$dataset_name %in% c('YuJ_2015'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'CHINACRCYU') 


crc.yuJ_md<- left_join(crc.yuJ, MD_crc.yuJ, by = 'host')


# FENG2015 CRC ----

crc.FengQ<- read.csv('~/Desktop/data_CC_04Feb/PRJEB7774_AUTCRCFENG_FengQ2015_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = paste0('SID', Alias))

MD_crc.FengQ<- combined_metadata[combined_metadata$dataset_name %in% c('FengQ_2015'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'AUTCRCFENG') 

crc.FengQ_md<- left_join(crc.FengQ, MD_crc.FengQ, by = 'host') %>% # 154 runs/rows, 154 hosts, 2/156 had missing phenotype
  filter(!is.na(phenotype))


# CRC Zeller ----

crc.zeller<- read.csv('~/Desktop/data_CC_04Feb/PRJEB6070_FGCRC_ZellerG2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED', 
         AssayType == 'WGS') %>%
  select(Run, Alias, host_subject_id, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  rename(country = geographic_location_country_sea_region) %>%
  mutate(host = host_subject_id) %>%
  arrange(host)

nrow(crc.zeller) # 645 rows
length(unique(crc.zeller$host)) # 199 unique hosts
length(unique(crc.zeller$Alias)) # 199 unique Alias
length(unique(crc.zeller$BioSample)) # 199 unique BioSample


MD_crc.zeller<- combined_metadata[combined_metadata$dataset_name %in% c('ZellerG_2014'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'FGCRC') 


crc.zeller_md<- left_join(crc.zeller, MD_crc.zeller, by = 'host') # 645 runs/rows, 199 hosts


# LiJ 2017 ----

hyper.liJ<- read.csv('~/Desktop/data_CC_04Feb/PRJEB13870_CHINAHYPER_LiJ2017_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = Alias)


MD_hyper.liJ<- combined_metadata[combined_metadata$dataset_name %in% c('LiJ_2017'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'CHINAHYPER') 


hyper.liJ_md<- left_join(hyper.liJ, MD_hyper.liJ, by = 'host') # 196 runs/rows, 196 hosts


# QinN_2014 ----

lc.qinN<- read.csv('~/Desktop/data_CC_04Feb/PRJEB6337_CHINALC_QinN2014_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = do.call('rbind', strsplit(Alias, '_'))[,1])


MD_lc.qinN<- combined_metadata[combined_metadata$dataset_name %in% c('QinN_2014'),] %>%
  select(subjectID, country, gender, age, BMI, study_condition) %>%
  rename(bmi = BMI, phenotype = study_condition, host = subjectID) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'CHINALC') 


MD_lc.qinN$host<- gsub('-', '', MD_lc.qinN$host, fixed = TRUE)
lc.qinN$host<- gsub('-', '', lc.qinN$host, fixed = TRUE)


lc.qinN_md<- left_join(lc.qinN, MD_lc.qinN, by = 'host') # 314 runs/rows, 237 hosts



# Zhang_2015 ----

md.zhang2015<- read_excel('~/Desktop/data_CC_04Feb/metadata_schmidt.xlsx') %>%
  filter(cohort == 'Zhang-RA' & material == 'stool' & timepoint == 0) %>%
  rename(host = BioSample_accession,
         gender = sex,
         phenotype = group) %>%
  mutate(country = 'china') %>%
  select(host, country, gender, age, bmi, phenotype) %>%
  mutate(phenotype = ifelse(phenotype == 'control', 'Healthy', phenotype),
         cohort = 'CHINARA')

Zhang.2015<- read.csv('~/Desktop/data_CC_04Feb/PRJEB6997_CHINARA_Zhang2015_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED',
         Title != 'oral sample from China') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  filter(BioSample %in% md.zhang2015$host) %>%
  mutate(host = BioSample)



Zhang.2015_md<- left_join(Zhang.2015, md.zhang2015, by = 'host') # 192 runs/rows, 192 hosts

nrow(Zhang.2015_md) # 192 rows
length(unique(Zhang.2015_md$host))



# Liu 2017 ----

ob.liu<- read.csv('~/Desktop/data_CC_04Feb/PRJEB12123_CHINAOB_Liu2017_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED') %>%
  select(Run, Alias, BioSample, SampleName, Sample_Name, BioProject, SRA_accession, environment_material, geographic_location_country_sea_region, project_name, Title) %>%
  mutate(host = BioSample) %>%
  mutate(check = do.call('rbind', strsplit(Alias, '-'))[,1])
length(unique(ob.liu$Alias))
length(unique(ob.liu$check))
length(unique(ob.liu$BioSample))


MD_ob.liu<- read_excel('~/Desktop/data_CC_04Feb/GMHI_supp1.xlsx', sheet = 'reformat') %>%
  filter(BioSample %in% ob.liu$BioSample) %>%
  mutate(country = 'China', cohort = 'CHINAOB') %>%
  rename(host = BioSample) %>%
  select(host, country, gender, age, bmi, phenotype, cohort)

ob.liu_md<- left_join(ob.liu, MD_ob.liu, by = 'host') %>%
  filter(!is.na(phenotype)) # 205 rows/runs, 205 hosts (52/257 had to phenotype data)



# Llyoy price ibd ----

ibd.lloyd_md<- read.csv('~/Desktop/data_CC_04Feb/PRJNA398089_iHMPIBD_Lloyd2019_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED',
         AssayType == 'WGS',
         Organism == 'human gut metagenome') %>%
  select(Run, BioSample, SampleName, host_subject_id, BioProject, visit_num, geo_loc_name, Host_Age, host_sex, host_disease) %>%
  mutate(host = host_subject_id)  %>%
  arrange(host_subject_id, SampleName) %>%
  mutate(first = !duplicated(host_subject_id)) %>%
  filter(first == TRUE) %>%
  select(-first) %>%
  mutate(phenotype = host_disease) %>%
  mutate(phenotype = gsub("Crohn''s disease", 'CD', fixed = TRUE, phenotype)) %>%
  mutate(phenotype = gsub("Ulcerative colitis", 'UC', fixed = TRUE, phenotype)) %>%
  mutate(phenotype = gsub("ulcerative colitis", 'UC', fixed = TRUE, phenotype)) %>%
  mutate(phenotype = ifelse(phenotype == '', 'Healthy', phenotype)) %>%
  mutate(cohort = 'iHMPIBD', bmi = NA) %>%
  rename(age = Host_Age,
         gender = host_sex,
         country = geo_loc_name) %>%
  select(Run, BioSample, SampleName, host_subject_id, BioProject, host, country, gender, age, bmi, phenotype, cohort)

table(ibd.lloyd_md$phenotype)
length(unique(ibd.lloyd_md$BioSample))
length(unique(ibd.lloyd_md$host))
nrow(ibd.lloyd_md)

# Complete with data from schirmer 2018

ibd.schirmer<- read.csv('~/Desktop/data_CC_04Feb/PRJNA389280_iHMPSchirmer_Schirmer2018_e.txt', header=TRUE, sep = '\t') %>%
  filter(LibraryLayout == 'PAIRED',
         data_type == 'MGX',
         Organism == 'human gut metagenome') %>%
  select(Run, BioSample, SampleName, host_subject_id, BioProject, geo_loc_name, host) %>%
  arrange(host_subject_id, SampleName)

gmhi<- read_excel('~/Desktop/data_CC_04Feb/GMHI_supp1.xlsx', sheet = 'reformat')

t<- left_join(ibd.schirmer, gmhi, by = 'BioSample') %>%
  filter(!is.na(status)) %>%
  filter(!host %in% ibd.lloyd_md$host)
unique(t$phenotype)
t$phenotype<- gsub("Crohn's disease", 'CD', fixed = TRUE, t$phenotype)
t$phenotype<- gsub("Ulcerative colitis", 'UC', fixed = TRUE, t$phenotype)
t$country = 'USA'
t$bmi = NA
t$cohort = 'iHMPIBD'

colnames(ibd.lloyd_md)

t<- t %>%
  select(Run, BioSample, SampleName, host_subject_id, BioProject, host, country, gender, age, bmi, phenotype, cohort)

ibd.lloyd_md<- rbind(ibd.lloyd_md, t)


nrow(ibd.lloyd_md)
length(unique(ibd.lloyd_md$host))
length(unique(ibd.lloyd_md$BioSample))

ibd.lloyd_md$host<- paste0('HMP', ibd.lloyd_md$host)



# He 2017 ----

ibd.HeJ<- read.csv('~/Desktop/data_CC_04Feb/PRJEB15371_CHINAIBD_He2017_e.txt', header=TRUE, sep = ',') %>%
  filter(LibraryLayout == 'PAIRED')


gmhi<- read_excel('~/Desktop/data_CC_04Feb/GMHI_supp1.xlsx', sheet = 'reformat')

t<- left_join(ibd.HeJ, gmhi, by = 'BioSample') %>%
  filter(!is.na(status))
t$phenotype[t$phenotype %in% c('Overweight', 'Underweight', 'Obesity')]<- 'Healthy'
t$phenotype<- gsub("Crohn's disease", 'CD', fixed = TRUE, t$phenotype)
t$bmi<- as.numeric(t$bmi)
t$age<- as.numeric(t$age)
t$country = 'China'
t$cohort = 'CHINAIBD'

t$host<- t$BioSample


ibd.HeJ_md<- t %>%
  select(Run, BioSample, SampleName, BioProject, host, country, gender, age, bmi, phenotype, cohort) %>% # 99 runs.rows, 99 samples, 24/123 samples had missing phentoype, but 14 of them also were EEN treatment samples
  arrange(host)

save.image("~/Desktop/data_CC_04Feb/assembly_session.RData")



# ASSEMBLE ALL OF THIS ----
swed.acvd_md$BioProject = 'PRJNA177201'
YeZ.china_md$BioProject = 'PRJNA431482'

mh.rich_md.trim<- mh.rich_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
mh.mgs_md.trim<- mh.mgs_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
mh.igc_md.trim<- mh.igc_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
costea.kazak_md.trim<- costea.kazak_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
karlsson.swedish_md.trim<- karlsson.swedish_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
QinJ.chinat2d_md.trim<- QinJ.chinat2d_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
jiez.china_md.trim<- jiez.china_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
swed.acvd_md.trim<- swed.acvd_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
YeZ.china_md.trim<- YeZ.china_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
crc.vogtman_md.trim<- crc.vogtman_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
crc.thomas_md.trim<- crc.thomas_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
crc.yuJ_md.trim<- crc.yuJ_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
crc.FengQ_md.trim<- crc.FengQ_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
crc.zeller_md.trim<- crc.zeller_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country.x, cohort) %>% rename(country = country.x)
hyper.liJ_md.trim<- hyper.liJ_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
lc.qinN_md.trim<- lc.qinN_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
Zhang.2015_md.trim<- Zhang.2015_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
ob.liu_md.trim<- ob.liu_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
ibd.lloyd_md.trim<- ibd.lloyd_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
ibd.HeJ_md.trim<- ibd.HeJ_md %>% select(Run, host, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)


data.all<- rbind(mh.rich_md.trim,
                  mh.mgs_md.trim,
                  mh.igc_md.trim,
                  costea.kazak_md.trim,
                  karlsson.swedish_md.trim,
                  QinJ.chinat2d_md.trim,
                  jiez.china_md.trim,
                  swed.acvd_md.trim,
                  YeZ.china_md.trim,
                  crc.vogtman_md.trim,
                  crc.thomas_md.trim,
                  crc.yuJ_md.trim,
                  crc.FengQ_md.trim,
                  crc.zeller_md.trim,
                  hyper.liJ_md.trim,
                  lc.qinN_md.trim,
                  Zhang.2015_md.trim,
                  ob.liu_md.trim,
                  ibd.lloyd_md.trim,
                  ibd.HeJ_md.trim)


data.all<- rename(data.all, run = Run)
data.all$cohort[data.all$cohort == 'ASCRC']<- 'USCRC'

nrow(data.all)
length(unique(data.all$run)) # check, all runs are unique, 5002
length(unique(data.all$host)) # check, all runs are unique, 3478
sum(is.na(data.all$phenotype)) # check, all rows have a phenotype


data.all %>%
  group_by(cohort) %>%
  summarise(nb_runs = length(unique(run)))

data.all %>%
  group_by(cohort) %>%
  summarise(nb_hosts = length(unique(host)))



data.all$host2<- data.all$host
data.all$host2<- gsub('.', '', data.all$host2, fixed = TRUE)
data.all$host2<- gsub('-', '', data.all$host2, fixed = TRUE)
data.all$host2<- gsub('_', '', data.all$host2, fixed = TRUE)

data.all$host2<- paste0(data.all$cohort, '_', data.all$host2)

#data.all %>%
#  select(phenotype, host, host2, cohort) %>%
#  unique() %>%
#  View()


data.all<- data.all %>%
  rename(host_original = host) %>%
  rename(host = host2) %>%
  select(run, host, host_original, BioSample, BioProject, gender, age, bmi, phenotype, country, cohort)
  

unique(data.all$phenotype)


write.table(data.all, '~/Desktop/data_CC_04Feb/assembled_CC_SRA.txt', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


