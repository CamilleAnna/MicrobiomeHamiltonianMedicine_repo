# On each table:

# replace geographic_location_(country_and/or_sea,region) by geographic_location_country_sea_region
# replace space by ''
# replace ( by ''
# replace ) by ''


cat PRJEB1220_MHMGS_NielsenB2014.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB1220_MHMGS_NielsenB2014_e.txt
cat PRJEB4336_MHRICHNESS_LechatelierE2013.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB4336_MHRICHNESS_LechatelierE2013_e.txt
cat PRJEB5224_MHIGC_LiJ2014.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB5224_MHIGC_LiJ2014_e.txt


cat PRJEB17632_KAZAK_Costea2017.txt | sed '1 s/geographic_location_(country_and\/or_sea)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB17632_KAZAK_Costea2017_e.txt


cat PRJEB1786_SWEDISH_Karlsson2013.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB1786_SWEDISH_Karlsson2013_e.txt
cat PRJNA422434_CHINAQINJ_QinJ2012.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA422434_CHINAQINJ_QinJ2012_e.txt


cat PRJEB21528_CHINAJIEZ_JieZ2017.txt | sed '1 s/geographic_location_(country_and\/or_sea)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB21528_CHINAJIEZ_JieZ2017_e.txt


cat PRJNA177201_SWEDACVD_Karlsson2012.txt  | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA177201_SWEDACVD_Karlsson2012_e.txt 
cat PRJNA431482_CHINAYEZ_YeZ2018.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA431482_CHINAYEZ_YeZ2018_e.txt

cat PRJEB12449_USCRC_VogtmannE2016.txt  | sed '1 s/geographic_location_(country_and\/or_sea)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB12449_USCRC_VogtmannE2016_e.txt 

cat PRJEB10878_CHINACRCYU_YuJ2015.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB10878_CHINACRCYU_YuJ2015_e.txt

cat PRJEB6997_CHINARA_Zhang2015.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB6997_CHINARA_Zhang2015_e.txt



cat PRJEB12123_CHINAOB_Liu2017.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB12123_CHINAOB_Liu2017_e.txt
cat PRJEB7774_CHINACRCFENG_FengQ2015.txt | sed '1 s/geographic_location_(country_and\/or_sea)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB7774_CHINACRCFENG_FengQ2015_e.txt
cat PRJEB13870_CHINAHYPER_LiJ2017.txt | sed '1 s/geographic_location_(country_and\/or_sea)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB13870_CHINAHYPER_LiJ2017_e.txt
cat PRJEB15371_CHINAIBD_He2017.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB15371_CHINAIBD_He2017_e.txt
cat PRJEB6070_FGCRC_WellerG2014.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB6070_FGCRC_WellerG2014_e.txt
cat PRJEB6337_CHINALC_QinN2014.txt | sed '1 s/geographic_location_(country_and\/or_sea,region)/geographic_location_country_sea_region/g' | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJEB6337_CHINALC_QinN2014_e.txt
cat PRJNA231909_EUROT1D_KosticAD2015.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA231909_EUROT1D_KosticAD2015_e.txt
cat PRJNA389280_iHMPSchirmer_Schirmer2018.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA389280_iHMPSchirmer_Schirmer2018_e.txt
cat PRJNA389927_USCRCHannigan_HanniganGD2017.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA389927_USCRCHannigan_HanniganGD2017_e.txt
cat SRP136711_ITCRC_ThomasAM2018.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > SRP136711_ITCRC_ThomasAM2018_e.txt
cat PRJNA398089_iHMPlloyd_Lloyd2019.txt | sed '1 s/ //g' | sed '1 s/(//g' | sed '1 s/)//g' > PRJNA398089_iHMPlloyd_Lloyd2019_e.txt
