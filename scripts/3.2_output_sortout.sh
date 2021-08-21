
# OVERVIEW of what ran/did not ran
cd ./logs
tail -n 1 * >> ../tails.txt # use sublime and excel to process



# GET HOST NAMES WITH SNPS OUTPUT
# Only if pipeline goes through all steps the endpoint  is copied.
# Gets coppied if "MIDAS snps ran fine", "No species satisfied filtering criteria for snps estimates"
# 3424/3478 jobs ran
# 3075 of them had snps
# 90 of them had species but no snps
# That's 3165 midas output that got copied to local directory
# I want to isolate the 3075 for which snps ran
ls ./*/snps/summary.txt >> ../hosts_with_SNPS_output.txt


# GET THEIR NUMBER OF READS
cat hosts_with_SNPS_output_edited.txt | while read line;
do
	grep 'nReads original' ./logs/$line\.log >> nb_reads_originals.txt
	grep 'nReads processed' ./logs/$line\.log >> nb_reads_processed.txt
done



# COPY THOSE HOSTS SPECIES ABUNDANCE TABLES TO LOCAL DIRECTORY
cat hosts_with_SNPS_output_edited.txt | while read line;
do
cat ./midas_output/$line\.midas_output/species/species_profile.txt | awk -v var="$line" '{ FS = OFS ="\t" }{print $0,var}' | grep -v 'species_id' >> assembled_species_profiles.txt ;
done
gzip assembled_species_profiles.txt
scp s1687811@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s1687811/feb6_caseControl_2021/processing/assembled_species_profiles.txt.gz .
gunzip assembled_species_profiles.txt.gz