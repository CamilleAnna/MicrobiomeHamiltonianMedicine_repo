#!/bin/sh 
#$ -o /exports/eddie/scratch/s1687811/joblogs
#$ -e /exports/eddie/scratch/s1687811/joblogs
#$ -N CC60221
#$ -l h_rt=48:00:00
#$ -l h_vmem=8G
#$ -P bio_MPGS_csimonet



# Initialise the environment modules
. /etc/profile.d/modules.sh


cd $TMPDIR
#cd $my_WORKING_DIR # use this if want to run job on a local directory to see output as processing happens. Specify variable my_WORKING_DIR in the qsub command, along with the $FINAL_DIR variable

# Modules import
export PATH=$PATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/sratoolkit.2.10.9-ubuntu64/bin


# Get host name corresponding to task and corresponding accession numbers
HOST=$(cat $FINAL_DIR/$INFILE | sed '1d'  | cut -f 2 | sort | uniq | awk "NR==$SGE_TASK_ID")
ACCESSIONS=$(grep -w $HOST $FINAL_DIR/$INFILE | cut -f 1)

echo '___________________________________________________________________' >> $FINAL_DIR/logs/$HOST.log

# Download fastq files in a host-specific directory, gzip files
mkdir $HOST
cd ./$HOST
echo "$HOST (task $SGE_TASK_ID): starting download" >> $FINAL_DIR/logs/$HOST.log
fastq-dump --split-3 --gzip $ACCESSIONS

file1=$(ls *_1.fastq.gz)
file2=$(ls *_2.fastq.gz)

if [ -n "$file1" ] && [ -n "$file2" ]; then echo "$HOST (task $SGE_TASK_ID): Download complete" >> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): Download failed. Exiting"  >>  $FINAL_DIR/logs/$HOST.log; exit; fi


# Rename files for compatibility with MOCAT
rename fastq fq *
rename _ . *

# Assemble files from different lanes in one, move into a host sub-directory for MOCAT, discard lanes files
cat *.1.fq.gz > $HOST.1.fq.gz
cat *.2.fq.gz > $HOST.2.fq.gz

mkdir $HOST
mv  $HOST.1.fq.gz ./$HOST/$HOST.1.fq.gz
mv  $HOST.2.fq.gz ./$HOST/$HOST.2.fq.gz
rm *.fq.gz

# Process reads with MOCAT (RTF). Generate command, then runs as sh script
cp $FINAL_DIR/MOCAT.cfg .
echo $HOST > sample_list

export PERL5LIB=$PYTHONPATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MOCAT/src
/exports/csce/eddie/biology/groups/mcnally/camille/programs/MOCAT//src/MOCAT.pl -sf sample_list -rtf -x
if [ -f MOCATJob_readtrimfilter.sample_list_*.1.sh ] ; then echo "$HOST (task $SGE_TASK_ID): MOCAT command generation ran" >> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): issue generating MOCAT command. Exiting">> $FINAL_DIR/logs/$HOST.log; exit; fi

sh *.1.sh


size_rtf=$(du -sh ./$HOST/reads.processed.solexaqa/*.pair.1.fq.gz | cut -f 1)
if [ $size_rtf != 0 ]; then echo "$HOST (task $SGE_TASK_ID): read trim filter OK" >> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): read trim filter left no read. Exiting" >> $FINAL_DIR/logs/$HOST.log; exit; fi


# HOUSEKEEPING
size_original=$(du -sh  ./$HOST/*.1.fq.gz | cut -f 1)
nbReads_original=$(cat ./$HOST/stats/$HOST.1.fq.gz.raw.reads.stats | cut -f 1 | sed 1d)
nbReads_processed_total=$(cat ./$HOST/stats/$HOST.readtrimfilter.solexaqa.stats | cut -f 1 | sed 1d)
nbReads_processed=$(expr $nbReads_processed_total / 2)

echo "$HOST (task $SGE_TASK_ID): original size: $size_original" >> $FINAL_DIR/logs/$HOST.log
echo "$HOST (task $SGE_TASK_ID): rtf size: $size_rtf" >> $FINAL_DIR/logs/$HOST.log
echo "$HOST (task $SGE_TASK_ID): nReads original: $nbReads_original" >> $FINAL_DIR/logs/$HOST.log
echo "$HOST (task $SGE_TASK_ID): nReads processed: $nbReads_processed" >> $FINAL_DIR/logs/$HOST.log

rm MOCAT*
rm -r logs
rm sample_list
rm ./$HOST/*.fq.gz

mv ./$HOST/reads.processed.solexaqa/* .
rm *.single.fq.gz
rm *qual_stats
rm -r $HOST


# 4) RUN MIDAS

if [ -f *pair.1.fq.gz ] && [ -f *pair.2.fq.gz ] ; then echo "$HOST (task $SGE_TASK_ID): running MIDAS" >> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): missing at least one processed read files. Exiting">> $FINAL_DIR/logs/$HOST.log; exit; fi


# Initialise the environment modules
# . /etc/profile.d/modules.sh # Already done at top of script

module load anaconda
source activate mypythonMIDAS


# Update environment
export PYTHONPATH=$PYTHONPATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS
export PATH=$PATH:/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/scripts
export MIDAS_DB=/exports/csce/eddie/biology/groups/mcnally/camille/programs/MIDAS/midas_db_v1.2


#cd /exports/eddie/scratch/s1687811/FINAL/midas_output
mkdir $HOST.midas_output

run_midas.py species ./$HOST.midas_output -1 *pair.1.fq.gz -2 *pair.2.fq.gz --remove_temp #-n 10000
if [ -f ./$HOST.midas_output/species/species_profile.txt ] ; then echo "$HOST (task $SGE_TASK_ID): MIDAS species ran fine" >> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): no MIDAS species output. Possible reason: rtf output near empty. Exiting" >> $FINAL_DIR/logs/$HOST.log;exit; fi


run_midas.py snps ./$HOST.midas_output -1 *pair.1.fq.gz -2 *pair.2.fq.gz --remove_temp #-n 100000
if [ -f ./$HOST.midas_output/snps/summary.txt ]; then echo "$HOST (task $SGE_TASK_ID): MIDAS snps ran fine">> $FINAL_DIR/logs/$HOST.log; else echo "$HOST (task $SGE_TASK_ID): No species satisfied filtering criteria for snps estimates">> $FINAL_DIR/logs/$HOST.log; fi


# 5) MOVE FINAL OUTPUT TO STORE DIRECTORY

mv ./$HOST.midas_output/ $FINAL_DIR/midas_output/












