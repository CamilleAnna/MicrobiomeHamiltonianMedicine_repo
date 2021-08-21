FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness'
my_WORKING_DIR='/exports/eddie/scratch/s1687811'









qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness',my_WORKING_DIR='$TMPDIR' -t 222 job_caseControl_2021_fastqdump.sh


qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1 job_caseControl_2021_fastqdump_quick.sh
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness',my_WORKING_DIR='$TMPDIR' -t 2 job_caseControl_2021_fastqdump_quick.sh




qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 222 job_caseControl_2021_fastqdump_quick_TMPDIR.sh
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 2-5 job_caseControl_2021_fastqdump_quick.sh



# Previous script failed because of environment loading modeuls issue
# I had forgot a line at the top to initialise environment
# Fixed now
# also add a long line output to log file at begining, to distinguish different runs on same sample
# Re-launching, with time = 5h and vmem lowered to 5G to hopefully enter queue faster
# should work, surely will for sample 3
# probably for sample 2
# Not sure for sample 222, might run out of time and/or memory

qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 1 job_caseControl_2021_fastqdump_quick_TMPDIR.sh # should have empty RTF
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 2 job_caseControl_2021_fastqdump_quick_TMPDIR.sh # should work
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 3 job_caseControl_2021_fastqdump_quick_TMPDIR.sh # should work, 2 sp in snaps summary
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 222 job_caseControl_2021_fastqdump_quick_TMPDIR.sh # should work

# Ok this seems to be working fine. Just added some last edits:
# for the horizontal line to actually appear
# for error in midas species output to be correctly reported in log file
# and in this event, exiting script

qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 1 job_caseControl_2021_fastqdump_quick_TMPDIR.sh


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OK, this is working, final command launch is:
# launching tasks 4-293 (there are 293 hosts, start at 4, I know 4 and 5 shouldnot give  any output. Just for double checking tasks  1-3 now I know output is  fine)
# COMMAND LAUNCH FOR METAHIT COHORT SAMPLES

qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 4-293 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh







# METAHIT  MGS
cd  /exports/eddie/scratch/s1687811/caseControl_2021
# TEST onmyWORKINGDIR with low runtime and vmem, task1-3
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-MGS',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh
# TEST on TMPDIR with low runtime and vmem, task1-3
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-MGS' -t 1-3 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh


# KARLSSON
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/karlsson',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh
# TEST on TMPDIR
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/karlsson' -t 1-3 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh
#(9480173)


# QINT2D
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/QinT2D',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh
# TEST on TMPDIR
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/QinT2D' -t 1-3 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh
#(9480196)


# FENGCRC
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/FengCRC',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh
# TEST on TMPDIR
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/FengCRC' -t 1-3 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh
#(9480205)



# WENRAAS
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/WenRAAS',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh



# QINLC
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/QinLC',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh


# Zeller CRC
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/ZellerCRC',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 44 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh



# iHMP
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/iHMP',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh


# MetS
cd  /exports/eddie/scratch/s1687811/caseControl_2021
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetS',my_WORKING_DIR='/exports/eddie/scratch/s1687811' -t 1-3 job_caseControl_2021_fastqdump_myWORKINGDIR_19012021.sh


# Ok so
# let the ones running finish. Check the  MetS output
# let the current testing TMPDIR go but skip it for the others, there is no reason for itnot to work
# Then write the set of all commands
# Thenlaunch them





# FINAL COMMANDS

# METAHIT RICHNESS
# DONE qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-Richness' -t 1-291 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh


cd  /exports/eddie/scratch/s1687811/caseControl_2021

# METAHIT  MGS
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetaHIT-MGS' -t 4-138 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# KARLSSON
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/karlsson' -t 1-145 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# QINT2D
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/QinT2D' -t 1-370 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# FENGCRC
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/FengCRC' -t 1-156 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# WENRAAS
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/WenRAAS' -t 1-211 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# QINLC
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/QinLC' -t 1-123 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# Zeller CRC
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/ZellerCRC' -t 1-157 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# iHMP
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/iHMP' -t 1-106 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh

# MetS
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/caseControl_2021/MetS' -t 1-86 job_caseControl_2021_fastqdump_TMPDIR_19012021.sh




# NEW LAUNCH ON FEBRUARY 6



# run other tasks on a long version of scripts. Might have to re-run tasks 1 and 2
qsub -v FINAL_DIR='/exports/eddie/scratch/s1687811/feb6_caseControl_2021/processing',INFILE='assembled_CC_SRA.txt' -t 1-3478 job_caseControl_2021_fastqdump_TMPDIR_6Feb_LONG.sh











