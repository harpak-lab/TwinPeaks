#!/bin/bash

#SBATCH -J ChengRepGnomAD      # Job name
#SBATCH -o ChengRepGnomAD.o%j   # Name of stdout output file
#SBATCH -e ChengRepGnomAD.e%j   # Name of stderr error file
#SBATCH -p normal             # Queue (partition) name
#SBATCH -N 5                  # Total # of nodes 
#SBATCH -n 21                  # Total # of mpi tasks
#SBATCH -t 01:00:00           # Run time (hh:mm:ss)
#SBATCH --mail-type=all       # Send email at begin and end of job
#SBATCH --mail-user=mjm8356@utexas.edu

module load Rstats
module load launcher

echo "Tissue is $1"

mkdir -p ancs/FIN/$1
mkdir -p ancs/ASJ/$1
mkdir -p ancs/NFE/$1

> ChengRepGnomAD_FIN.$1.tasks
> ChengRepGnomAD_ASJ.$1.tasks
> =ChengRepGnomAD_NFE.$1.tasks

for i in {1..21};do

echo "Rscript /scratch/08312/mjm8356/ChengRep/gnomAD/ChengRepTwinPeaks_gnomad.R FIN $1 $i" >> ChengRepGnomAD_FIN.$1.tasks
echo "Rscript /scratch/08312/mjm8356/ChengRep/gnomAD/ChengRepTwinPeaks_gnomad.R ASJ $1 $i" >> ChengRepGnomAD_ASJ.$1.tasks
echo "Rscript /scratch/08312/mjm8356/ChengRep/gnomAD/ChengRepTwinPeaks_gnomad.R NFE $1 $i" >> ChengRepGnomAD_NFE.$1.tasks

done

export LAUNCHER_WORKDIR=/scratch/08312/mjm8356/ChengRep/gnomAD
export LAUNCHER_JOB_FILE=ChengRepGnomAD_FIN.$1.tasks

${LAUNCHER_DIR}/paramrun

export LAUNCHER_WORKDIR=/scratch/08312/mjm8356/ChengRep/gnomAD
export LAUNCHER_JOB_FILE=ChengRepGnomAD_ASJ.$1.tasks

${LAUNCHER_DIR}/paramrun

export LAUNCHER_WORKDIR=/scratch/08312/mjm8356/ChengRep/gnomAD
export LAUNCHER_JOB_FILE=ChengRepGnomAD_NFE.$1.tasks

${LAUNCHER_DIR}/paramrun

cat ./ancs/FIN/$1/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr1.txt > ./ancs/FIN/$1/TwinPeaksDF.gnomADv3_FIN.GTExv8.AllChr.txt
cat ./ancs/ASJ/$1/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr1.txt > ./ancs/ASJ/$1/TwinPeaksDF.gnomADv3_ASJ.GTExv8.AllChr.txt
cat ./ancs/NFE/$1/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr1.txt > ./ancs/NFE/$1/TwinPeaksDF.gnomADv3_NFE.GTExv8.AllChr.txt

for i in {2..21};do

sed '1d' ./ancs/FIN/$1/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr$i.txt >> ./ancs/FIN/$1/TwinPeaksDF.gnomADv3_FIN.GTExv8.AllChr.txt
sed '1d' ./ancs/ASJ/$1/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr$i.txt >> ./ancs/ASJ/$1/TwinPeaksDF.gnomADv3_ASJ.GTExv8.AllChr.txt
sed '1d' ./ancs/NFE/$1/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr$i.txt >> ./ancs/NFE/$1/TwinPeaksDF.gnomADv3_NFE.GTExv8.AllChr.txt

done

for i in FIN ASJ NFE;do

mv ancs/$i/Gonads/TwinPeaksDF.gnomADv3_$i.GTExv8.AllChr.txt ancs/$i/Gonads/TwinPeaksDF.gnomADv3_$i.GTExv8_Gonads.AllChr.txt
mv ancs/$i/Brain/TwinPeaksDF.gnomADv3_$i.GTExv8.AllChr.txt ancs/$i/Brain/TwinPeaksDF.gnomADv3_$i.GTExv8_Brain.AllChr.txt

done

# TO RUN
# sbatch ChengRepGnomAD.slurm Gonads
# sbatch ChengRepGnomAD.slurm Brain

# cat ./ancs/FIN/Gonads/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr1.txt > ./ancs/FIN/Gonads/TwinPeaksDF.gnomADv3_FIN_Gonads.GTExv8.AllChr.txt
# cat ./ancs/ASJ/Gonads/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr1.txt > ./ancs/ASJ/Gonads/TwinPeaksDF.gnomADv3_ASJ_Gonads.GTExv8.AllChr.txt
# cat ./ancs/NFE/Gonads/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr1.txt > ./ancs/NFE/Gonads/TwinPeaksDF.gnomADv3_NFE_Gonads.GTExv8.AllChr.txt

# cat ./ancs/FIN/Brain/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr1.txt > ./ancs/FIN/Brain/TwinPeaksDF.gnomADv3_FIN_Brain.GTExv8.AllChr.txt
# cat ./ancs/ASJ/Brain/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr1.txt > ./ancs/ASJ/Brain/TwinPeaksDF.gnomADv3_ASJ_Brain.GTExv8.AllChr.txt
# cat ./ancs/NFE/Brain/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr1.txt > ./ancs/NFE/Brain/TwinPeaksDF.gnomADv3_NFE_Brain.GTExv8.AllChr.txt

# for i in {2..21};do

# echo $i

# sed '1d' ./ancs/FIN/Gonads/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr$i.txt >> ./ancs/FIN/Gonads/TwinPeaksDF.gnomADv3_FIN_Gonads.GTExv8.AllChr.txt
# sed '1d' ./ancs/ASJ/Gonads/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr$i.txt >> ./ancs/ASJ/Gonads/TwinPeaksDF.gnomADv3_ASJ_Gonads.GTExv8.AllChr.txt
# sed '1d' ./ancs/NFE/Gonads/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr$i.txt >> ./ancs/NFE/Gonads/TwinPeaksDF.gnomADv3_NFE_Gonads.GTExv8.AllChr.txt

# sed '1d' ./ancs/FIN/Brain/TwinPeaksDF.gnomADv3_FIN.GTExv8.chr$i.txt >> ./ancs/FIN/Brain/TwinPeaksDF.gnomADv3_FIN_Brain.GTExv8.AllChr.txt
# sed '1d' ./ancs/ASJ/Brain/TwinPeaksDF.gnomADv3_ASJ.GTExv8.chr$i.txt >> ./ancs/ASJ/Brain/TwinPeaksDF.gnomADv3_ASJ_Brain.GTExv8.AllChr.txt
# sed '1d' ./ancs/NFE/Brain/TwinPeaksDF.gnomADv3_NFE.GTExv8.chr$i.txt >> ./ancs/NFE/Brain/TwinPeaksDF.gnomADv3_NFE_Brain.GTExv8.AllChr.txt

# done

# rsync -avrW ancs/FIN/Gonads/TwinPeaksDF.gnomADv3_FIN.GTExv8_Gonads.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/
# rsync -avrW ancs/ASJ/Gonads/TwinPeaksDF.gnomADv3_ASJ.GTExv8_Gonads.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/
# rsync -avrW ancs/NFE/Gonads/TwinPeaksDF.gnomADv3_NFE.GTExv8_Gonads.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/
# rsync -avrW ancs/FIN/Brain/TwinPeaksDF.gnomADv3_FIN.GTExv8_Brain.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/
# rsync -avrW ancs/ASJ/Brain/TwinPeaksDF.gnomADv3_ASJ.GTExv8_Brain.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/
# rsync -avrW ancs/NFE/Brain/TwinPeaksDF.gnomADv3_NFE.GTExv8_Brain.AllChr.txt mjm8356@kirkcomp01.ccbb.utexas.edu:/stor/work/Harpak/mjm8356/ChengRep/gnomAD/


