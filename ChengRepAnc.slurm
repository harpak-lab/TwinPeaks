#!/bin/bash

#SBATCH -J ChengRepAnc      # Job name
#SBATCH -o ChengRepAnc.o%j   # Name of stdout output file
#SBATCH -e ChengRepAnc.e%j   # Name of stderr error file
#SBATCH -p normal             # Queue (partition) name
#SBATCH -N 1                  # Total # of nodes 
#SBATCH -n 22                  # Total # of mpi tasks
#SBATCH -t 6:00:00           # Run time (hh:mm:ss)
#SBATCH --mail-type=all       # Send email at begin and end of job
#SBATCH --mail-user=mjm8356@utexas.edu

export PATH="$HOME/bcftools-1.16:$PATH"
module load launcher
module load Rstats

cd /scratch/08312/mjm8356/ChengRep/ChengRep/phasedVCF

# export anclist="GBR FIN CHS PUR CDX CLM IBS PEL PJL KHV ACB GWD ESN BEB MSL STU ITU CEU YRI CHB JPT LWK ASW MXL TSI GIH"
export anclist="$1 $2 $3 $4 $5 $6"
# export anclist="GBR FIN"
# for anc in $anclist;
# do 
# 	echo $anc;
# 	grep -w 'male' 1kgSampleSexesFilt.Anc.txt | grep -w "$anc" | awk -F '\t' '{print $1}' > ${anc}VCF/1kgSampleMales.${anc}.txt;
# 	grep -w 'female' 1kgSampleSexesFilt.Anc.txt | grep -w "$anc" | awk -F '\t' '{print $1}' > ${anc}VCF/1kgSampleFemales.${anc}.txt;
# done

for anc in $anclist;do

# grep -w 'male' 1kgSampleSexesFilt.Anc.txt | grep -w "$anc" | awk -F '\t' '{print $1}' > ${anc}VCF/1kgSampleMales.${anc}.txt
# grep -w 'female' 1kgSampleSexesFilt.Anc.txt | grep -w "$anc" | awk -F '\t' '{print $1}' > ${anc}VCF/1kgSampleFemales.${anc}.txt

# bcftools view -S GBRVCF/1kgSampleMales.GBR.txt -O z -o GBRVCF/chr1_males.GBR.vcf.gz ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
# bcftools view -S GBRVCF/1kgSampleFemales.GBR.txt -O z -o GBRVCF/chr1_females.GBR.vcf.gz ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz

# for i in {1..22};do echo "Getting Male VCF for Chr ${i} in anc ${anc}";bcftools view -S ${anc}VCF/1kgSampleMales.${anc}.txt -O z -o ${anc}VCF/chr${i}_males.${anc}.vcf.gz ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz;echo "Getting Female VCF for Chr ${i} in anc ${anc}";bcftools view -S ${anc}VCF/1kgSampleFemales.${anc}.txt -O z -o ${anc}VCF/chr${i}_females.${anc}.vcf.gz ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz;done

# > AncGetVCF.${anc}.tasks
# for i in {1..22};do echo "bcftools view -S ${anc}VCF/1kgSampleMales.${anc}.txt -O z -o ${anc}VCF/chr${i}_males.${anc}.vcf.gz ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz;bcftools view -S ${anc}VCF/1kgSampleFemales.${anc}.txt -O z -o ${anc}VCF/chr${i}_females.${anc}.vcf.gz ALL.chr${i}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz" >> AncGetVCF.${anc}.tasks;done

# export LAUNCHER_WORKDIR=/scratch/08312/mjm8356/ChengRep/ChengRep/phasedVCF
# export LAUNCHER_JOB_FILE=AncGetVCF.${anc}.tasks

# ${LAUNCHER_DIR}/paramrun

# for i in {1..22};do echo "Getting Male Info for Chr ${i}";bcftools view --types snps -Ou ${anc}VCF/chr${i}_males.${anc}.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' -o ${anc}VCF/chr${i}_males.${anc}.INFO;sed -i '1 i\CHROM\tPOS\tREF\tALT\tAC\tAN\tAF' ${anc}VCF/chr${i}_males.${anc}.INFO;echo "Getting Female Info for Chr ${i}";bcftools view --types snps -Ou ${anc}VCF/chr${i}_females.${anc}.vcf.gz | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' -o ${anc}VCF/chr${i}_females.${anc}.INFO;sed -i '1 i\CHROM\tPOS\tREF\tALT\tAC\tAN\tAF' ${anc}VCF/chr${i}_females.${anc}.INFO;done

cd /scratch/08312/mjm8356/ChengRep/ChengRep

> AncChrFST.${anc}.tasks
for i in {1..22};do echo "Rscript ChengRepTwinPeaks_Anc.R ${anc} ${i}" >> AncChrFST.${anc}.tasks;done

export LAUNCHER_WORKDIR=/scratch/08312/mjm8356/ChengRep/ChengRep
export LAUNCHER_JOB_FILE=AncChrFST.${anc}.tasks

${LAUNCHER_DIR}/paramrun

cd /scratch/08312/mjm8356/ChengRep/ChengRep/phasedVCF

> ${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt

for i in {1..22};do cat ${anc}VCF/TwinPeaksDF_chr${i}.${anc}.77v3.txt >> ${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt;done
sed -i '1 i\name\tdelta\tavgFST\tmaxFST\tavgFST_HUD\tmaxFST_HUD' ${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt

done

# export anclist="GBR FIN CHS PUR CDX CLM IBS PEL PJL KHV ACB GWD ESN BEB MSL STU ITU CEU YRI CHB JPT LWK ASW MXL TSI GIH"
# for i in $anclist;do echo $i;cp ./phasedVCF/${i}VCF/TwinPeaksDF_chrALL.${i}.txt ./AncTP/;done

# for anc in $anclist;do echo $anc;> ${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt;for i in {1..22};do cat ${anc}VCF/TwinPeaksDF_chr${i}.${anc}.77v3.txt >> ${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt;done;done

for anc in $anclist;do echo $anc;sed -i '1 i\name\tdelta\tavgFST\tmaxFST\tavgFST_HUD\tmaxFST_HUD' phasedVCF/${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt;cp phasedVCF/${anc}VCF/TwinPeaksDF_chrALL.${anc}.77v3.txt ./AncTP77v3/;done
