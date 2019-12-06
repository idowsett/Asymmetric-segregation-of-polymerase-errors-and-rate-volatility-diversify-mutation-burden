#!/bin/bash
set -e

#make subdirectories to keep things clean.  We usually run this script from within the directory that contains the fasta files for the strains listed after "strain_name".
mkdir tmp
mkdir discordant
mkdir split

REFPATH=~/Bioinformatics/genomes
SLICENUM=1
#export $REFPATH
strain_name='AH164_2_2 AH164_2_3 AH164_10_1 AH164_10_2 AH164_10_4 AH164_11_1 AH164_11_2 AH164_11_3 AH164_12_1 AH164_12_3 AH164_12_4'

#Alignment of strain genomes against entire yeast genome

count=0
for strain in $strain_name; do
bwa mem -M $REFPATH/yeast/S288C-masked-genome.fasta ${strain}_seq1_fixed.fq.gz ${strain}_seq2_fixed.fq.gz |samblaster -M -d discordant/${strain}.disc.sam -s split/${strain}.split.sam | grep -v SA:Z |grep -v XA: | samtools view -Sb -q 10 - |samtools sort -o tmp/${strain}_pe.unique_sorted.bam - &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/picard.jar AddOrReplaceReadGroups INPUT=tmp/${strain}_pe.unique_sorted.bam OUTPUT=tmp/${strain}_pe.fixedhdr.bam RGID=1 RGLB=1 RGPL=illumina   RGPU=TTAATA RGSM=${strain} VALIDATION_STRINGENCY=LENIENT &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    samtools index tmp/${strain}_pe.fixedhdr.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T RealignerTargetCreator -I tmp/${strain}_pe.fixedhdr.bam -R ${REFPATH}/yeast/S288C-masked-genome.fasta -o ${strain}.intervals &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T IndelRealigner -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe.fixedhdr.bam -targetIntervals ${strain}.intervals -o tmp/${strain}_pe-realigned.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T LeftAlignIndels -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe-realigned.bam -o tmp/${strain}_pe-leftaligned.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T BaseRecalibrator -I tmp/${strain}_pe-leftaligned.bam -R ${REFPATH}/yeast/S288C-masked-genome.fasta -knownSites ${REFPATH}/yeast/BY4741-diploid_snp_sorted_final.vcf -o tmp/${strain}_pe-recalibrated.grp &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

for strain in $strain_name; do
    java -Xmx1g -jar ~/Bioinformatics/programs/GATK3/GenomeAnalysisTK.jar -T PrintReads -R ${REFPATH}/yeast/S288C-masked-genome.fasta -I tmp/${strain}_pe-leftaligned.bam -BQSR tmp/${strain}_pe-recalibrated.grp -o ${strain}_pe-recalibrated.bam &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
	samtools mpileup -Bf $REFPATH/yeast/S288C-masked-genome.fasta ${strain}_pe-recalibrated.bam > ${strain}_pe.mpileup &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
    java -jar ~/Bioinformatics/programs/VarScan.v2.3.9.jar mpileup2snp ${strain}_pe.mpileup --min-coverage 18 --min-var-freq .1 --strand-filter 1 |python ~/Bioinformatics/programs/variant_deSNPer2013.py -s ~/Bioinformatics/genomes/yeast/AH0401.snp > ${strain}.noSNP.var &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
java -jar ~/Bioinformatics/programs/VarScan.v2.3.9.jar mpileup2indel ${strain}_pe.mpileup --min-coverage 18 --min-var-freq .22 --strand-filter 1 |python ~/Bioinformatics/programs/variant_deSNPer2013.py -s ~/Bioinformatics/genomes/yeast/AH0401indel.snp > ${strain}.no-indelSNP.var &
        let count+=1
        [[ $((count%${SLICENUM})) -eq 0 ]] && wait
done

count=0
for strain in $strain_name; do
samtools view -Sb discordant/${strain}.disc.sam|samtools sort -o ${strain}.disc.sort.bam -
samtools view -Sb split/${strain}.split.sam|samtools sort -o ${strain}.split.sort.bam -
samtools index ${strain}.disc.sort.bam
samtools index ${strain}.split.sort.bam
rm *.intervals
rm -r tmp/
rm -r discordant
rm -r split
done
