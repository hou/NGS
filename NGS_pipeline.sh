#!/usr/bin/env bash

# Liping Hou 
# liping.hou@nih.gov
# Feb, 2016

# bwa, samtools, picard, and gatk need to be installed before you run this script.
if [ $# -ne 1 ] && [ $# -ne 2 ]; then
    echo
    echo "Usage: NGS_pipeline.sh <BAM_file/FASTQ_Files>"
    echo "  For example:" 
    echo "    NGS_pipeline.sh id1.bam"
    echo "  OR:"
    echo "    NGS_pipeline.sh id1_R1.fastq.gz id1_R2.fastq.gz"
    echo
    echo "This script implements the GATK Best Practices (gatk v3.5). This is what this script will do:"
    echo
    echo "If the input file is a BAM file, it will be converted back to FASTQ files then followed by:"
    echo " 1) mapping (bwa mem)" 
    echo " 2) mark duplicates & sort (Picard)" 
    echo " 3) Base recalibration" 
    echo " 4) Variant calling (HaplotypeCaller GVCF mode)"
    echo "If the input files are two FASTQ files, the mapping will be started right away."
    echo
    exit
fi

#modify the following lines as needed
data_dir=/data/houl3/Amish/DRIFT-NIMH/Freeze2/BAM
gatk_bundle_dir=/data/houl3/BSC_FAM/gatk_bundle
ref=human_g1k_v37_decoy.fasta
dbSNP=dbsnp_138.b37.vcf
indel_1000=1000G_phase1.indels.b37.vcf
indel_Mills=Mills_and_1000G_gold_standard.indels.b37.vcf
gatk_dir=/data/houl3/Softwares/gatk_3.5
picard_dir=/data/houl3/Softwares/picard-tools-2.0.1
#target_regions=xgen-exome-research-panel-targets_Padded.v37.interval_list
target_regions=/data/houl3/Amish/DRIFT-NIMH/Freeze1/BAM/xgen-exome-research-panel-targets_100bp_Padded.v37.interval_list

#https://wiki.dnanexus.com/Scientific-Notes/human-genome
#http://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it
#https://www.broadinstitute.org/gatk/guide/article?id=1247
# Download GATK's resource bundle
if [ ! -f "$gatk_bundle_dir/$ref" ]; then
    echo "Warning: $gatk_bundle_dir/$ref not found!"
    echo "---> wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37_decoy.*"
    wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37_decoy.*
    gunzip $gatk_bundle_dir/human_g1k_v37_decoy.*.gz
fi

if [ ! -f "$gatk_bundle_dir/human_g1k_v37_decoy.fasta.bwt" ]; then
    echo "Warning: can not find the index files for bwa!"
    echo "---> bwa index -a bwtsw human_g1k_v37_decoy.fasta"
    bwa index -a bwtsw human_g1k_v37_decoy.fasta
fi

if [ ! -f "$gatk_bundle_dir/$dbSNP" ]; then
    echo "Warning: $gatk_bundle_dir/$dbSNP not found!"
    echo "---> wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.*"
    wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/dbsnp_138.b37.*
    gunzip $gatk_bundle_dir/dbsnp_138.b37.*.gz
fi

if [ ! -f "$gatk_bundle_dir/$indel_1000" ]; then
    echo "Warning: $gatk_bundle_dir/$indel_1000 not found!"
    echo "---> wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.indels.b37.vcf.*"
    wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.indels.b37.vcf.*
    gunzip $gatk_bundle_dir/1000G_phase1.indels.b37.vcf.*.gz
fi

if [ ! -f "$gatk_bundle_dir/$indel_Mills" ]; then
    echo "Warning: $gatk_bundle_dir/$indel_Mills not found!"
    echo "---> wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.*"
    wget -P $gatk_bundle_dir ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.*
    gunzip $gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.b37.vcf.*.gz
fi

#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf.*
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/hapmap_3.3.b37.vcf.*
#wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/b37/1000G_omni2.5.b37.vcf.*
#http://gatkforums.broadinstitute.org/wdl/discussion/4133/when-should-i-use-l-to-pass-in-a-list-of-intervals

echo "Analysis started: $(date)"
if [ $# -eq 1 ]; then
    # Step 0: convert BAM to fastq files
    id=$(basename $1 .realign.bam)
    echo "Converting BAM file to FASTQ files"
    echo "---> samtools fastq -1 ${id}_R1.fastq -2 ${id}_R2.fastq -n -O $1"
    samtools fastq -1 ${id}_R1.fastq -2 ${id}_R2.fastq -n -O $1
    echo "---> python3 fixFastq.py ${id}_R1.fastq ${id}_R2.fastq ${id}_fixed"
    python3 $data_dir/fixFastq.py ${id}_R1.fastq ${id}_R2.fastq ${id}_fixed
fi

if [ $# -eq 2 ]; then
    echo "Checking $1 and $2 to make sure the encoding is correct (recode them to illumina-1.8 otherwise)"
    n=$(basename $1 .fastq.gz)
    id=${n%_*}
    echo "---> python3 fixFastq.py $1 $2 ${id}_fixed --checkEncodingOnly"
    python3 fixFastq.py $1 $2 ${id}_fixed --checkEncodingOnly
fi

# Step 1: Mapping

if [ -e "${id}_fixed_R1.fastq" ]; then
    echo "Aligning ${id}_fixed_R1.fastq and ${id}_fixed_R2.fastq in paired analysis mode"
    echo "---> bwa mem -t 16 -M -R \"@RG\tID:BSC\tSM:$id\tPL:Illumina\" $gatk_bundle_dir/$ref ${id}_fixed_R1.fastq ${id}_fixed_R2.fastq | gzip >${id}.sam.gz"
    bwa mem -t 16 -M -R "@RG\tID:BSC\tSM:$id\tPL:Illumina" $gatk_bundle_dir/$ref ${id}_fixed_R1.fastq ${id}_fixed_R2.fastq | gzip >${id}.sam.gz

    if [ $? -eq 0 ]; then
        rm ${id}_*.fastq
    fi

else
    echo "Aligning $1 and $2 in paired analysis mode"
    echo "---> bwa mem -t 16 -M -R \"@RG\tID:BSC\tSM:$id\tPL:Illumina\" $gatk_bundle_dir/$ref $1 $2 | gzip >${id}.sam.gz"
    bwa mem -t 16 -M -R "@RG\tID:BSC\tSM:$id\tPL:Illumina" $gatk_bundle_dir/$ref $1 $2 | gzip >${id}.sam.gz
fi

echo "Converting ${id}.sam.gz to ${id}_aln.bam"
echo "---> samtools view -bT $gatk_bundle_dir/$ref -o ${id}_aln.bam ${id}.sam.gz"
samtools view -bT $gatk_bundle_dir/$ref -o ${id}_aln.bam ${id}.sam.gz

if [ $? -eq 0 ]; then
    rm ${id}.sam.gz
fi

echo "Sorting ${id}_aln.bam"

echo "---> java -Xmx16g -jar $picard_dir/picard.jar SortSam I=${id}_aln.bam O=${id}_sorted.bam SO=coordinate CREATE_INDEX=true"
java -Xmx16g -jar $picard_dir/picard.jar SortSam I=${id}_aln.bam O=${id}_sorted.bam SO=coordinate CREATE_INDEX=true

if [ $? -eq 0 ]; then
    rm ${id}_aln.bam
fi

## Step 2: Duplicate marking
echo "Marking duplicate reads"

echo "---> java -Xmx8g -jar $picard_dir/picard.jar MarkDuplicates I=${id}_sorted.bam O=${id}_mkdup.bam M=${id}_markdup_metrics.txt VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true"
java -Xmx8g -jar $picard_dir/picard.jar MarkDuplicates I=${id}_sorted.bam O=${id}_mkdup.bam M=${id}_markdup_metrics.txt VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true

if [ $? -eq 0 ]; then
    rm ${id}_sorted.bam
fi

# Step 3: Base quality recalibration
echo "Base quality recalibration"

echo "---> java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -I ${id}_mkdup.bam -R $gatk_bundle_dir/$ref -knownSites $gatk_bundle_dir/$dbSNP -knownSites $gatk_bundle_dir/$indel_Mills -knownSites $gatk_bundle_dir/$indel_1000 -o ${id}_mkdup_recal.table -nct 8 -L $target_regions"
java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -I ${id}_mkdup.bam -R $gatk_bundle_dir/$ref -knownSites $gatk_bundle_dir/$dbSNP -knownSites $gatk_bundle_dir/$indel_Mills -knownSites $gatk_bundle_dir/$indel_1000 -o ${id}_mkdup_recal.table -nct 8 -L $target_regions

echo "---> java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -I ${id}_mkdup.bam -R $gatk_bundle_dir/$ref -knownSites $gatk_bundle_dir/$dbSNP -knownSites $gatk_bundle_dir/$indel_Mills -knownSites $gatk_bundle_dir/$indel_1000 -BQSR ${id}_mkdup_recal.table -o ${id}_mkdup_after_recal.table -nct 8 -L $target_regions"
java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator -I ${id}_mkdup.bam -R $gatk_bundle_dir/$ref -knownSites $gatk_bundle_dir/$dbSNP -knownSites $gatk_bundle_dir/$indel_Mills -knownSites $gatk_bundle_dir/$indel_1000 -BQSR ${id}_mkdup_recal.table -o ${id}_mkdup_after_recal.table -nct 8 -L $target_regions

echo "---> java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $gatk_bundle_dir/$ref -before ${id}_mkdup_recal.table -after ${id}_mkdup_after_recal.table -plots ${id}_recal_plots.pdf"
java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T AnalyzeCovariates -R $gatk_bundle_dir/$ref -before ${id}_mkdup_recal.table -after ${id}_mkdup_after_recal.table -plots ${id}_recal_plots.pdf

echo "---> java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads -R $gatk_bundle_dir/$ref -I ${id}_mkdup.bam -BQSR ${id}_mkdup_recal.table -o ${id}_mkdup_recal.bam -nct 8"
java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads -R $gatk_bundle_dir/$ref -I ${id}_mkdup.bam -BQSR ${id}_mkdup_recal.table -o ${id}_mkdup_recal.bam -nct 8

if [ $? -eq 0 ]; then
    rm ${id}_mkdup.bam
fi

# Step 4: Variant calling
echo "Variant calling by HaplotypeCaller in GVCF mode"

echo "---> java -Xmx16g -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -R $gatk_bundle_dir/$ref -I ${id}_mkdup_recal.bam --emitRefConfidence GVCF --dbsnp $gatk_bundle_dir/$dbSNP -o ${id}.raw.snps.indels.g.vcf -nct 8 -L $target_regions"
java -Xmx16g -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -R $gatk_bundle_dir/$ref -I ${id}_mkdup_recal.bam --emitRefConfidence GVCF --dbsnp $gatk_bundle_dir/$dbSNP -o ${id}.raw.snps.indels.g.vcf -nct 8 -L $target_regions

echo "Analysis finished: $(date)"
