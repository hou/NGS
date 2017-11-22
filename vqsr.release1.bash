gatk_bundle_dir=/data/houl3/BSC_FAM/gatk_bundle
ref=human_g1k_v37_decoy.fasta
dbSNP=dbsnp_138.b37.vcf
indel_1000=1000G_phase1.indels.b37.vcf
indel_Mills=Mills_and_1000G_gold_standard.indels.b37.vcf
gatk_dir=/data/houl3/Softwares/gatk_3.5
picard_dir=/data/houl3/Softwares/picard-tools-2.0.1
#target_regions=BSC_FAM.release1_target_regions.interval_list
omni_1000=1000G_omni2.5.b37.vcf
snps_1000=1000G_phase1.snps.high_confidence.b37.vcf
hapmap=hapmap_3.3.b37.vcf

java -Xmx4g -jar $gatk_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $gatk_bundle_dir/$ref -input BSC_FAM.release1.raw.snps.indels.vcf.gz -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatk_bundle_dir/$hapmap  -resource:omni,known=false,training=true,truth=true,prior=12.0 $gatk_bundle_dir/$omni_1000 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $gatk_bundle_dir/$snps_1000 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_bundle_dir/$dbSNP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile BSC_FAM.release1.raw.snps.vcf.vqsr.recal -tranchesFile BSC_FAM.release1.raw.snps.vcf.vqsr.tranches -rscriptFile BSC_FAM.release1.raw.snps.vcf.vqsr.R.plot #-an DP -an InbreedingCoeff

java -Xmx4g -jar $gatk_dir/GenomeAnalysisTK.jar -T VariantRecalibrator -R $gatk_bundle_dir/$ref -input BSC_FAM.release1.raw.snps.indels.vcf.gz --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 $gatk_bundle_dir/$indel_Mills -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatk_bundle_dir/$dbSNP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 90.0 -recalFile BSC_FAM.release1.raw.indels.vcf.vqsr.recal -tranchesFile BSC_FAM.release1.raw.indels.vcf.vqsr.tranches -rscriptFile BSC_FAM.release1.raw.indels.vcf.vqsr.R.plot #-an DP -an InbreedingCoeff

java -Xmx4g -jar $gatk_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $gatk_bundle_dir/$ref -input BSC_FAM.release1.raw.snps.indels.vcf.gz -recalFile BSC_FAM.release1.raw.snps.vcf.vqsr.recal -tranchesFile BSC_FAM.release1.raw.snps.vcf.vqsr.tranches --ts_filter_level 99.5 -mode SNP -o BSC_FAM.release1.snps.vqsr.vcf

java -Xmx4g -jar $gatk_dir/GenomeAnalysisTK.jar -T ApplyRecalibration -R $gatk_bundle_dir/$ref -input BSC_FAM.release1.snps.vqsr.vcf -recalFile BSC_FAM.release1.raw.indels.vcf.vqsr.recal -tranchesFile BSC_FAM.release1.raw.indels.vcf.vqsr.tranches --ts_filter_level 99.0 -mode INDEL -o BSC_FAM.release1.snps.indels.vqsr.vcf

bgzip BSC_FAM.release1.snps.indels.vqsr.vcf -f
tabix -p vcf -f BSC_FAM.release1.snps.indels.vqsr.vcf.gz 
bcftools stats  -f PASS,. -I -v -d 0,500,1 -s - BSC_FAM.release1.snps.indels.vqsr.vcf.gz >BSC_FAM.release1.snps.indels.vqsr.vcf.gz.stats
#vi BSC_FAM.release1.snps.indels.vqsr.vcf.gz.stats
#plot-vcfstats -t Ambigen -p Ambigen_union_VQSR BSC_FAM.release1.snps.indels.vqsr.vcf.gz.stats
