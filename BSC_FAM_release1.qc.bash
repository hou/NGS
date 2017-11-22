bcftools reheader -s convert_ids_release1.txt -o BSC_FAM.release1.snps.indels.vqsr.final.vcf.gz BSC_FAM.release1.snps.indels.vqsr.vcf.gz
bcftools reheader -s convert_ids_release1.txt -o BSC_FAM.release1_intersect.regions.snps.indels.vqsr.final.vcf.gz BSC_FAM.release1_intersect.regions.snps.indels.vqsr.vcf.gz

bcftools stats  -f PASS,. -I -v -d 0,500,1 -s - BSC_FAM.release1.snps.indels.vqsr.final.vcf.gz >BSC_FAM.release1.snps.indels.vqsr.final.vcf.gz.stats
bcftools stats  -f PASS,. -I -v -d 0,500,1 -s - BSC_FAM.release1_intersect.regions.snps.indels.vqsr.final.vcf.gz >BSC_FAM.release1_intersect.regions.snps.indels.vqsr.final.vcf.gz.stats

#vi BSC_FAM.release1_intersect.regions.snps.indels.vqsr.vcf.gz.stats
#vi BSC_FAM.release1.snps.indels.vqsr.vcf.gz.stats
plot-vcfstats -p BSC_FAM_Intersect_VQSR BSC_FAM.release1_intersect.regions.snps.indels.vqsr.final.vcf.gz.stats
plot-vcfstats -p BSC_FAM_Union_VQSR BSC_FAM.release1.snps.indels.vqsr.final.vcf.gz.stats

python3 /data/houl3/Softwares/NGS/vcfQC.py --GQ 20 --DP 10 --AB 0.25 --MISSING 1.0 BSC_FAM.release1.snps.indels.vqsr.final.vcf.gz BSC_FAM.release1.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf
python3 /data/houl3/Softwares/NGS/vcfQC.py --GQ 20 --DP 10 --AB 0.25 --MISSING 0.25 BSC_FAM.release1_intersect.regions.snps.indels.vqsr.final.vcf.gz BSC_FAM.release1_intersect.regions.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf
bgzip BSC_FAM.release1_intersect.regions.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf -f
tabix -p vcf BSC_FAM.release1_intersect.regions.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf.gz -f
bgzip -f BSC_FAM.release1.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf
tabix -p vcf -f BSC_FAM.release1.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf.gz
bcftools view -k -v 'snps' -f "PASS,." -O z -o BSC_FAM.release1_intersect.regions.snps.vqsr.GQ20.DP10.AB0.25.known.pass.vcf.gz BSC_FAM.release1_intersect.regions.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf.gz
plink1.9 --vcf BSC_FAM.release1_intersect.regions.snps.vqsr.GQ20.DP10.AB0.25.known.pass.vcf.gz --vcf-filter --keep-allele-order --make-bed --out BSC_FAM.vqsr.known.passOnly --double-id
awk '$1>22' BSC_FAM.vqsr.known.passOnly.bim >BSC_FAM.vqsr.known.passOnly.nonAutosomes.snps

#GWAS_QC.sh BSC_FAM.vqsr.known.passOnly hg19
cp BSC_FAM.vqsr.known.passOnly.fam BSC_FAM.vqsr.known.passOnly.fam.archive
R --vanilla <BSC_FAM.release1.fam.file.R
plink --bfile BSC_FAM.vqsr.known.passOnly --bmerge /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --make-bed --out BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged --geno 0.05 --noweb
plink --bfile BSC_FAM.vqsr.known.passOnly --exclude BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged.missnp --make-bed --out BSC_FAM.vqsr.known.passOnly.noBadSNPs --noweb
plink --bfile BSC_FAM.vqsr.known.passOnly.noBadSNPs --bmerge /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bed /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.bim /data/houl3/Share/GWAS_QC/hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps.fam --make-bed --out BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged --geno 0.05 --noweb
plink --bfile BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged --keep /data/houl3/Share/GWAS_QC/CEU.subjects --exclude /data/houl3/Share/GWAS_QC/high-LD-regions-hg19.txt --range --indep-pairwise 50 5 0.2 --out BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged --noweb
plink --bfile BSC_FAM.vqsr.known.passOnly --extract BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged.prune.in --exclude BSC_FAM.vqsr.known.passOnly.nonAutosomes.snps --maf 0.05 --geno 0.05 --genome --out BSC_FAM.vqsr.known.passOnly --noweb
R --vanilla <BSC_FAM_release1_relcheck_plot.R

python3 /data/houl3/Softwares/NGS/vcfPedcheck.py -z BSC_FAM.release1.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf.gz BSC_FAM_release1.modified.ped BSC_FAM.release1.snps.indels.vqsr.QCed.pedCheck
python3 /data/houl3/Softwares/NGS/vcfPedcheck.py -z BSC_FAM.release1_intersect.regions.snps.indels.vqsr.GQ20.DP10.AB0.25.vcf.gz BSC_FAM_release1.modified.ped BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.pedCheck
python3 /data/houl3/Softwares/NGS/vcfQC.py --GQ 20 --DP 10 --AB 0.25 --MISSING 0.25 BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.pedCheck.vcf BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf
python3 /data/houl3/Softwares/NGS/vcfQC.py --GQ 20 --DP 10 --AB 0.25 --MISSING 1.0 BSC_FAM.release1.snps.indels.vqsr.QCed.pedCheck.vcf BSC_FAM.release1.snps.indels.vqsr.QCed.final.vcf
bgzip -f BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf
tabix -p vcf -f BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf.gz
bgzip -f BSC_FAM.release1.snps.indels.vqsr.QCed.final.vcf
tabix -p vcf -f BSC_FAM.release1.snps.indels.vqsr.QCed.final.vcf.gz
bcftools stats  -f PASS,. -I -v -d 0,500,1 -s - BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf.gz >BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf.gz.stats
plot-vcfstats -p BSC_FAM_Intersect_Final -t BSC_FAM BSC_FAM.release1_intersect.regions.snps.indels.vqsr.QCed.final.vcf.gz.stats
bcftools stats  -f PASS,. -I -v -d 0,500,1 -s - BSC_FAM.release1.snps.indels.vqsr.QCed.final.vcf.gz >BSC_FAM.release1.snps.indels.vqsr.QCed.final.vcf.gz.stats
plink --bfile BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged --extract BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged.prune.in --exclude BSC_FAM.vqsr.known.passOnly.nonAutosomes.snps --maf 0.05 --geno 0.05 --make-bed --out BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged.pruned --noweb
king -b BSC_FAM.vqsr.known.passOnly.hapmap3r2.merged.pruned.bed --kinship --prefix BSC_FAM.hapmap3r2.merged
