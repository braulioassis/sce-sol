# Separate VCF into populations
vcftools --vcf \
HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf
--keep SDsamples.txt \
--recode --out
SD-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf

vcftools --vcf \
HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf
--keep SFsamples.txt \
--recode --out
SF-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf

vcftools --vcf \
HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf
--keep EEsamples.txt \
--recode --out
EE-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf

# Separate vcf into a single vcf per contig

cat 
EE-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep "##" > EEheader.txt
 
cat 
EE-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" > EEno_header.vcf

for chr in {1..24}; do awk -v chr=scaffold_$chr '$1 == chr' EEno_header.vcf 
> EE.chrom$chr.vcf; done


cat 
SD-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep "##" > SDheader.txt
 
cat 
SD-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" > SDno_header.vcf

for chr in {1..24}; do awk -v chr=scaffold_$chr '$1 == chr' SDno_header.vcf 
> SD.chrom$chr.vcf; done


cat 
SF-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep "##" > SFheader.txt
 
cat 
SF-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" > SFno_header.vcf

for chr in {1..24}; do awk -v chr=scaffold_$chr '$1 == chr' SFno_header.vcf 
> SF.chrom$chr.vcf; done

# Add header back to each vcf

cat 
EE-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" | awk 'NR==1{print}' > EEcolumn_header.txt

cat EEheader.txt EEcolumn_header.txt > EEmain_header.vcf

cat 
SD-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" | awk 'NR==1{print}' > SDcolumn_header.txt

cat SDheader.txt SDcolumn_header.txt > SDmain_header.vcf

cat 
SF-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf 
| grep -v "##" | awk 'NR==1{print}' > SFcolumn_header.txt

cat SFheader.txt SFcolumn_header.txt > SFmain_header.vcf

# automate addition of header to each chrom. then rename for lassi.
for i in {1..24}; do cat EEmain_header.txt EE.chrom$i.vcf > 
EE.chrom$i.vcf.head; done

for i in {1..24}; do mv EE.chrom$i.vcf.head EE.chrom$i.vcf; done

for i in {1..24}; do cat SFmain_header.txt SF.chrom$i.vcf > 
SF.chrom$i.vcf.head; done

for i in {1..24}; do mv SF.chrom$i.vcf.head SF.chrom$i.vcf; done

for i in {1..24}; do cat SDmain_header.txt SD.chrom$i.vcf > 
SD.chrom$i.vcf.head; done

for i in {1..24}; do mv SD.chrom$i.vcf.head SD.chrom$i.vcf; done

#################################################################################

for i in {1..24}
do
  lassip --vcf 
EE.chrom$i.vcf 
--hapstats --winsize 201 --k 20 --calc-spec --winstep 100 --out 
EE.chrom$i 
--pop 
ids.pop.txt 
done

for i in {1..24}
do
  lassip --vcf 
SF.chrom$i.vcf 
--hapstats --winsize 201 --k 20 --calc-spec --winstep 100 --out 
SF.chrom$i 
--pop 
ids.pop.txt 
done

for i in {1..24}
do
  lassip --vcf 
SD.chrom$i.vcf 
--hapstats --winsize 201 --k 20 --calc-spec --winstep 100 --out 
SD.chrom$i 
--pop 
ids.pop.txt 
done

###################################################################################

# Calc average spectra to acquire null background
lassip --spectra 
*.EE.lassip.hap.spectra.gz 
--avg-spec --out 
EE.chromAll

lassip --spectra 
*.SF.lassip.hap.spectra.gz 
--avg-spec --out 
SF.chromAll

lassip --spectra 
*.SD.lassip.hap.spectra.gz 
--avg-spec --out 
SD.chromAll

###################################################################################
lassip --spectra 
EE.chrom${CHR}.EE.lassip.hap.spectra.gz 
--salti --out 
EE.chrom${CHR} 
--null-spec 
EE.chromAll.lassip.null.spectra.gz 

lassip --spectra 
SF.chrom${CHR}.SF.lassip.hap.spectra.gz 
--salti --out 
SF.chrom${CHR} 
--null-spec 
SF.chromAll.lassip.null.spectra.gz 

lassip --spectra 
SD.chrom${CHR}.SD.lassip.hap.spectra.gz 
--salti --out 
SD.chrom${CHR} 
--null-spec 
SD.chromAll.lassip.null.spectra.gz 

###################################################################################

# automate removal of column head from lassioutput files to concatenate
for i in {1..22}; do awk 'NR>1' EE.chrom$i.lassip.hap.out > 
EE.chrom$i.lassip.hap.out.nohead; done

for i in {1..22}; do awk 'NR>1' SF.chrom$i.lassip.hap.out > 
SF.chrom$i.lassip.hap.out.nohead; done

for i in {1..22}; do awk 'NR>1' SD.chrom$i.lassip.hap.out > 
SD.chrom$i.lassip.hap.out.nohead; done

# now concatenate
ls -v EE.chrom** | xargs cat > EE_ALL.saltilassi.out
ls -v SF.chrom** | xargs cat > SF_ALL.saltilassi.out
ls -v SD.chrom** | xargs cat > SD_ALL.saltilassi.out
