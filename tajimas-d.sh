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

# Add AC-Info line to VCF file
bcftools +fill-tags \
SD-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf > 
SD-AC.vcf

bcftools +fill-tags \
EE-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf > 
EE-AC.vcf

bcftools +fill-tags \
SF-HC-gtALL.PE.bwa_mem.Sceloporus_undulatus_funannotate.passed.RG.sorted_passedSNPs.vcftoolsFiltered-2alleles-noIndels-hwe-geno-mind.vcf.recode.vcf > 
SF-AC.vcf
# Then run vcf-kit

vk tajima 100000 20000 \
SD-AC.vcf >
SD_output.100kb.20kb.TajimaD

vk tajima 100000 20000 \
EE-AC.vcf >
EE_output.100kb.20kb.TajimaD

vk tajima 100000 20000 \
SF-AC.vcf >
SF_output.100kb.20kb.TajimaD