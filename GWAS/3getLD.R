for i in {1..7}
do
vcftools --vcf maf.imputed.all.vcf --chr scaffold_$i --positions Chr${i}GWASe-6.txt --recode --out Chr${i}GWASe-6.LD
done

for i in {1..7}
do
vcftools --vcf Chr{i}GWASe-6.LD.recode.vcf --geno-r2 --ld-window-bp 50000 --out Chr${i}-GWASe-6-50Kb
done