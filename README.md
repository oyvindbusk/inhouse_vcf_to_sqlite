# inhouse_vcf_to_sqlite
Collects inhouse vcf files and reads/writes to a sqlite


##### This is how the reported vcf should look:
```
##fileformat=VCFv4.1
##fileDate=20170821
##source=CentoMD
##reference=GRChr37
##INFO=<ID=IHSAMPLES,Number=1,Type=String,Description="Sample names from inhouse db that contains this variant and het state in parentheses">
##INFO=<ID=IHHET,Number=1,Type=String,Description="Number of inhouse hets for this variant">
##INFO=<ID=IHHOM,Number=1,Type=String,Description="Number of inhouse homozygotes for this variant">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
1       13281   .       C       G       .       .       IHSAMPLES=401_18(HET)_287_12(HET);IHHET=2;IHHOM=0
```