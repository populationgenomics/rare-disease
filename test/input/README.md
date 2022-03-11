

used to turn a subset of the multi-sample VCF into a single sample extract

```bash
bcftools view -H hail_classes.vcf.bgz.vcf | head | awk '{ print $1, $2, $3, $4, $5, $6, $7, $8, $9, "0/1:15,17:32:99:.:./.:.:603,0,432:.:.:." }'  > test_vcf_vars
```