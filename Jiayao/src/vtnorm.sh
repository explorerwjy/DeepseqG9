vt decompose /share/data/resources/ANNOVAR_DATA/gnomAD/WES/gnomad.exomes.r2.0.2.sites.vcf.gz -s -o gnomad_decomp.vcf
vt normalize gnomad_decomp.vcf -o gnomad_decomp_norm.vcf -r /share/data/resources/hg19/references/hg19.fasta
