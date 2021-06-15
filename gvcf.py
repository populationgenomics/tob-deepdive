# deep dive in the firsts teps of the pipeline: GVCFs


import hail as hl;

# All datasets in TOB-WGS are using GRCh38
hl.init(default_reference='GRCh38');

# read a gvcf to explore structure. Chose TOB1520 randomly.
gvcf = hl.import_vcf('gs://cpg-tob-wgs-test/gvcf/batch1/TOB1520.g.vcf.gz', min_partitions=12, force_bgz=True)

gvcf.describe()
gvcf.show()
gvcf.info.END.show()

hl.summarize_variants(gvcf)
# 4464040 SNPs
# as expected most variants are homrefs blocks. as expected the peseudo allele 'NONREF' is present

# the first step of the pipeline is to apply gatk ReblockGVCF. Let's check what that does
# open the reblocked gvcf for TOB1520. Got the URL from Vlad, seems to be a storage for intermediate results
rb_gvcf = hl.import_vcf('gs://cpg-tob-wgs-test-tmp/joint-calling/v2/hail/batch/728b87/1/output_gvcf.g.vcf.gz', force_bgz=True, min_partitions=12)
rb_gvcf.describe();
hl.summarize_variants(rb_gvcf);
# 4094043 SNPs
# we lost a lot of variants. Not only homref blocks, but also many SNPs

# investitate what variants were lost on SNPs
# extract SNPs
snps = hl.filter_alleles(gvcf, lambda a,i:hl.is_snp(gvcf.alleles[0], a))
snps.show()
hl.summarize_variants(snps)

# command below does not works because the alleles where changed by filter_alleles
lost_snps = snps.filter_rows(~ hl.is_defined(rb_gvcf.rows()[snps.row_key]))

# instead we have to use the old locus and alleles in the join expression
lost_snps = snps.filter_rows(~ hl.is_defined(rb_gvcf.rows()[snps.old_locus, snps.old_alleles]))

# some exploratory coed to figure out whether a variant is a SNP
allele_pairs = hl.range(1, hl.len(ht.alleles)).map(lambda i: (ht.alleles[0], ht.alleles[i]))
hl.any(allele_pairs.map(lambda pair: hl.is_snp(pair[0], pair[1]))).show()

# this is not a satisfactory way to filter out homref blocks
# a better way would be to filter out variants with only 2 alleles: ref and <NONREF>

no_block = gvcf.filter_rows(hl.len(gvcf.alleles)>2)
hl.summarize_variants(no_block)
lost_variants = no_block.filter_rows(~ hl.is_defined(rb_gvcf.rows()[no_block.row_key]))
hl.summarize_variants(lost_variants)
# we lost 565,345 variants with ReblockGVCF. Why ? not all of these have GQ==0
# ReblockGVCF: what criteria are used to remove variants ?

