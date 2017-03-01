



rule freebayes_call_single_chrom:
	input:
		bam=expand("{sample}.sorted.rg.realigned.recal.FIXEDrg.bam", sample=samples),
		bai=expand("{sample}.sorted.rg.realigned.recal.FIXEDrg.bam.bai", sample=samples),
		ref=genome_path
		target_chrom = lambda wildcards: config["chromosomes"][wildcards.chrom]
	output:
		"calls/all.{chrom}.raw.vcf"
	params:
		region={chrom}
	threads: 4
	shell:
		"freebayes -f {input.ref} --region {params.region} --pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities {input.bam} > {output}"
