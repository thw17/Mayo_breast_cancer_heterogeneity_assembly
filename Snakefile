configfile: "breast_cancer_config.json"

# hg19_ref_path =
# hg38_ref_path =

orig_bam_directory = "/mnt/storage/SAYRES/MAYO/"

# bwa_path =
# samtools_path =
# freebayes_path =
# gatk_path =
# picard_path =
xyalign_path = "/scratch/thwebste/xyalign_test/XYalign/xyalign/xyalign.py"

rule all:
	input:
		expand("xyalign/logfiles/{sample}_strip_reads_xyalign.log", sample=config["sample_list"])

rule strip_reads:
	input:
		bam = lambda wildcards: orig_bam_directory + config["sample_bams"][wildcards.sample]
	output:
		logfile = "xyalign/logfiles/{sample}_strip_reads_xyalign.log"
	params:
		xyalign = xyalign_path,
		sample_id = "{sample}_strip_reads"
	shell:
		"source activate xyalign_env && python {params.xyalign} --STRIP_READS --ref null --bam {input.bam} --sample_id {params.sample_id} --output_dir xyalign --chromosomes ALL"

# rule prepare_reference_hg19:
# 	input:
# 		hg19_ref_path
# 	output:
# 		fai = hg19_ref_path + ".fai",
# 		amb = hg19_ref_path + ".amb",
# 		dict = hg19_ref_prefix + ".dict"
# 	params:
# 		samtools = samtools_path,
# 		bwa = bwa_path
# 	run:
# 		# faidx
# 		shell("{params.samtools} faidx {input}")
# 		# .dict
# 		shell("{params.samtools} dict -o {output.dict} {input}")
# 		# bwa
# 		shell("{params.bwa} index {input}")
#
# rule prepare_reference_hg38:
# 	input:
# 		hg38_ref_path
# 	output:
# 		fai = hg38_ref_path + ".fai",
# 		amb = hg38_ref_path + ".amb",
# 		dict = hg38_ref_prefix + ".dict"
# 	params:
# 		samtools = samtools_path,
# 		bwa = bwa_path
# 	run:
# 		# faidx
# 		shell("{params.samtools} faidx {input}")
# 		# .dict
# 		shell("{params.samtools} dict -o {output.dict} {input}")
# 		# bwa
# 		shell("{params.bwa} index {input}")
#
# rule map_and_process_trimmed_reads_hg19:
# 	input:
# 		fq1 =
# 		fq2 =
# 		fai = hg19_ref_path + ".fai",
# 		ref = hg19_ref_path
# 	output:
# 		"processed_bams/{sample}.hg19.sorted.bam"
# 	params:
# 		id = lambda wildcards: config[wildcards.sample]["ID"],
# 		sm = lambda wildcards: config[wildcards.sample]["SM"],
# 		lb = lambda wildcards: config[wildcards.sample]["LB"],
# 		pu = lambda wildcards: config[wildcards.sample]["PU"],
# 		pl = lambda wildcards: config[wildcards.sample]["PL"],
# 		bwa = bwa_path,
# 		samtools = samtools_path
# 	threads: 4
# 	shell:
# 		" {params.bwa} mem -t {threads} -R "
# 		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
# 		"{input.ref} {input.fq1} {input.fq2}"
# 		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
# 		"-O bam -o {output}"
#
# rule map_and_process_trimmed_reads_hg38:
# 	input:
# 		fq1 =
# 		fq2 =
# 		fai = hg38_ref_path + ".fai",
# 		ref = hg38_ref_path
# 	output:
# 		"processed_bams/{sample}.hg38.sorted.bam"
# 	params:
# 		id = lambda wildcards: config[wildcards.sample]["ID"],
# 		sm = lambda wildcards: config[wildcards.sample]["SM"],
# 		lb = lambda wildcards: config[wildcards.sample]["LB"],
# 		pu = lambda wildcards: config[wildcards.sample]["PU"],
# 		pl = lambda wildcards: config[wildcards.sample]["PL"],
# 		bwa = bwa_path,
# 		samtools =  samtools_path
# 	threads: 4
# 	shell:
# 		" {params.bwa} mem -t {threads} -R "
# 		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
# 		"{input.ref} {input.fq1} {input.fq2}"
# 		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
# 		"-O bam -o {output}"
#
# rule freebayes_call_single_chrom:
# 	input:
# 		bam = expand("{sample}.sorted.rg.realigned.recal.FIXEDrg.bam", sample=samples),
# 		bai = expand("{sample}.sorted.rg.realigned.recal.FIXEDrg.bam.bai", sample=samples),
# 		ref = genome_path,
# 		target_chrom = lambda wildcards: config["chromosomes"][wildcards.chrom]
# 	output:
# 		"calls/all.{chrom}.raw.vcf"
# 	params:
# 		region = {chrom},
# 		freebayes = freebayes_path
# 	threads: 4
# 	shell:
# 		"{params.freebayes} -f {input.ref} --region {params.region} --pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities {input.bam} > {output}"
