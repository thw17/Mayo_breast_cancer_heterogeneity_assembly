configfile: "breast_cancer_config.json"

# Removed "PS13-585-B3-P3-rep" from the sample list due to inconsistency
# file name and read group name inside.


hg19_ref_path = "reference/ucsc.hg19.fasta"
hg19_ref_prefix = "reference/ucsc.hg19"
hg38_ref_path = "reference/Homo_sapiens_assembly38.fasta"
hg38_ref_prefix = "reference/Homo_sapiens_assembly38"

dbsnp_138_hg19_path = "misc/dbsnp_138.hg19.vcf"
dbsnp_138_hg38_path = "misc/dbsnp_138.hg38.vcf"
dpsnp_146_hg38_path = "misc/dbsnp_146.hg38.vcf"
mills_1kg_indels_hg19 = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
mills_1kg_indels_hg38 = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf"

orig_bam_directory = "/mnt/storage/SAYRES/MAYO/"
temp_dir_path = "temp/"

bwa_path = "bwa"
samtools_path = "samtools"
samblaster_path = "samblaster"
bgzip_path = "bgzip"
tabix_path = "tabix"
# freebayes_path =
gatk_path = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"
# picard_path =
xyalign_path = "/scratch/thwebste/xyalign_test/XYalign/xyalign/xyalign.py"

rule all:
	input:
		expand("xyalign/logfiles/{sample}_strip_reads_xyalign.log", sample=config["sample_list"]),
		expand("xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz", sample=config["sample_list"]),
		expand("xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz", sample=config["sample_list"]),
		expand("processed_bams/{sample}.hg38.sorted.bam", sample=config["sample_list"]),
		expand("processed_bams/{sample}.hg19.sorted.bam", sample=config["sample_list"]),
		expand("stats/{sample}.hg19.mkdup.sorted.bam.stats", sample=config["sample_list"]),
		expand("stats/{sample}.hg38.mkdup.sorted.bam.stats", sample=config["sample_list"])

rule strip_reads:
	input:
		bam = lambda wildcards: orig_bam_directory + config["sample_bams"][wildcards.sample]
	output:
		logfile = "xyalign/logfiles/{sample}_strip_reads_xyalign.log",
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq"
	params:
		xyalign = xyalign_path,
		sample_id = "{sample}_strip_reads",
		xmx = "16g",
		cpus = "4"
	shell:
		"source activate xyalign_env && python {params.xyalign} --STRIP_READS --ref null --bam {input.bam} --cpus {params.cpus} --xmx {params.xmx} --sample_id {params.sample_id} --output_dir xyalign --chromosomes ALL"

rule gzip_stripped_reads:
	# After stripping reads, they seemed to simply follow a "single read group
	# per sample, with the sample name as the readgroup name" format.  So,
	# should be able to simply gzip based on xyalign naming output fastqs
	# after the read groups, which happen to be sample names
	input:
		logfile = "xyalign/logfiles/{sample}_strip_reads_xyalign.log",
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq"
	output:
		gz_fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
		gz_fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz"
	run:
		shell("gzip {input.fq1}")
		shell("gzip {input.fq2}")

rule prepare_reference_hg19:
	input:
		hg19_ref_path
	output:
		fai = hg19_ref_path + ".fai",
		amb = hg19_ref_path + ".amb",
		dict = hg19_ref_prefix + ".dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule prepare_reference_hg38:
	input:
		hg38_ref_path
	output:
		fai = hg38_ref_path + ".fai",
		amb = hg38_ref_path + ".amb",
		dict = hg38_ref_prefix + ".dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

rule map_and_process_trimmed_reads_hg19:
	input:
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
		fai = hg19_ref_path + ".fai",
		ref = hg19_ref_path
	output:
		"processed_bams/{sample}.hg19.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools = samtools_path,
		samblaster = samblaster_path
	threads: 4
	shell:
		" {params.bwa} mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} "
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule map_and_process_trimmed_reads_hg38:
	input:
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
		fai = hg38_ref_path + ".fai",
		ref = hg38_ref_path
	output:
		"processed_bams/{sample}.hg38.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samtools =  samtools_path,
		samblaster = samblaster_path
	threads: 4
	shell:
		" {params.bwa} mem -t {threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} "
		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
		"-O bam -o {output}"

rule index_initial_bams:
	input:
		bam_19 = "processed_bams/{sample}.hg19.sorted.bam",
		bam_38 = "processed_bams/{sample}.hg38.sorted.bam"
	output:
		bam_19_idx = "processed_bams/{sample}.hg19.sorted.bam.bai",
		bam_38_idx = "processed_bams/{sample}.hg38.sorted.bam.bai"
	params:
		samtools = samtools_path
	run:
		shell("{params.samtools} index {input.bam_19}")
		shell("{params.samtools} index {input.bam_38}")

rule base_quality_recalibration_hg19_step1:
	input:
		bam = "processed_bams/{sample}.hg19.sorted.bam",
		bam_idx = "processed_bams/{sample}.hg19.sorted.bam.bai",
		dbsnp_gz = "misc/dbsnp_138.hg19.vcf.gz",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
		ref = hg19_ref_path
	output:
		recal = "stats/{sample}_recalibration_report_hg19.grp"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -o {output.recal} -knownSites {input.dbsnp_gz} -knownSites {input.mills_gz}"

rule base_quality_recalibration_hg19_step2:
	input:
		bam = "processed_bams/{sample}.hg19.sorted.bam",
		ref = hg19_ref_path,
		recal = "stats/{sample}_recalibration_report_hg19.grp"
	output:
		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"

rule base_quality_recalibration_hg38_step1:
	input:
		bam = "processed_bams/{sample}.hg38.sorted.bam",
		bam_idx = "processed_bams/{sample}.hg38.sorted.bam.bai",
		dbsnp_gz = "misc/dbsnp_138.hg38.vcf.gz",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
		ref = hg38_ref_path
	output:
		recal = "stats/{sample}_recalibration_report_hg38.grp"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -o {output.recal} -knownSites {input.dbsnp_gz} -knownSites {input.mills_gz}"

rule base_quality_recalibration_hg38_step2:
	input:
		bam = "processed_bams/{sample}.hg38.sorted.bam",
		ref = hg38_ref_path,
		recal = "stats/{sample}_recalibration_report_hg38.grp"
	output:
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"

rule indel_targetcreator_hg19:
	input:
		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
		ref = hg19_ref_path
	output:
		target_list = "stats/{sample}_target_list_hg19.txt"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -known {input.mills_gz}"

rule indel_realignment_hg19:
	input:
		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
		ref = hg19_ref_path,
		target_list = "stats/{sample}_target_list_hg19.txt"
	output:
		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.indelrealigned.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} -o {output.bam} -known {input.mills_gz} -targetIntervals {input.target_list}"

rule indel_targetcreator_hg38:
	input:
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
		ref = hg38_ref_path
	output:
		target_list = "stats/{sample}_target_list_hg38.txt"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -known {input.mills_gz}"

rule indel_realignment_hg38:
	input:
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam",
		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
		ref = hg38_ref_path,
		target_list = "stats/{sample}_target_list_hg38.txt"
	output:
		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.indelrealigned.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} -o {output.bam} -known {input.mills_gz} -targetIntervals {input.target_list}"

rule bam_stats:
	input:
		"processed_bams/{sample}.{chrom}.sorted.mkdup.recal.indelrealigned.bam"
	output:
		"stats/{sample}.{chrom}.mkdup.sorted.indel_realigned.bam.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"
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
