configfile: "breast_cancer_config.json"

normal_1750 = "PS13-1750-RightBreast-N-P3"
normal_585 = "PS13-585-Normal"

tumor_1750 = [x for x in config["PS13-1750"] if x != normal_1750]
tumor_585 = [x for x in config["PS13-585"] if x != normal_585]

# hg19_ref_path = "reference/ucsc.hg19.fasta"
# hg19_ref_prefix = "reference/ucsc.hg19"
# hg38_ref_path = "reference/Homo_sapiens_assembly38.fasta"
# hg38_ref_prefix = "reference/Homo_sapiens_assembly38"

# dbsnp_138_hg19_path = "misc/dbsnp_138.hg19.vcf"
# dbsnp_138_hg38_path = "misc/dbsnp_138.hg38.vcf"
# dpsnp_146_hg38_path = "misc/dbsnp_146.hg38.vcf"
# mills_1kg_indels_hg19 = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
# mills_1kg_indels_hg38 = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf"
# cosmic_coding_hg19 = "/home/thwebste/Data/COSMIC/CosmicCodingMuts_v81_hg19.vcf.gz"
# cosmic_noncoding_hg19 = "/home/thwebste/Data/COSMIC/CosmicNonCodingVariants_v81_hg19.vcf.gz"
# cosmic_coding_hg38 = "/home/thwebste/Data/COSMIC/CosmicCodingMuts_v81_hg38.vcf.gz"
# cosmic_noncoding_hg38 = "/home/thwebste/Data/COSMIC/CosmicNonCodingVariants_v81_hg38.vcf.gz"

orig_bam_directory = "/mnt/storage/SAYRES/MAYO/"
temp_dir_path = "temp/"

bwa_path = "bwa"
samtools_path = "samtools"
samblaster_path = "samblaster"
bgzip_path = "bgzip"
tabix_path = "tabix"
freebayes_path = "freebayes"
gatk_path = "/home/thwebste/Tools/GenomeAnalysisTK_37.jar"
mutect_path = "/home/thwebste/Tools/mutect-1.1.7.jar"
# picard_path =
xyalign_path = "/scratch/thwebste/xyalign_test/XYalign/xyalign/xyalign.py"
xyalign_env_name = "breast_cancer_xyalign"

rule all:
	input:
		expand(
			"xyalign/logfiles/{sample}_strip_reads_xyalign.log",
			sample=config["sample_list"]),
		expand(
			"xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
			sample=config["sample_list"]),
		expand(
			"xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
			sample=config["sample_list"]),
		expand(
			"stats/{individual}.{chrom}.{assembly}.compare_genotypes.txt",
			individual=config["individuals"], chrom=config["chromosomes"],
			assembly=config["reference_list"]),
		expand(
			"stats/{sample}.{assembly}.mkdup.sorted.indelrealigned.bam.stats",
			sample=config["sample_list"], assembly=config["reference_list"])
		# expand(
		# 	"processed_bams/{sample}.hg38.sorted.bam",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"processed_bams/{sample}.hg19.sorted.bam",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"processed_bams/{sample}.hg19.sorted.mkdup.recal.indelrealigned.bam",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"processed_bams/{sample}.hg38.sorted.mkdup.recal.indelrealigned.bam",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"stats/{sample}.hg19.mkdup.sorted.indel_realigned.bam.stats",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"stats/{sample}.hg38.mkdup.sorted.indel_realigned.bam.stats",
		# 	sample=config["sample_list"]),
		# expand(
		# 	"calls/PS13-1750.{chrom}.hg19.raw.vcf.gz.tbi",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"calls/PS13-585.{chrom}.hg19.raw.vcf.gz.tbi",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"calls/PS13-1750.{chrom}.hg38.raw.vcf.gz.tbi",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"calls/PS13-585.{chrom}.hg38.raw.vcf.gz.tbi",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"stats/PS13-1750.{chrom}.hg19.compare_genotypes.txt",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"stats/PS13-585.{chrom}.hg19.compare_genotypes.txt",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"stats/PS13-1750.{chrom}.hg38.compare_genotypes.txt",
		# 	chrom=config["chromosomes"]),
		# expand(
		# 	"stats/PS13-585.{chrom}.hg38.compare_genotypes.txt",
		# 	chrom=config["chromosomes"]),
		# expand("vcf/{sample}.585.hg19.mutect2.raw.vcf.gz", sample=tumor_585)

rule strip_reads:
	input:
		bam = lambda wildcards: orig_bam_directory + config["sample_bams"][wildcards.sample]
	output:
		logfile = "xyalign/logfiles/{sample}_strip_reads_xyalign.log",
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq"
	params:
		xyalign = xyalign_path,
		xyalign_env = xyalign_env_name,
		sample_id = "{sample}_strip_reads",
		xmx = "16g",
		cpus = "4"
	shell:
		"source activate {params.xyalign_env} && python {params.xyalign} --STRIP_READS --ref null --bam {input.bam} --cpus {params.cpus} --xmx {params.xmx} --sample_id {params.sample_id} --output_dir xyalign --chromosomes ALL"

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

rule xyalign_create_xx_reference:
	input:
		ref = lambda wildcards: config["reference_genomes"][wildcards.assembly]
	output:
		"xyalign/reference/{assembly}.XXonly.fasta"
	params:
		xyalign = xyalign_path,
		xyalign_env = xyalign_env_name,
		gen_assembly = "{assembly}"
	threads:
		4
	shell:
		"source activate {params.xyalign_env} && python {params.xyalign} --PREPARE_REFERENCE --ref {input.ref} --bam null --xx_ref_out {params.gen_assembly}.XXonly.fasta --output_dir xyalign --x_chromosome chrX --y_chromosome chrY"

rule prepare_reference:
	input:
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		fai = "xyalign/reference/{assembly}.XXonly.fasta.fai",
		amb = "xyalign/reference/{assembly}.XXonly.fasta.amb",
		dict = "xyalign/reference/{assembly}.XXonly.dict"
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

rule map_and_process_trimmed_reads:
	input:
		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
		fai = "xyalign/reference/{assembly}.XXonly.fasta.fai",
		amb = "xyalign/reference/{assembly}.XXonly.fasta.amb",
		dict = "xyalign/reference/{assembly}.XXonly.dict",
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		"processed_bams/{sample}.{assembly}.sorted.bam"
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

rule index_initial_bams:
	input:
		bam = "processed_bams/{sample}.{assembly}.sorted.bam",
	output:
		bam_idx = "processed_bams/{sample}.{assembly}.sorted.bam.bai",

	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input.bam}"

rule base_quality_recalibration_step1:
	input:
		bam = "processed_bams/{sample}.{assembly}.sorted.bam",
		bam_idx = "processed_bams/{sample}.{assembly}.sorted.bam.bai",
		dbsnp_gz = lambda wildcards: config["dbsnp"][wildcards.assembly],
		mills_gz = lambda wildcards: config["mills"][wildcards.assembly],
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		recal = "stats/{sample}_recalibration_report_{assembly}.grp"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -o {output.recal} -knownSites {input.dbsnp_gz} -knownSites {input.mills_gz}"

rule base_quality_recalibration_step2:
	input:
		bam = "processed_bams/{sample}.{assembly}.sorted.bam",
		ref = "xyalign/reference/{assembly}.XXonly.fasta",
		recal = "stats/{sample}_recalibration_report_{assembly}.grp"
	output:
		bam = "processed_bams/{sample}.{assembly}.sorted.mkdup.recal.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"

rule indel_targetcreator:
	input:
		bam = "processed_bams/{sample}.{assembly}.sorted.mkdup.recal.bam",
		mills_gz = lambda wildcards: config["mills"][wildcards.assembly],
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		target_list = "stats/{sample}_target_list_{assembly}.list"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -known {input.mills_gz}"

rule indel_realignment:
	input:
		bam = "processed_bams/{sample}.{assembly}.sorted.mkdup.recal.bam",
		mills_gz = lambda wildcards: config["mills"][wildcards.assembly],
		ref = "xyalign/reference/{assembly}.XXonly.fasta",
		target_list = "stats/{sample}_target_list_{assembly}.list"
	output:
		bam = "processed_bams/{sample}.{assembly}.sorted.mkdup.recal.indelrealigned.bam"
	params:
		gatk = gatk_path,
		temp_dir = temp_dir_path,
		java_mem = "-Xmx16g"
	shell:
		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} -o {output.bam} -known {input.mills_gz} -targetIntervals {input.target_list}"

rule freebayes_call_single_chrom:
	input:
		bam = lambda wildcards: expand("processed_bams/{sample}.{genome}.sorted.mkdup.recal.indelrealigned.bam", sample=config[wildcards.individual], genome=wildcards.assembly),
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		"calls/{individual}.{chrom}.{assembly}.raw.vcf"
	params:
		region = "{chrom}",
		freebayes = freebayes_path
	threads: 4
	shell:
		"{params.freebayes} -f {input.ref} --region {params.region} --pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities {input.bam} > {output}"

rule bam_stats:
	input:
		"processed_bams/{sample}.{assembly}.sorted.mkdup.recal.indelrealigned.bam"
	output:
		"stats/{sample}.{assembly}.mkdup.sorted.indelrealigned.bam.stats"
	shell:
		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"

rule zip_cosmic_vcf:
	input:
		vcf = "misc/Cosmic{type}_v81_{genome}.vcf"
	output:
		"misc/Cosmic{type}_v81_{genome}.vcf.gz"
	params:
		bgzip = bgzip_path
	shell:
		"{params.bgzip} {input.vcf}"

rule index_zipped_cosmic_vcf:
	input:
		vcf = "misc/Cosmic{type}_v81_{genome}.vcf.gz"
	output:
		"misc/Cosmic{type}_v81_{genome}.vcf.gz.tbi"
	params:
		tabix = tabix_path
	shell:
		"{params.tabix} -p vcf {input}"

rule combine_cosmic_vcfs_hg19:
	input:
		noncoding = "misc/CosmicNonCodingVariants_v81_{assembly}.vcf.gz",
		noncoding_idx = "misc/CosmicNonCodingVariants_v81_{assembly}.vcf.gz.tbi",
		coding = "misc/CosmicCodingMuts_v81_{assembly}.vcf.gz",
		coding_idx = "misc/CosmicCodingMuts_v81_{assembly}.vcf.gz.tbi",
		ref = "xyalign/reference/{assembly}.XXonly.fasta"
	output:
		"misc/cosmic_{assembly}_combined.vcf.gz"
	params:
		temp_dir = temp_dir_path,
		gatk = gatk_path
	shell:
		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T CombineVariants -R {input.ref} --variant {input.coding} --variant {input.noncoding} -o {output}"

rule zip_vcfs:
	input:
		vcf = "calls/{individual}.{chrom}.{assembly}.raw.vcf"
	output:
		"calls/{individual}.{chrom}.{assembly}.raw.vcf.gz"
	shell:
		"bgzip {input.vcf}"

rule index_zipped_vcf:
	input:
		"calls/{individual}.{chrom}.{assembly}.raw.vcf.gz"
	output:
		"calls/{individual}.{chrom}.{assembly}.raw.vcf.gz.tbi"
	shell:
		"tabix -p vcf {input}"

rule compare_genotypes:
	input:
		vcf = "calls/{individual}.{chrom}.{assembly}.raw.vcf.gz",
		idx = "calls/{individual}.{chrom}.{assembly}.raw.vcf.gz.tbi"
	output:
		"stats/{individual}.{chrom}.{assembly}.compare_genotypes.txt"
	shell:
		"python scripts/Compare_genotypes.py --input_vcf {input.vcf} --output {output} --mapq 0"

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

# rule map_and_process_trimmed_reads_hg19:
# 	input:
# 		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
# 		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
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
# 		samtools = samtools_path,
# 		samblaster = samblaster_path
# 	threads: 4
# 	shell:
# 		" {params.bwa} mem -t {threads} -R "
# 		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
# 		"{input.ref} {input.fq1} {input.fq2}"
# 		"| {params.samblaster} "
# 		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
# 		"-O bam -o {output}"
#
# rule map_and_process_trimmed_reads_hg38:
# 	input:
# 		fq1 = "xyalign/fastq/{sample}_strip_reads_{sample}_1.fastq.gz",
# 		fq2 = "xyalign/fastq/{sample}_strip_reads_{sample}_2.fastq.gz",
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
# 		samtools =  samtools_path,
# 		samblaster = samblaster_path
# 	threads: 4
# 	shell:
# 		" {params.bwa} mem -t {threads} -R "
# 		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
# 		"{input.ref} {input.fq1} {input.fq2}"
# 		"| {params.samblaster} "
# 		"| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
# 		"-O bam -o {output}"
#
# rule index_initial_bams:
# 	input:
# 		bam_19 = "processed_bams/{sample}.hg19.sorted.bam",
# 		bam_38 = "processed_bams/{sample}.hg38.sorted.bam"
# 	output:
# 		bam_19_idx = "processed_bams/{sample}.hg19.sorted.bam.bai",
# 		bam_38_idx = "processed_bams/{sample}.hg38.sorted.bam.bai"
# 	params:
# 		samtools = samtools_path
# 	run:
# 		shell("{params.samtools} index {input.bam_19}")
# 		shell("{params.samtools} index {input.bam_38}")
#
# rule base_quality_recalibration_hg19_step1:
# 	input:
# 		bam = "processed_bams/{sample}.hg19.sorted.bam",
# 		bam_idx = "processed_bams/{sample}.hg19.sorted.bam.bai",
# 		dbsnp_gz = "misc/dbsnp_138.hg19.vcf.gz",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
# 		ref = hg19_ref_path
# 	output:
# 		recal = "stats/{sample}_recalibration_report_hg19.grp"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -o {output.recal} -knownSites {input.dbsnp_gz} -knownSites {input.mills_gz}"
#
# rule base_quality_recalibration_hg19_step2:
# 	input:
# 		bam = "processed_bams/{sample}.hg19.sorted.bam",
# 		ref = hg19_ref_path,
# 		recal = "stats/{sample}_recalibration_report_hg19.grp"
# 	output:
# 		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"
#
# rule base_quality_recalibration_hg38_step1:
# 	input:
# 		bam = "processed_bams/{sample}.hg38.sorted.bam",
# 		bam_idx = "processed_bams/{sample}.hg38.sorted.bam.bai",
# 		dbsnp_gz = "misc/dbsnp_138.hg38.vcf.gz",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
# 		ref = hg38_ref_path
# 	output:
# 		recal = "stats/{sample}_recalibration_report_hg38.grp"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T BaseRecalibrator -R {input.ref} -I {input.bam} -o {output.recal} -knownSites {input.dbsnp_gz} -knownSites {input.mills_gz}"
#
# rule base_quality_recalibration_hg38_step2:
# 	input:
# 		bam = "processed_bams/{sample}.hg38.sorted.bam",
# 		ref = hg38_ref_path,
# 		recal = "stats/{sample}_recalibration_report_hg38.grp"
# 	output:
# 		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T PrintReads -R {input.ref} -I {input.bam} -BQSR {input.recal} -o {output.bam}"
#
# rule indel_targetcreator_hg19:
# 	input:
# 		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
# 		ref = hg19_ref_path
# 	output:
# 		target_list = "stats/{sample}_target_list_hg19.list"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -known {input.mills_gz}"
#
# rule indel_realignment_hg19:
# 	input:
# 		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.bam",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz",
# 		ref = hg19_ref_path,
# 		target_list = "stats/{sample}_target_list_hg19.list"
# 	output:
# 		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.indelrealigned.bam"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} -o {output.bam} -known {input.mills_gz} -targetIntervals {input.target_list}"
#
# rule indel_targetcreator_hg38:
# 	input:
# 		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
# 		ref = hg38_ref_path
# 	output:
# 		target_list = "stats/{sample}_target_list_hg38.list"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T RealignerTargetCreator -R {input.ref} -I {input.bam} -o {output.target_list} -known {input.mills_gz}"
#
# rule indel_realignment_hg38:
# 	input:
# 		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.bam",
# 		mills_gz = "misc/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
# 		ref = hg38_ref_path,
# 		target_list = "stats/{sample}_target_list_hg38.list"
# 	output:
# 		bam = "processed_bams/{sample}.hg38.sorted.mkdup.recal.indelrealigned.bam"
# 	params:
# 		gatk = gatk_path,
# 		temp_dir = temp_dir_path,
# 		java_mem = "-Xmx16g"
# 	shell:
# 		"java {params.java_mem} -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T IndelRealigner -R {input.ref} -I {input.bam} -o {output.bam} -known {input.mills_gz} -targetIntervals {input.target_list}"
#
# rule bam_stats:
# 	input:
# 		"processed_bams/{sample}.{genome}.sorted.mkdup.recal.indelrealigned.bam"
# 	output:
# 		"stats/{sample}.{genome}.mkdup.sorted.indel_realigned.bam.stats"
# 	shell:
# 		"samtools stats {input} | grep ^SN | cut -f 2- > {output}"
#
# rule freebayes_call_single_chrom_hg19:
# 	input:
# 		bam = lambda wildcards: expand("processed_bams/{sample}.hg19.sorted.mkdup.recal.indelrealigned.bam", sample=config[wildcards.individual]),
# 		ref = hg19_ref_path
# 	output:
# 		"calls/{individual}.{chrom}.hg19.raw.vcf"
# 	params:
# 		region = "{chrom}",
# 		freebayes = freebayes_path
# 	threads: 4
# 	shell:
# 		"{params.freebayes} -f {input.ref} --region {params.region} --pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities {input.bam} > {output}"
#
# rule freebayes_call_single_chrom_hg38:
# 	input:
# 		bam = lambda wildcards: expand("processed_bams/{sample}.hg38.sorted.mkdup.recal.indelrealigned.bam", sample=config[wildcards.individual]),
# 		ref = hg38_ref_path
# 	output:
# 		"calls/{individual}.{chrom}.hg38.raw.vcf"
# 	params:
# 		region = "{chrom}",
# 		freebayes = freebayes_path
# 	threads: 4
# 	shell:
# 		"{params.freebayes} -f {input.ref} --region {params.region} --pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities {input.bam} > {output}"
#
# rule zip_cosmic_vcf:
# 	input:
# 		vcf = "misc/Cosmic{type}_v81_{genome}.vcf"
# 	output:
# 		"misc/Cosmic{type}_v81_{genome}.vcf.gz"
# 	params:
# 		bgzip = bgzip_path
# 	shell:
# 		"{params.bgzip} {input.vcf}"
#
# rule index_zipped_cosmic_vcf:
# 	input:
# 		vcf = "misc/Cosmic{type}_v81_{genome}.vcf.gz"
# 	output:
# 		"misc/Cosmic{type}_v81_{genome}.vcf.gz.tbi"
# 	params:
# 		tabix = tabix_path
# 	shell:
# 		"{params.tabix} -p vcf {input}"
#
# rule combine_cosmic_vcfs_hg19:
# 	input:
# 		noncoding = "misc/CosmicNonCodingVariants_v81_hg19.vcf.gz",
# 		noncoding_idx = "misc/CosmicNonCodingVariants_v81_hg19.vcf.gz.tbi",
# 		coding = "misc/CosmicCodingMuts_v81_hg19.vcf.gz",
# 		coding_idx = "misc/CosmicCodingMuts_v81_hg19.vcf.gz.tbi",
# 		ref = hg19_ref_path
# 	output:
# 		"misc/cosmic_hg19_combined.vcf.gz"
# 	params:
# 		temp_dir = temp_dir_path,
# 		gatk = gatk_path
# 	shell:
# 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T CombineVariants -R {input.ref} --variant {input.coding} --variant {input.noncoding} -o {output}"
#
#
# # rule mutect2_single_chrom_585_hg19:
# # 	input:
# # 		bam = "processed_bams/{sample}.hg19.sorted.mkdup.recal.indelrealigned.bam",
# # 		normal = "processed_bams/PS13-585-Normal.hg19.sorted.mkdup.recal.indelrealigned.bam",
# # 		ref = hg19_ref_path,
# # 		dbsnp_gz = "misc/dbsnp_138.hg19.vcf.gz",
# # 		cosmic = "misc/cosmic_hg19_combined.vcf.gz"
# # 	output:
# # 		"vcf/{sample}.585.hg19.mutect2.raw.vcf.gz"
# # 	params:
# # 		temp_dir = temp_dir_path,
# # 		gatk = gatk_path
# # 	shell:
# # 		"java -Xmx12g -Djava.io.tmpdir={params.temp_dir} -jar {params.gatk} -T MuTect2 -R {input.ref} -I:tumor {input.bam} -I:normal {input.normal} --dbsnp {input.dbsnp_gz} --cosmic {input.cosmic} -o {output}"
#
# # rule mutect2_single_chrom_1750:
#
# rule zip_vcfs:
# 	input:
# 		vcf = "calls/{individual}.{chrom}.{genome}.raw.vcf"
# 	output:
# 		"calls/{individual}.{chrom}.{genome}.raw.vcf.gz"
# 	shell:
# 		"bgzip {input.vcf}"
#
# rule index_zipped_vcf:
# 	input:
# 		"calls/{individual}.{chrom}.{genome}.raw.vcf.gz"
# 	output:
# 		"calls/{individual}.{chrom}.{genome}.raw.vcf.gz.tbi"
# 	shell:
# 		"tabix -p vcf {input}"
#
# rule compare_genotypes:
# 	input:
# 		vcf = "calls/{individual}.{chrom}.{genome}.raw.vcf.gz",
# 		idx = "calls/{individual}.{chrom}.{genome}.raw.vcf.gz.tbi"
# 	output:
# 		"stats/{individual}.{chrom}.{genome}.compare_genotypes.txt"
# 	shell:
# 		"python scripts/Compare_genotypes.py --input_vcf {input.vcf} --output {output} --mapq 0"
