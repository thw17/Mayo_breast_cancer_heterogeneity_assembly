# Breast_cancer_population_analysis

This repo contains code used in the initial stages of analyses for:

Phung et al. *Submitted*. Unique evolutionary trajectories of breast cancers with distinct genomic and spatial heterogeneity.

Specifically, read processing, read mapping, and BAM processing occur here, while code format
all downstream analyses can be found here: https://github.com/SexChrLab/Mayo_breast_regional_heterogeneity

## Quick note about running the pipeline

Note that there are two master rules in the Snakefile. The first, "rule all",
will run the steps described above, as well as some experimental variant calling with
Freebayes and generate some PCAs and cluster analyses. The second, "rule all_no_genotyping",
just runs the steps to produce the processed BAM files that are used in the downstream analyses
in the repo described above. The first option is default and what will run with any version of
Snakemake's default command lines. To run the second, use it as a target with a command like:

` $ snakemake all_no_genotyping `

## Steps in the pipeline

### 1. Strip reads
- The Mayo Clinic stores their data in BAM files, rather than as raw reads, so we
initially use `--STRIP_READS` in `XYalign` to strip the reads and sort them into
fastq files.

### 2. Prepare reference
- We used `XYalign` to initially prepare XX versions of the reference genome.
- We then used `BWA` and `samtools` to generate required indices.

### 3. Map reads
- We mapped reads with `BWA MEM`, adding read group info along the way.
- We used `samblaster` to mark duplicate reads as they were mapped.
- The output of `samblaster` is passed to `samtools` to fix mate pairing and sort the
mapped reads.
- We then index the resulting BAM file with `samtools`.

### 4. Base quality recalibration and indel realignment
- We used `GATK` for base quality recalibration and indel realignment following default
protocols.
- We used the Mills_and_1000G_gold_standard indels and dbSNP as reference panels for these steps.

#### These recalibrated, and indel realigned BAM files are then passed to the other pipeline (described above) for downstream analyses

------------

## Supplemental steps
If you run the default `Snakemake` command, the following additional steps are also run.

### s5. Freebayes variant calling
- Experimental run of `Freebayes` to apply "population" joint genotyping across samples.
- Uses parameters `--pooled-continuous --pooled-discrete -F 0.03 -C 2 --allele-balance-priors-off --genotype-qualities`.
- We zipped and indexed the output VCF file with `bgzip` and `tabix`, respectively.

### s6. BAM stats
- We used `samtools` to calculate BAM file statistics.

### s7. Filter VCF
- We wrote a custom Python script (`scripts/filter_vcf.py`) to process and filter the VCF file.
- We again used `bgzip` and `tabix` to zip and index this file.

### s8. Genotype comparisons, PCAs, and hierarchical cluster analyses
- In the final four rules, a couple of Python scripts and an R script to compare genotypes,
prepare genotypes for analyses, and run PCAs and cluster analyses
