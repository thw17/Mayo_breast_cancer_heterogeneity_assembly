from __future__ import print_function
import argparse
import cyvcf2
import numpy as np


def parse_args():
	""" Parse command line arguments """

	parser = argparse.ArgumentParser(
		description="This program takes an input VCF file and outputs a table "
		"for downstream analysis. This program does not filter or split "
		"sites by chromosome. Of particular note is that its designed to "
		"take a VCF with no missing data.  Currently only supports Freebayes output")

	parser.add_argument(
		"--vcf", required=True,
		help="Input VCF file.  Must be bgzipped and tabix-indexed.")

	parser.add_argument(
		"--output", required=True,
		help="Output file. Will overwrite if exists.")

	parser.add_argument(
		"--normal_sample", default=None,
		help="The name of the sample (must match the name in the VCF exactly) "
		"to use to standardize allele descriptions (which is 'minor' and which "
		"is 'major').  For example, a control or germline sample. If blank, "
		"will use the first sample listed in the header as default.")

	args = parser.parse_args()


def main():
	args = parse_args()

	vcf = cyvcf2.VCF(args.vcf)

	samples = vcf.samples

	if args.normal_sample is None:
		normal_idx = 0
	else:
		try:
			normal_idx = np.where(samples == args.normal_sample)[0][0]
		except IndexError:
			normal_idx = 0

	geno_dict = {}
	for i in samples:
		geno_dict[i] = [i]

	balance_dict = {}
	for i in samples:
		balance_dict[i] = [i]

	for variant in vcf:
		gen_list = variant.genotype
		if len(set([tuple(x) for x in gen_list])) > 1 or (len(set([tuple(x) for x in gen_list])) == 1 and gen_list[0][0] != gen_list[0][1]):
			ref_depth_list = variant.format("RO")
			alt_depth_list = variant.format("AO")
			if alt_depth_list.shape[1] > 1:
				multiple_alts = True
			else:
				multiple_alts = False
	###################### Finish this
	##################################
			normal_genotype = gen_list[normal_idx][:2]
			a1, a2 = normal_genotype
			if a1 == 0:
				a1_freq = ref_depth_list[normal_idx]
			else:
				if multiple_alts is False:
					a1_freq = alt_depth_list[normal_idx]
				else:
					a1_freq = alt_depth_list[normal_idx][a1 - 1]

			if a2 == 0:
				a2_freq = ref_depth_list[normal_idx]
			else:
				if multiple_alts is False:
					a1_freq = alt_depth_list[normal_idx]
				else:
					a1_freq = alt_depth_list[normal_idx][a1 - 1]

			for idx, i in enumerate(samples):
				genotype = gen_list[idx]
				if genotype[:2] == [0, 0]:
					geno_dict[i].append(0)
				elif genotype[:2] == [0, 1]:
					geno_dict[i].append(1)
				elif genotype[:2] == [1, 1]:
					geno_dict[i].append(2)
				elif genotype[:2] == [1, 2]:
					geno_dict[i].append(3)
				elif genotype[:2] == [2, 2]:
					geno_dict[i].append(4)
				else:
					geno_dict[i].append(5)

# Output genotype table
output_list = []
for indv in samples:
	output_list.append(geno_dict[indv])
outfile = args.output
with open(outfile, "w") as f:
	w = csv.writer(f, dialect="excel-tab")
	w.writerows(output_list)
