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
		"take a VCF with no missing data.")

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

	data_dict = {}
	for i in samples:
		data_dict[i] = [i]
	for variant in vcf:
		
