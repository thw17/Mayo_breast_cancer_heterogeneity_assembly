from __future__ import print_function
import argparse
import cyvcf2
import numpy as np
import pandas as pd


def parse_args():
	parser = argparse.ArgumentParser(
		description="Compares genotypes among samples in a vcf.  Currently "
		"designed to handle Freebayes output")

	parser.add_argument(
		"--input_vcf", required=True,
		help="Bgzipped and tabix-indexed vcf file")

	parser.add_argument(
		"--output", required=True,
		help="Name of output file for distance matrix")

	parser.add_argument(
		"--qual", type=int, default=30,
		help="Minimum quality (int) for a site to be considered")

	parser.add_argument(
		"--mapq", type=int, default=20,
		help="Minimum mapq (int) for a site to be considered. Same value "
		"applied separately to MQM and MQMR")

	parser.add_argument(
		"--depth", type=int, default=10,
		help="Minimum depth (int) per sample for a site to be considered. "
		"*ALL* samples much be greater than or equal to this value for a site "
		"to be included.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	qual_threshold = args.qual
	mapq_threshold = args.mapq
	depth_threshold = args.depth

	vcf = cyvcf2.VCF(args.input_vcf)

	num_samples = len(vcf.samples)

	sample_array = np.zeros((num_samples, num_samples), dtype=np.int)

	counter = 0
	for variant in vcf:
		if variant.QUAL >= qual_threshold and variant.INFO('NS') == num_samples:
			if all([variant.INFO('MQM'), variant.INFO('MQM')]) >= mapq_threshold:
				dp_array = variant.format('DP')
				if len(dp_array[dp_array >= depth_threshold]):
					genotypes = variant.gt_types
					for idx, val in enumerate(genotypes):
						for k in range(idx + 1, len(genotypes)):
							if val != genotypes[k]:
								sample_array[val][k] += 1
		counter += 1
		if counter % 10000 == 0:
			print("{} records processed...".format(counter))

	df = pd.DataFrame(data=sample_array, index=vcf.samples, columns=vcf.samples)
	df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
	main()
