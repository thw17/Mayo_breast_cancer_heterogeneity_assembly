from future import print_function
import argparse
import cyvcf2


def parse_args():
	""" Parse command line arguments """
	parser.argparse.ArgumentParser(
		description="This program filters VCF files based on user defined "
		"thresholds")

	parser.add_argument(
		"--vcf", required=True,
		help="Input VCF file.  Must be bgzipped and tabix-indexed.")

	parser.add_argument(
		"--output_vcf", required=True,
		help="Output VCF file.")

	parser.add_argument(
		"--variant_caller", type=str.lower(),
		default="freebayes", choices=["freebayes"],
		help="Variant caller used.  Currently only supports freebayes.  Default "
		"is freebayes.")

	parser.add_argument(
		"--QUAL", type=int, default=0
		help="Minimum phred-scaled *site* quality (QUAL column of vcf). "
		"Default is 0.")

	parser.add_argument(
		"--sample_depth", type=int, default=0,
		help="Minimum depth per sample required to retain site.  Default is 0.")

	parser.add_argument(
		"--type", type=str.upper, nargs=*, default=["ALL"], choices=[
			"ALL", "SNP", "INDEL", "MNP"],
		help="Type of variant to retain. Default is to retain all.  Choices "
		"include (can be either all uppcase or all lower case): all, snp, "
		"indel, mnp.  Multiple options can be selected, e.g. '--type snp indel'")

	parser.add_argument(
		"--min_samples", type=int, default=0,
		help="Minimum number of samples passing all filters required to retain "
		"sites. Integer. Default is 0.")

	args = parser.parse_args()

	return args


def main():
	args = parse_args()

	vcf = cyvcf2.VCF(args.input_vcf)

	vcf.add_to_header(
		"Filter_vcf_CMD=python Filter_vcf.py "
		"--vcf {} "
		"--output_vcf {} "
		"--variant_caller {}"
		"--min_samples {} "
		"--QUAL {} "
		"--sample_depth {}"
		"--type {} ".format(
			args.vcf,
			args.output_vcf,
			args.variant_caller,
			args.min_samples,
			args.QUAL,
			args.sample_depth,
			args.type))

	out_vcf = cyvcf2(vcf, args.output_vcf)

	for variant in vcf:
