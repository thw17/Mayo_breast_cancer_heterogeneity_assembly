library("ggfortify")
library("ggplot2")
library("optparse")

option_list <- list(
	make_option(c("-w", "--working_directory"), action="store", type="character",
		help="Working Directory"),
	make_option(c("-c", "--chromosome",), action="store", type="character",
		help="Name of chromosome"),
	make_option(c("-i", "--input"), action="store", type="character",
		help="Input text file"),
	make_option(c("-o", "--output"), action="store", type="character",
		help="Output file"),
	make_option(c("-n", "--sample_names"), action="store", type="character",
		help="File containing sample names"),
	make_option(c("-t", "--sample_types"), action="store", type="character",
		help="File containing sample types. Must be in same order as sample names")
)

opt = parse_args(OptionParser(option_list=option_list))

# Processed parsed flags
setwd(opt$working_directory)

chrom <- opt$chromosome
input_file <- read.table(opt$input, sep="\t", header=
output_file <- opt$output
sample_names <- readLines(opt$sample_names)
sample_types <- readLines(opt$sample_types)

#if ((length(sample_names) == length(sample_types)) == FALSE) {stop("Length of sample names and lengths are not equal. Exiting")}

# Read input
input_table <- input_file[,-1]
rownames(input_table) < 
