library("ggfortify")
library("ggplot2")
library("optparse")

option_list <- list(
	make_option(c("--working_directory"), action="store", type="character",
		help="Working Directory"),
	make_option(c("--chromosome"), action="store", type="character",
		help="Name of chromosome"),
	make_option(c("--input"), action="store", type="character",
		help="Input text file"),
	make_option(c("--output_pca"), action="store", type="character",
		help="Output file for PCA"),
	make_option(c("--output_hclust"), action="store", type="character",
		help="Output file for cluster analysis"),
	make_option(c("--sample_names"), action="store", type="character",
		help="File containing sample names"),
	make_option(c("--sample_types"), action="store", type="character",
		help="File containing sample types/locations. Must be in same order as sample names")
)

opt = parse_args(OptionParser(option_list=option_list))

# Processed parsed flags
setwd(opt$working_directory)

chrom <- opt$chromosome
input_file <- read.table(opt$input, sep="\t", header=
output_pca <- opt$output_pca
output_hclust <- opt$output_hclust
sample_names <- readLines(opt$sample_names)
sample_types <- readLines(opt$sample_types)

#if ((length(sample_names) == length(sample_types)) == FALSE) {stop("Length of sample names and lengths are not equal. Exiting")}

# Read input
input_table <- input_file[,-1]
rownames(input_table) <- sample_names
input_table$location <- sample_types
df2 <- input_table[, names(input_table) != "location"]

input.pca <- prcomp(df2, center = TRUE)

color_palette <- c("firebrick3","blue","black","green4")

pdf(output_pca)
pca_plot <- ggplot(input.pca,aes(x=PC1, y=PC2, color = input_table$location)) + geom_text(aes(label= rownames(input_table)), size =5) + ggtitle(paste(chrom, "PCA - PC1 and PC2", sep=" ") + scale_colour_manual(values=color_palette[1:length(unique(sample_types))]) + theme_bw()
dev.off()

dist_matrix <- dist(df2, method="euclidean")
fit_hclust <- hclust(dist_matrix, method="ward.D")

pdf(output_hclust)
plot(fit_hclust, main = chrom, xlab= "")
dev.off()
