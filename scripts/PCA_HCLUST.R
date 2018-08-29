library("ggfortify")
library("ggplot2")
library("optparse")
library("dendextend")

option_list <- list(
	make_option(c("-w", "--working_directory"), action="store", type="character",
		help="Working Directory"),
	make_option(c("-c", "--chromosome"), action="store", type="character",
		help="Name of chromosome"),
	make_option(c("-i", "--input"), action="store", type="character",
		help="Input text file"),
	make_option(c("-p", "--output_pca"), action="store", type="character",
		help="Output file for PCA"),
	make_option(c("-d", "--output_hclust"), action="store", type="character",
		help="Output file for cluster analysis"),
	make_option(c("-n", "--sample_names"), action="store", type="character",
		help="File containing sample names"),
	make_option(c("-t", "--sample_types"), action="store", type="character",
		help="File containing sample types/locations. Must be in same order as sample names")
)

opt = parse_args(OptionParser(option_list=option_list))

# Processed parsed flags
setwd(opt$working_directory)

chrom <- opt$chromosome
input_file <- read.table(opt$input, sep="\t", header=FALSE)
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
ggplot(input.pca, aes(x=PC1, y=PC2, color = input_table$location)) + geom_text(aes(label= rownames(input_table)), size =5) + ggtitle(paste(chrom, "PCA - PC1 and PC2", sep=" ")) + scale_colour_manual(values=color_palette[1:length(unique(sample_types))]) + theme_bw()
dev.off()

dist_matrix <- dist(df2, method="euclidean")
fit_hclust <- hclust(dist_matrix, method="ward.D")
hcd <- as.dendrogram(fit_hclust)

pdf(output_hclust)
# Note "A2-P5" and "A2-P4" are only found in sample 1750, while "Normal" is found in both 1750 and 585
# "B1-P3", "B2-P3", "D1-P3", "F1-P3", and "F2-P3" are the adjacent and distant nodes in 585
# This code colors "Normal" blue and increases "Normal" line weight to 8 in both 1750 and 585
# In 1750, it will also make "A2-P5" and "A2-P4" red, increase their line weight to 4, and make "A2-P5" a dashed line
# In 585, it will increase all adjacent and distant node line weight to 4, change their color to red, and further make all distant nodes a dashed line
hcd %>%
	set("by_labels_branches_col", value = c("Normal"), TF_values = c("blue",1)) %>%
	set("by_labels_branches_col", value = c(
		"A2-P5", "A2-P4", "B1-P3", "B2-P3", "B3-P3", "D1-P3", "F1-P3", "F2-P3"), TF_values = c("red", Inf)) %>%
	set("by_labels_branches_lwd", value = c(
		"A2-P5", "A2-P4", "B1-P3", "B2-P3", "B3-P3", "D1-P3", "F1-P3", "F2-P3"), TF_values = c(4,Inf)) %>%
	set("by_labels_branches_lwd", value = c("Normal"), TF_values = c(8,Inf)) %>%
	set("by_labels_branches_lty", value = c("A2-P5", "B3-P3", "D1-P3", "F1-P3", "F2-P3"), TF_values = c(3,Inf)) %>%
	set("labels_cex", c(1.5)) %>% plot(main = chrom)
#plot(fit_hclust, main = chrom, xlab= "")
dev.off()
