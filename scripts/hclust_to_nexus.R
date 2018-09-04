
# Command line should be:
# Rscript hclust_to_nexus.R output.nex input1 input2 input3 input4 ...
# where:
# output.nex is the name of your desired output file
# input1, input2, etc. are the files output by Process_vcf_output_table.py for each of your samples, chroms, etc.
# For both output and inputs, include full paths

library("ape")

args = commandArgs(trailingOnly=TRUE)

output_file <- args[1]
print(output_file)
input_files <- args[-1]
print(input_files)

phylo_list <- vector("list", length(input_files))

i <- 1
while(i <= length(input_files)) {
	in1 <- read.table(input_files[i], sep="\t", header=FALSE)
	input_table <- in1[, -1]
	rownames(input_table) <- in1[,1]
	dist_matrix <- dist(input_table, method="euclidean")
	fit_hclust <- hclust(dist_matrix, method="ward.D")
	phylo_list[[i]] <- as.phylo(fit_hclust)
	i <- i + 1
	print(i)
}

write.nexus(phylo_list, file = output_file)
