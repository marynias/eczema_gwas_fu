library(tools)
library(readxl)
library(rio)

args = commandArgs(trailingOnly=TRUE)

my_file <- args[1]

# Read each file and write it to csv
my_path <- my_file
z <- lapply(excel_sheets(my_path), read_excel, path = my_path)
for (sheet in seq_len(length(z)))
{
output_file <- paste(file_path_sans_ext(my_file), ".", sheet , ".tsv", sep="")
write.table(z[[sheet]], output_file, row.names=FALSE, sep="\t", quote=FALSE)
}