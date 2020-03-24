library(biomaRt)
mart <- useMart("ensembl")
mart <- useDataset("hsapiens_gene_ensembl",mart)
string <- read.table("9606.protein.links.full.v10.5.conf.txt", header=T, stringsAsFactors=FALSE)
res1 <- getBM(attributes = c('ensembl_transcript_id', 
                            'ensembl_gene_id', 
                            'ensembl_peptide_id', 
                            'hgnc_symbol',
                            'external_gene_name',
                            'description',
                            'start_position',
                            'end_position',
                            'transcript_start',
                            'transcript_end'),
             filters = 'ensembl_peptide_id', 
             values = string$protein1,
             mart = mart)

res2 <- getBM(attributes = c('ensembl_transcript_id', 
                             'ensembl_gene_id', 
                             'ensembl_peptide_id', 
                             'hgnc_symbol',
                             'external_gene_name',
                             'description',
                             'start_position',
                             'end_position',
                             'transcript_start',
                             'transcript_end'),
              filters = 'ensembl_peptide_id', 
              values = string$protein2,
              mart = mart)

gene_name1 <- res1$hgnc_symbol[match(string$protein1, res1$ensembl_peptide_id)]
string$protein1_hugo <- gene_name1 

gene_name2 <- res2$hgnc_symbol[match(string$protein2, res2$ensembl_peptide_id)]
string$protein2_hugo <- gene_name2

#Remove pairs if either does not have a gene name assigned.
my_keep <- string[string$protein1_hugo!="" & string$protein2_hugo!="",]
my_keep <- my_keep[!is.na(string$protein1_hugo),]
my_keep <- my_keep[!is.na(string$protein2_hugo),]

write.table(my_keep, file="9606.protein.links.full.v10.5.conf.ensembl", quote=F, sep="\t", col.names=T)