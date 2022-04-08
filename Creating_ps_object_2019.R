library(qiime2R)
library(phyloseq)

# Importing each element
otus = qza_to_phyloseq("WB2019-table4.qza")
samdat = read.csv("WB2019-metadata.csv", header=TRUE)
taxtab = read_qza("WB2019-taxonomy.qza")
taxtab = parse_taxonomy(taxtab$data)

# Creating otu table
otus = otu_table(otus,taxa_are_rows = TRUE)

# Creating sample data
samdat2 = sample_data(samdat)
sample_names(samdat2)=samdat$X.SampleID

# Creating taxonomy table
taxtab2 = tax_table(taxtab)
taxa_names(taxtab2)=row.names(taxtab)

# Creating phyloseq object
ps = phyloseq(otus,samdat2,taxtab2)
ps

saveRDS(ps,"WB2019.ps")
