# Load packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(DECIPHER)
library(Biostrings)

# Set working directory
setwd("~/Box/WhitmanLabMaster/WhitmanLab/Projects/WoodBuffalo/NWT-Fire-Microbes/WB2019/code")

# Import 2019 dataset
ps.5 = readRDS("../data/Seqs-processing/WB2019.ps")
ps.5

# Remove Blanks and multicommunity samples
# May want to remove U06-M because we don't have information on it)
# First checking names
sample_names(ps.5)
# Making list of keeper samples
# Note: WB19-38-B was the new site established because the previous site had been logged between sampling years. (15S-NT-38)
# If we check the ordinations, the M and O horizons for the B and non-B sites are similar to each other.
# We'll use the B site going forward, since the first site has likely been disturbed
# We also have WB19-27-O-1 and WB19-27-O-2.
# It seems like WB19-27-O-2 is the correct sample for that site.
SamplesToKeep = c("WB19-01-M", "WB19-01-O",  "WB19-02-O", "WB19-04-M",  "WB19-05-M", "WB19-05-O", "WB19-06-M", "WB19-06-O", "WB19-07-M",  
                  "WB19-07-O",  "WB19-08-M", "WB19-08-O", "WB19-09-M",  "WB19-09-O", "WB19-10-M", "WB19-10-O", "WB19-11-M", "WB19-11-O", "WB19-12-M",  
                  "WB19-12-O", "WB19-13-O", "WB19-14-M", "WB19-14-O", "WB19-15-M", "WB19-15-O", "WB19-16-M", "WB19-16-O", "WB19-17-M", "WB19-17-O", 
                  "WB19-18-M", "WB19-19-O", "WB19-22-M", "WB19-22-O", "WB19-23-M", "WB19-23-O", "WB19-27-M","WB19-27-O-2", "WB19-28-M",  
                  "WB19-28-O", "WB19-30-M", "WB19-30-O", "WB19-33-M", "WB19-33-O", "WB19-34-M", "WB19-34-O", "WB19-36-O",  "WB19-38-B-M", "WB19-38-B-O",
                  "WB19-39-M", "WB19-39-O", "WB19-40-M",  "WB19-40-O",  "WB19-41-M", "WB19-41-O", "WB19-42-M",  "WB19-42-O", 
                  "WB19-U01-M", "WB19-U01-O", "WB19-U03-M", "WB19-U03-O", "WB19-U04-M", "WB19-U04-O", "WB19-U05-M", "WB19-U06-M",
                  "WB19-U05-O", "WB19-U06-M", "WB19-U06-O", "WB19-U07-O", "WB19-U08-M", "WB19-U08-O", "WB19-U09-M", "WB19-U09-O", "WB19-U10-O") 

# Not keeping blanks and multicommunity
ps.5.pruned =  prune_samples(SamplesToKeep, ps.5)

# Remove chloroplasts and mitochondria
colnames(tax_table(ps.5.pruned))=c("Domain","Phylum","Class","Order","Family","Genus","Species")
ps.5.pruned = subset_taxa(ps.5.pruned, Genus != "Chloroplast")
ps.5.pruned = subset_taxa(ps.5.pruned, Genus != "Mitochondria")
ps.5 = ps.5.pruned

# Bringing in the 1 year post fire dataset
ps.1 = readRDS("../../WB2015/data/Seq_data/16S/CHTC/Dada2_Results_Full/ps.merged")
# Keep only sites that are in the 2019 dataset
ps.1 = subset_samples(ps.1,Site_ID %in% sample_data(ps.5.pruned)$Site_ID)
# Remove two unburned mineral samples that we only collected in 2015 (but do have O horizons for)
ps.1 =  subset_samples(ps.1, !(Sample_ID %in% c("15S-NT-U07M","15S-NT-U10M")))
# Remove mitochondria
ps.1 = subset_taxa(ps.1, Family != "Mitochondria")
# Make sure all 2019 sites are also in the 2015 sites
ps.5 = subset_samples(ps.5,Site_ID %in% sample_data(ps.1)$Site_ID)

# Make sure no zero-count samples are present
ps.5 = subset_taxa(ps.5, taxa_sums(ps.5)>0)
ps.1 = subset_taxa(ps.1, taxa_sums(ps.1)>0)

# Formatting severity
sample_data(ps.5)$Severity_Class = sample_data(ps.5)$Severity_Class = factor(sample_data(ps.5)$Severity_Class,
                                                                             levels = c("High", "Moderate", "Low", 'Unburned'), ordered = TRUE)
# Fix false zeros (NAs)
sample_data(ps.5.pruned)$pH[sample_data(ps.5.pruned)$pH==0]=NA

# Add year burned
sample_data(ps.1)$Years_Since_Fire = "1"
sample_data(ps.5)$Years_Since_Fire = "5"

############## First approach is to take the fasta files and merge based on identical sequences

# Because the sequenced regions were identical between the two years,
# we should be able merge all OTUs by their sequence alone

# Bring in the DNA sequences
dna.5 = readDNAStringSet("../data/Seqs-processing/Exported_Data/WB2019-rep-seqs4.fasta")
dna.1 = readDNAStringSet("../../WB2015/data/Seq_data/16S/CHTC/Dada2_Results_Full/DADA2_seqs_nochim.fasta")

# Extract DNA sequences
seqString.5 = c()
for (i in 1:length(dna.5)){
  seq = paste(dna.5[[i]])
  seqString.5[i]=seq
}
head(seqString.5)
length(seqString.5)

# Make table with ID and fasta
df.5 = data.frame(OTU=dna.5@ranges@NAMES,Sequence=seqString.5)
dim(df.5)
head(df.5)

# Subset it so only OTUs that remain in the phyloseq object are present
df.5 = df.5[df.5$OTU %in% taxa_names(ps.5),]
dim(df.5)

# Let's rename our OTUs.
ps.5.fasta = ps.5
sum(taxa_names(ps.5.fasta)==df.5$OTU)

# They are in the same order, so we can do the following
taxa_names(ps.5.fasta)=df.5$Sequence

### Follow the same procedure for 2015 dataset

# Extract DNA sequences
seqString.1 = c()
for (i in 1:length(dna.1)){
  seq = paste(dna.1[[i]])
  seqString.1[i]=seq
}
length(seqString.1)
# Have a string of the DNA sequences

# Make table with ID and fasta
df.1 = data.frame(OTU=dna.1@ranges@NAMES,Sequence=seqString.1)
dim(df.1)
head(df.1)

# 2015 dataset had number of sequences in addition to taxa names
# Need to separate these to they match our phyloseq object (names only)
m = data.frame(matrix(unlist(strsplit(as.character(df.1$OTU),";")),ncol=2,byrow=T))
colnames(m)=c("OTU","size")
head(m)

df.1$OTU = m$OTU
head(df.1$OTU)

# Subset it so only OTUs that remain in the phyloseq object are present
df.1 = df.1[df.1$OTU %in% taxa_names(ps.1),]
dim(df.1)

# Let's rename our OTUs.
ps.1.fasta = ps.1
sum(taxa_names(ps.1.fasta)==df.1$OTU)

# They are in the same order, so we can do the following
taxa_names(ps.1.fasta)=df.1$Sequence

### Now we have two phyloseq objects where the names of the OTUs are just their DNA sequences.

# Ok, we now should have taxa names are the fasta sequence for each phyloseq object.
# Now we want to merge them.
ps.merged.fasta = merge_phyloseq(ps.1.fasta,ps.5.fasta)
ps.merged.fasta
ps.1.fasta
ps.5.fasta
# All the samples made it through.

# How many OTUs were there in the two datasets alone?
ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)

# After merging, how many OTUs were there?
ntaxa(ps.merged.fasta)
# There are 25,788 unique OTUs in the combined datasets

# What was the difference?
ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)-ntaxa(ps.merged.fasta)

# So, 6,597 OTUs were present in the separate datasets
# So, the number of shared OTUs in the combined dataset should be half that
(ntaxa(ps.1.fasta)+ntaxa(ps.5.fasta)-ntaxa(ps.merged.fasta))/2
# Roughly 3,299 shared OTUs.
# Seems pretty good - other taxa likely rare.

################ Now combining them based on alignment and clustering
## This approach accounts for if, for some reason, there was slightly different trimming, etc.

# Running clustering on full fasta file
# Set processors
nproc=4

# Join the DNA sequences
dna = c(dna.1,dna.5)

### Using a chained guide tree can help speed up the alignment
# Recommended in the DECIPHER manual
# Form a chained guide tree

gT <- lapply(order(width(dna), decreasing=TRUE),
             function(x) {
               attr(x, "height") <- 0
               attr(x, "label") <- names(dna)[x]
               attr(x, "members") <- 1L
               attr(x, "leaf") <- TRUE
               x
             })

attr(gT, "height") <- 0.5
attr(gT, "members") <- length(dna)
class(gT) <- "dendrogram"

# Use the guide tree as input for alignment
# Align the sequences
# Started at 2:05 PM - finished at midnight

#aln <- DECIPHER::AlignSeqs(dna, processors = nproc, guideTree = gT, iterations = 0, refinements = 0)
##saveRDS(aln,"alignment.rds")
#aln = readRDS("alignment.rds")

# Find similarity between all seqs
# By ensuring includeTerminalGaps is false, we can be sure that
# differences in sequence length alone are not counted.

#d <- DECIPHER::DistanceMatrix(aln, includeTerminalGaps=FALSE, processors = nproc, type="dist")

# Cluster at identical similarity
#clusters <- DECIPHER::IdClusters(
#  d, 
#  method = "single",
#  cutoff = 0.0, # must be identical (except ends); could use `cutoff = 0.03` for a 97% OTU 
#  processors = nproc)

#saveRDS(clusters,"clusters.rds")
clusters = readRDS("clusters.rds")
head(clusters)
length(unique(clusters$cluster))
dim(clusters)
# 33183 unique OTUs even after merging
# This number still includes pruned taxa, though (e.g., mitochondria, multicomm)

# Let's see how common the clusters are
clust.num = clusters%>%
  dplyr::group_by(cluster)%>%
  dplyr::summarize(TotalOTUs=n())%>%
  dplyr::arrange(-TotalOTUs)
head(clust.num)
# There are only two clusters that have more than 2 OTUs (3).
# Doesn't seem like an issue - probably a minor end-of seq difference.
# It actually turns out they were mitochondria, so won't make it into the phyloseq objects
# We can use this, along with the subsequent code, to merge the OTUs.

## Use dplyr to merge the columns of the seqtab matrix for ASVs in the same OTU
# prep by adding sequences to the `clusters` data frame
clusters <- clusters %>%
  dplyr::mutate(OTUID = row.names(clusters))

head(clusters)
tail(clusters)
# We do need to split out the sq3;size= format from the first dataset

# Where the two datasets split
clusters[20020:20025,]

# Create temporary naming cluster
m = clusters[1:20020,]

# pull out the OTU ids as they exist in the ps object
m = data.frame(matrix(unlist(strsplit(as.character(m$OTUID),";")),ncol=2,byrow=T))
colnames(m)=c("OTU","size")
head(m)

# Make proper name column (OTU)
clusters$OTU = clusters$OTUID
clusters$OTU[1:20020]=paste(m$OTU)
head(clusters)
tail(clusters)
length(levels(as.factor(clusters$cluster)))

# Make separate naming dataframes for each year
clusters.1 = clusters[1:20020,]
clusters.5 = clusters[20021:dim(clusters)[1],]

# All the OTU names for the cluster IDs are in the phyloseq objects.
ps.1
sum(taxa_names(ps.1) %in% clusters.1$OTU)

ps.5
sum(taxa_names(ps.5) %in% clusters.5$OTU)
# All are present.

# Now we need to rename the OTUs to match each other between datasets
names.5 = data.frame(OTU=taxa_names(ps.5))
head(names.5)
names.5 =  plyr::join(names.5,clusters.5,by="OTU")
dim(names.5)

names.1 = data.frame(OTU=taxa_names(ps.1))
head(names.1)
names.1 = plyr::join(names.1,clusters.1,by="OTU")
dim(names.1)

# How many unique taxa in the phyloseq objects?
allnames = c(paste(names.1$cluster),paste(names.5$cluster))
head(allnames)
tail(allnames)
length(unique(allnames))
# 25870 unique OTU "clusters"
ps.merged.fasta
# Using the fastas alone, we had 25,788 unique OTUs.
# Fewer OTUs is "better", because it suggests the merging was more complete
# This suggests the alignment wasn't perfect or lost a few taxa (~100 i.e. 50)
# Out of 3,299 shared OTUs, missing 50 isn't too bad.
# Will need to check that the counts are the same - if we didn't lose any.

ps.5
sum(taxa_names(ps.5)==names.5$OTU)
ps.1
sum(taxa_names(ps.1)==names.1$OTU)
# Yes, both are in the right order.

# Running merge OTUs
ps.1.aln = ps.1
ps.5.aln = ps.5
# Get the taxonomy table and add the clustering data
taxtab.1 = data.frame(as(tax_table(ps.1.aln),"matrix"))
taxtab.1$Merge = names.1$cluster
taxtab.1$OTU = row.names(taxtab.1)
taxtab.1.temp = tax_table(taxtab.1)
taxa_names(taxtab.1.temp) = taxtab.1$OTU
tax_table(ps.1.aln)=tax_table(taxtab.1.temp)
colnames(tax_table(ps.1.aln)) = c("Domain","Phylum","Class","Order","Family",
                              "Genus","Species","Merge","OTU")

taxtab.5 = data.frame(as(tax_table(ps.5.aln),"matrix"))
taxtab.5$Merge = names.5$cluster
taxtab.5$OTU = row.names(taxtab.5)
taxtab.5.temp = tax_table(taxtab.5)
taxa_names(taxtab.5.temp)=taxtab.5$OTU
tax_table(ps.5.aln)=tax_table(taxtab.5.temp)
colnames(tax_table(ps.5.aln)) = c("Domain","Phylum","Class","Order","Family",
                              "Genus","Species","Merge","OTU")

# Now we have the merge column added in the OTU tables.
# Can work with these taxa as they were, using the merge reference as the join-up factor.
# Not a bad idea.
# If we do want to work with them together, do need to merge.
taxa_names(ps.5.aln) = paste(tax_table(ps.5.aln)[,"Merge"])
taxa_names(ps.1.aln) = paste(tax_table(ps.1.aln)[,"Merge"])


# So, we should be able to glom taxa based on that, after merging phyloseq objects
ps.merged.aln = merge_phyloseq(ps.5.aln,ps.1.aln)
ps.merged.aln

# How many total observations do we have in each dataset?
sum(sample_sums(ps.merged.aln))
sum(sample_sums(ps.merged.fasta))
# The same number. So neither merging approach dropped anything.
# Thus, we should likely take the dataset with the fewest OTUs - "best" merging
# Especially since this was the dataset that was just merged simply by identical sequences
ps.merged.fasta
ps.merged.aln

# Final merged phyloseq object
ps.merged = ps.merged.fasta
saveRDS(ps.merged,"ps.merged")
