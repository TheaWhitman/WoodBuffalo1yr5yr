# Load packages
library(phyloseq)
library(ggplot2)
library(dplyr)
library(reshape2)
library(vegan)
library(wesanderson)
library(metagenomeSeq)
library(biomformat)
library(reshape)
library(RColorBrewer)
library(ggforce)
library(Biostrings)

# Set working directory
setwd("~/Box/WhitmanLabMaster/WhitmanLab/Projects/WoodBuffalo/NWT-Fire-Microbes/WB2019/code")

# Load ps object
ps.merged = readRDS("ps.merged")
ps.merged

# Importing new pH values for 2019 soils
pH.2019 = read.csv("../data/Soils_data/Summaries/WB2019 pH means_04-30-21.csv")
# Make sure same sample IDs are present
as.factor(pH.2019$Sample_ID) %in% sample_data(ps.merged)$Sample_ID
sample_data(ps.merged)$Sample_ID %in% as.factor(pH.2019$Sample_ID)
# Yes, good, all sample IDs in our phyloseq object are found in the pH table
# we have a few measured pH values that aren't in our dataset - that's fine

# Check out the ones in the pH table that are not in the phyloseq object and vice versa
pH.2019[!(as.factor(pH.2019$Sample_ID) %in% sample_data(ps.merged)$Sample_ID),]

# 15S-NT-U06M wasn't in the 1-year dataset
  # maybe soil sampling didn't hit the mineral horizon then
# 15S-NT-29O isn't in the 5-year sequencing dataset - site 29
# Note: WB19-38-B was the new site established because the previous site had been logged between sampling years. (15S-NT-38)
  # So, for the purposes of this paper, 38 might still be better comparison.
as.factor(pH.2019$Sample_ID)
sample_data(ps.merged)$Sample_ID

#### Note: U06 / U07-O - U06-O is written on the soil bag used for pH, U07-O for DNA

# Updating pH data
add.pH = data.frame(sample_data(ps.merged)[(sample_data(ps.merged)$Sample_ID %in% pH.2019$Sample_ID) & (sample_data(ps.merged)$Years_Since_Fire==5),])
add.pH = merge(add.pH,pH.2019,by="Sample_ID")
colnames(add.pH)
row.names(add.pH) = add.pH$X.SampleID
add.pH$pH = add.pH$Mean_pH
add.pH = add.pH[,1:110]
remaining = data.frame(sample_data(ps.merged)[!((sample_data(ps.merged)$Sample_ID %in% pH.2019$Sample_ID) & (sample_data(ps.merged)$Years_Since_Fire==5)),])
dim(remaining)
dim(add.pH)
samdat = sample_data(rbind(remaining,add.pH))
sample_data(ps.merged)=samdat

# Similarly, update soil C values
# Importing new pH values for 2019 soils
C.2019 = read.csv("../data/Soils_data/Summaries/WB2019_MandOhorizon_TCandTN_Summary_06-16-2021.csv")

# Make sure same sample IDs are present
as.factor(C.2019$Sample_ID) %in% sample_data(ps.merged)$Sample_ID
sample_data(ps.merged)$Sample_ID %in% as.factor(C.2019$Sample_ID)
# As with pH, all sample IDs in our phyloseq object are found in the pH table

# Same samples in the CN table that are not in the phyloseq object
C.2019[!(as.factor(C.2019$Sample_ID) %in% sample_data(ps.merged)$Sample_ID),]

# Updating C data
add.C = data.frame(sample_data(ps.merged)[(sample_data(ps.merged)$Sample_ID %in% C.2019$Sample_ID) & (sample_data(ps.merged)$Years_Since_Fire==5),])
add.C = merge(add.C,C.2019,by="Sample_ID")
colnames(add.C)
row.names(add.C) = add.C$X.SampleID
add.C$TC_pct = add.C$TC_dry_pct
add.C$Total_N_pct = add.C$TN_dry_pct
# Keep just the original columns
add.C = add.C[,1:110]
remaining = data.frame(sample_data(ps.merged)[!((sample_data(ps.merged)$Sample_ID %in% C.2019$Sample_ID) & (sample_data(ps.merged)$Years_Since_Fire==5)),])
dim(remaining)
dim(add.C)
samdat = sample_data(rbind(remaining,add.C))
sample_data(ps.merged)=samdat

# Transform to relative abundance
ps.merged.norm = transform_sample_counts(ps.merged, function(x) x/sum(x))
ps.merged.hell = transform_sample_counts(ps.merged, function(x) (x/sum(x))^0.5)

# Create ordinations
#ord.hell.PCoA = ordinate(ps.merged.hell,method="PCoA",distance="bray")
#ord.norm.PCoA = ordinate(ps.merged.norm,method="PCoA",distance="bray")
ord.hell.nmds = ordinate(ps.merged.hell,method="NMDS",distance="bray",k=3,maxit=1000)
#ord.norm.nmds = ordinate(ps.merged.norm,method="NMDS",distance="bray")

# Plot ordinations
#ord = ord.hell.PCoA
#ord = ord.norm.PCoA
ord = ord.hell.nmds
#ord = ord.norm.nmds

ps = ps.merged.hell
#ps = ps.merged.norm



palette = brewer.pal(9, "YlGnBu")[3:9]

# Quick ordination
plot_ordination(ps,ord,color="pH",shape="Burned_Unburned",axes=c(1,2))
plot_ordination(ps,ord,color="pH",shape="Burned_Unburned",axes=c(1,2),label="Site_ID")+facet_wrap(~sample_data(ps)$Org_or_Min~sample_data(ps)$Years_Since_Fire)
plot_ordination(ps,ord,color="pH",shape="Burned_Unburned",axes=c(1,2))+facet_grid(~sample_data(ps)$Years_Since_Fire) + theme_bw() + scale_colour_gradientn(colours=rev(palette))

# I want to connect points that are the same site/sample
# For PCoA
#x = data.frame(ord$vectors)$Axis.1
#y = data.frame(ord$vectors)$Axis.2

# For NMDS
x = data.frame(ord$points)$MDS1
y = data.frame(ord$points)$MDS2

sample_data(ps)$x = x
sample_data(ps)$y = y

# These are the co-ordinates of each point.
# We need to know how to join them.

sample_data(ps)$Site_ID
sample_data(ps)$Org_or_Min
sample_data(ps)$Pairs = paste(sample_data(ps)$Site_ID,sample_data(ps)$Org_or_Min,sep="_")
sample_data(ps)$Pairs


# We will use geom_segment, which takes x, y, and xend, yend.
Segments = data.frame(sample_data(ps)) %>%
     dplyr::select(Pairs,Years_Since_Fire,Burned_Unburned,Org_or_Min,Burn_Severity_Index,x,y)%>%
    arrange(Pairs)
head(Segments)

Segments = melt(Segments,id=c("Pairs","Years_Since_Fire","Burned_Unburned","Org_or_Min","Burn_Severity_Index"))
head(Segments)

Segments = dcast(Segments,Pairs+Burned_Unburned+Org_or_Min+Burn_Severity_Index~Years_Since_Fire+variable,value.var = "value",fun=mean)
colnames(Segments) = c("Pair","Burned_Unburned","Org_or_Min","Burn_Severity_Index","One_x","One_y","Five_x","Five_y")
head(Segments)
sample_names(ps.merged.hell)

# Creating a plot with the two years of data, linked by lines.
d = data.frame(sample_data(ps))
d$BSI.plot = d$Burn_Severity_Index
d$BSI.plot[d$Years_Since_Fire=="1"]=1

# Set up labels for facets
Labs = c("Top Mineral Horizon", "Organic Horizon")
names(Labs) = c("M", "O")

p = ggplot(d) + theme_bw()
p = p + geom_segment(data=Segments,aes(x=One_x,y=One_y,xend=Five_x,yend=Five_y),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="lightgrey")
p = p + geom_point(data=d,aes(x=x,y=y,shape=Burned_Unburned,fill=Years_Since_Fire),
                   size=3, alpha=0.8)
p = p + scale_fill_manual(values=c("lightblue","darkblue"))
p = p + scale_shape_manual(values=c(21,24))
p = p + facet_wrap(~Org_or_Min,labeller = labeller(Org_or_Min = Labs))
p = p + xlab("Axis 1") + ylab("Axis 2")
p = p + guides(shape=guide_legend(title="Burned or\nUnburned"),fill=guide_legend(title="Years Since\nFire",override.aes = list(shape = 21)))
p

############# Plotting BSI as colour and facet
split1 = 3.5
split2 = 4
d$BSI.cat = as.factor(ifelse(d$Burn_Severity_Index == 1,"Unburned", ifelse((d$Burn_Severity_Index>1 & d$Burn_Severity_Index< split1),"Low", ifelse(d$Burn_Severity_Index>split1 & d$Burn_Severity_Index<split2,"Medium","High"))))
Segments$BSI.cat = as.factor(ifelse(Segments$Burn_Severity_Index == 1,"Unburned", ifelse((Segments$Burn_Severity_Index>1 & Segments$Burn_Severity_Index< split1),"Low", ifelse(Segments$Burn_Severity_Index>split1 & Segments$Burn_Severity_Index<split2,"Medium","High"))))
d$BSI.cat = ordered(d$BSI.cat,levels=c("Unburned","Low","Medium","High"))
Segments$BSI.cat = ordered(Segments$BSI.cat,levels=c("Unburned","Low","Medium","High"))
palette = wes_palette("Darjeeling1")

d.plot=d
Segments.plot=Segments

p = ggplot(d.plot) + theme_bw()
p = p + geom_segment(data=Segments.plot,aes(x=One_x,y=One_y,xend=Five_x,yend=Five_y),arrow = arrow(type="closed",length = unit(0.25,"cm")), color="lightgrey")
p = p + geom_point(data=d.plot,aes(x=x,y=y,shape=Burned_Unburned,fill=Years_Since_Fire),
                   size=3, alpha=0.9)
p = p + scale_shape_manual(values=c(21,24))
#p = p + scale_colour_gradient(low = "lightgrey",high="black")
#p = p + scale_colour_gradient2(low = palette[2],mid = palette[3],high=palette[1],midpoint=2)
#p = p + scale_colour_manual(values=palette[c(1,3,4,2)])
#p = p + scale_fill_manual(values=palette[c(1,3,4,2)])
p = p + scale_fill_manual(values=c("lightblue","darkblue"))
p = p + facet_grid(~BSI.cat)
p = p + xlab("Axis 1") + ylab("Axis 2")
p = p + guides(shape=guide_legend(title="Burned or\nUnburned"),color=guide_colorbar(title="Burn Severity\nIndex"))
p

############################################
# We're probably interested in seeing how the same site in year 1 differs from that site in year 5
ps = ps.merged.hell
sample_data(ps)$Pairs = paste(sample_data(ps)$Site_ID,sample_data(ps)$Org_or_Min,sep="_")

Dist.mb = as.matrix(vegdist(otu_table(ps), method="bray", type="samples"))
Dist.mb[upper.tri(Dist.mb)] = NA
df = data.frame(melt(Dist.mb))
head(df)

#Generates a dataframe with each contrast and the dissimilarity for the mb comm
colnames(df) = c("Sample_ID_1","Sample_ID_2","Mb_dist")
df = df[df$Mb_dist!=0,]
df = df[!is.na(df$Mb_dist),]
head(df)

# Making a matrix with one data entry for each site
SamDat = sample_data(ps)
SamDat$Sample_ID = row.names(SamDat)

# Need to add datasets for each site type, and then whether they are the same or not.
# Let's start with wetland vs. upland
for (i in df$Sample_ID_1){
  df$Pairs_1[df$Sample_ID_1==i] = paste(SamDat$Pairs[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$Pairs_2[df$Sample_ID_2==i] = paste(SamDat$Pairs[SamDat$Sample_ID==i])
}
df$Pairs = ifelse(df$Pairs_1==df$Pairs_2,df$Pairs_1,"Different")

# Selecting only comparisons between unburned and burned sites, same horizon, same veg comm
df = df[df$Pairs != "Different",]
dim(df)

# Want to add back in site info
siteinfo = data.frame(SamDat)%>%
  dplyr::group_by(Pairs,Burned_Unburned,Veg_Comm,Burn_Severity_Index,Org_or_Min)%>%
  summarize(N=n())
head(siteinfo)

# Add the data and clean it up
df = merge(df,siteinfo,by="Pairs")%>%
  select(Pairs,Burned_Unburned,Veg_Comm,Burn_Severity_Index,Org_or_Min,Mb_dist)
dim(df)
head(df)

# Check whether relationship is significant
# May want to test with burnedonly, too
df.burnedonly = df %>%
  filter(Burned_Unburned=="Burned")
# Kruskal Wallis test - one way anova by ranks
kruskal.test(Mb_dist~Burn_Severity_Index,data=df.burnedonly)

lm = lm(Mb_dist~Burn_Severity_Index,data=df)
summary(lm)
lm$coefficients[2]

# Non-parametric - independed Mann-Whitney U test
wilcox.test(Mb_dist~Burned_Unburned,data=df)
# Kruskal Wallis test - one way anova by ranks
kruskal.test(Mb_dist~Burn_Severity_Index,data=df)
kruskal.test(Mb_dist~Burned_Unburned,data=df)
# Both sig.

p = ggplot(df,aes(x=Burn_Severity_Index,y=Mb_dist,color=Veg_Comm, shape=Org_or_Min))
p = p + geom_point(size=2)
p = p + scale_color_manual(values=wes_palette("Darjeeling1")[c(1,4,3,2)])
p = p + scale_fill_manual(values=wes_palette("Darjeeling1")[c(1,4,3,2)])
p = p + scale_shape_manual(values=c(16,17,23,25,15))
linedata = data.frame(Land_Class=c("Upland","Wetland"), slope=c(lm$coefficients[2],NA), intercept=c(lm$coefficients[1],NA))
p = p + geom_abline(data=linedata, aes(slope=slope,intercept=intercept),linetype=2)
p = p + ylim(values=c(0,1))
p = p + theme_bw()
#p = p + facet_grid(~Veg_Comm~Org_or_Min)
p = p + ylab("Dissimilarity between 1 and 5 year-\ncommunities (Bray-Curtis)") + xlab("Burn Severity Index")
p = p + guides(color=guide_legend(title="Vegetation\nCommunity"),shape=guide_legend(title="Organic or\nMineral Horizon"))
p

# There is a significant positive relationship betwwen 1-5 year dissimilarity and burn severity
# So, we might think that the sites that had more severe burns changed more between sampling times
# However, when testing only the burned sites, this pattern did not hold.
# Thus it is likely more drive by the unburned sites remaining similar to the previous timepoint
# while all the burned sites shifted similarly.
# We can test this:

ps.1.hell = subset_samples(ps.merged.hell,Years_Since_Fire==1)
ps.5.hell = subset_samples(ps.merged.hell,Years_Since_Fire==5)

ps.1.hell = prune_taxa(taxa_sums(ps.1.hell)>0, ps.1.hell)
ps.5.hell = prune_taxa(taxa_sums(ps.5.hell)>0, ps.5.hell)

# Adjust this for running for 1 or 5 years post fire
physeq = ps.5.hell

Dist.mb = as.matrix(phyloseq::distance(physeq, method="bray", type="samples"))
Dist.mb[upper.tri(Dist.mb)] = NA
df = data.frame(melt(Dist.mb))
head(df)

#Generates a dataframe with each contrast and the dissimilarity for the mb comm
colnames(df) = c("Sample_ID_1","Sample_ID_2","Mb_dist")
df = df[df$Mb_dist!=0,]
df = df[!is.na(df$Mb_dist),]
head(df)

# Making a matrix with one data entry for each site
SamDat = sample_data(physeq)
SamDat$Sample_ID = row.names(SamDat)

# Need to add datasets for each site type, and then whether they are the same or not.
# Let's start with wetland vs. upland
for (i in df$Sample_ID_1){
  df$Land_Class_1[df$Sample_ID_1==i] = paste(SamDat$Land_Class[SamDat$Sample_ID==i])
}

for (i in df$Sample_ID_2){
  df$Land_Class_2[df$Sample_ID_2==i] = paste(SamDat$Land_Class[SamDat$Sample_ID==i])
}
df$Land_Class = ifelse(df$Land_Class_1==df$Land_Class_2,df$Land_Class_1,"Different")

# Let's add Veg_Comm
for (i in df$Sample_ID_1){
  df$Veg_Comm_1[df$Sample_ID_1==i] = paste(SamDat$Veg_Comm[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$Veg_Comm_2[df$Sample_ID_2==i] = paste(SamDat$Veg_Comm[SamDat$Sample_ID==i])
}
df$Veg_Comm = ifelse(df$Veg_Comm_1==df$Veg_Comm_2,df$Veg_Comm_1,"Different")

# Let's add Org_or_Min
for (i in df$Sample_ID_1){
  df$Org_or_Min_1[df$Sample_ID_1==i] = paste(SamDat$Org_or_Min[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$Org_or_Min_2[df$Sample_ID_2==i] = paste(SamDat$Org_or_Min[SamDat$Sample_ID==i])
}
df$Org_or_Min = ifelse(df$Org_or_Min_1==df$Org_or_Min_2,df$Org_or_Min_1,"Different")

# Let's add Burned_Unburned
for (i in df$Sample_ID_1){
  df$Burned_Unburned_1[df$Sample_ID_1==i] = paste(SamDat$Burned_Unburned[SamDat$Sample_ID==i])
}
for (i in df$Sample_ID_2){
  df$Burned_Unburned_2[df$Sample_ID_2==i] = paste(SamDat$Burned_Unburned[SamDat$Sample_ID==i])
}
df$Burned_Unburned = ifelse(df$Burned_Unburned_1==df$Burned_Unburned_2,df$Burned_Unburned_1,"Different")

# And let's add Burn_Severity_Index
for (i in df$Sample_ID_1){
  df$Burn_Severity_Index_1[df$Sample_ID_1==i] = SamDat$Burn_Severity_Index[SamDat$Sample_ID==i]
}
for (i in df$Sample_ID_2){
  df$Burn_Severity_Index_2[df$Sample_ID_2==i] = SamDat$Burn_Severity_Index[SamDat$Sample_ID==i]
}
# And calculate the differences between BSI for the two samples
df$Burn_Severity_Index_1 = as.numeric(df$Burn_Severity_Index_1)
df$Burn_Severity_Index_2 = as.numeric(df$Burn_Severity_Index_2)
df$Burn_Severity_Index_Diff = abs(as.numeric(df$Burn_Severity_Index_1) - as.numeric(df$Burn_Severity_Index_2))
df$Burn_Severity_Index_Mean = (as.numeric(df$Burn_Severity_Index_1) + as.numeric(df$Burn_Severity_Index_2))/2
df$Burn_Severity_Index_Ratio = as.numeric(df$Burn_Severity_Index_1) / as.numeric(df$Burn_Severity_Index_2)
df$Burn_Severity_Index_Class = ifelse(df$Burn_Severity_Index_1>mean(df$Burn_Severity_Index_1) &
                                        df$Burn_Severity_Index_2>mean(df$Burn_Severity_Index_1),"high",
                                      ifelse(df$Burn_Severity_Index_1<mean(df$Burn_Severity_Index_1) &
                                               df$Burn_Severity_Index_2<mean(df$Burn_Severity_Index_1),"low","mixed"))
df$Burn_Severity_Index_Class = ordered(df$Burn_Severity_Index_Class, levels = c("high", "mixed", "low"))

# Selecting only comparisons between unburned and burned sites, same horizon, same veg comm
df = df[(df$Burned_Unburned_1=="Unburned" | df$Burned_Unburned_2=="Unburned"),]
df = df[df$Org_or_Min != "Different",]
df = df[df$Veg_Comm != "Different",]
dim(df)

# Is dissimilarity to unburned significantly related to burn severity index? (yes)
lm = lm(Mb_dist~Burn_Severity_Index_Diff,data=df)
summary(lm)
lm$coefficients[2]

# We want to compare the 1 year and 5 year datasets (run above code for each)
#df.5=df
#df.1=df

# Check normality first
shapiro.test(df.5$Mb_dist)
shapiro.test(df.1$Mb_dist)

# df.1 not normally distributed
# So, use Mann-Whitney U test
wilcox.test(df.5$Mb_dist,df.1$Mb_dist)

# There is a significant difference between the two years.

# If we wanted to test only the burned sites, we can do that, too
df.5.burned = df.5[!(df.5$Burned_Unburned_1=="Unburned" & df.5$Burned_Unburned_2=="Unburned"),]
df.1.burned = df.1[!(df.1$Burned_Unburned_1=="Unburned" & df.1$Burned_Unburned_2=="Unburned"),]
# Check normality first
shapiro.test(df.5.burned$Mb_dist)
shapiro.test(df.1.burned$Mb_dist)

#Mann-Whitney U test
wilcox.test(df.5.burned$Mb_dist,df.1.burned$Mb_dist)
# Still significant.

##### Next, we want to see whether the mean predicted 16S copy number 
##### is higher in the two datasets, and how it corresponds to burn severity.

# Adding RDP Classifier data

DNA = Biostrings::DNAStringSet(merged$OTU)
names(DNA)=merged$OTU
Biostrings::writeXStringSet(DNA,"mergedResponders.fasta")

# We already ran rrnDB classifier RDP Classifier version 2.12 https://rrndb.umms.med.umich.edu/estimate/run_classifier
# Using 0.8 cutoff.
# Importing the results of that run

# Sequence counts
# This file has what they classified each read as
RDP = read.csv("../data/Seqs-processing/rrnDB/mergedResponders/mergedResponders.tsv",header=FALSE,sep=";")

# We want to extract the genus they assigned for each OTU.
# Create function to extract genus, if present
GenusGenerator <- function(taxon) {
  Genus = strsplit(gsub(".*family\\s*|genus*", "", taxon),split="\t")[[1]][2]
  return(Genus)
}

# Extract the genus
RDP$GenusRDP = sapply(RDP$V1,GenusGenerator)
head(RDP)

# Might as well pull out OTU ID to be certain
OTUGenerator <- function(taxon) {
  OTU = strsplit(paste(taxon),split="\t")[[1]][1]
  return(OTU)
}
RDP$OTU = sapply(RDP$V1,OTUGenerator)
head(RDP)

# Ok, we've got what we need.
# Trim it down
RDP = RDP[,c(2,3)]

# Can now pull data from RRNDB.
# Reading in the rrnDB v5.5 file
rrnDB = read.csv("../../../../../../../../../../../../Volumes/SeagateBackupPlusDrive/Databases/rrnDB-5.5/rrnDB-5.5_pantaxa_stats_RDP.tsv",sep="\t")
rrnDB = read.csv("../../../../../../../../../../Volumes/SeagateBackupPlusDrive/Databases/rrnDB-5.5/rrnDB-5.5_pantaxa_stats_RDP.tsv",sep="\t")
head(rrnDB)

# Creating a list of genera in the DB
rrnDBGenera = as.character(rrnDB[rrnDB$rank=="genus",]$name)

# Matching up genus name with predicted copy number
for (i in 1:length(RDP$GenusRDP)){
  GenusRDP = paste(RDP$GenusRDP[i])
  CopyNum = ifelse(GenusRDP %in% rrnDBGenera, rrnDB[rrnDB$name==GenusRDP,9],"")
  RDP$CopyNum[i] = CopyNum
}
head(RDP)

mdf = psmelt(ps.merged.norm)

# Add the rrnDB copy number data to the melted phyloseq object
mdf = plyr::join(mdf,RDP,by="OTU")
mdf$CopyNum = as.numeric(mdf$CopyNum)
mdf$Abundance = as.numeric(mdf$Abundance)

# From Nemergut et al. (2016) - OTU data were then normalized (standardized) for copy number 
# by dividing by copy number. For each sample, we calculated the community aggregated trait value 
# (weighted mean) by taking the product of the estimated operon copy number and the relative abundance 
# for each OTU, and summing this value across all OTUs in a sample. 

# So, first, we divide abundance by copy number
# Then, we re-calculate the relative abundnace, now adjusted for copy number
# The risk there, is, for any organisms without assigned copy numbers, they are excluded from this calculation.
# However, I think we have a pretty good fraction of the community with copy numbers
# To check:

d = mdf %>%
  dplyr::group_by(Sample)%>%
  dplyr::filter(is.na(CopyNum))%>%
  dplyr::summarize(NoCopyNum = sum(Abundance))
hist(d$NoCopyNum)
# Not too bad, but a wide range
# Could assume that unassigned taxa have mean copy number across dataset:
meanCopyNum = RDP%>%
  group_by(GenusRDP,CopyNum)%>%
  summarize(N=n())
meanCopyNum = mean(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)

# Calculating weighted mean copy numbers (with optional mean replacement; not used):
df = mdf %>%
  dplyr::mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",meanCopyNum,CopyNum))%>%
  dplyr::filter(!is.na(CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::group_by(Sample,Site_ID,Org_or_Min,Veg_Comm,Land_Class,pH,TC_pct,CBI,Understory_CBI,CFSI,Burn_Severity_Index,RBR,Burned_Unburned,Years_Since_Fire)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))
hist(df$WtMeanCopyNum)

### Now we want to know how the weighted mean copy number
# changed between the two years
# We need to rearrange the data table
df.sm = data.frame(df) %>% group_by() %>% dplyr::select(Site_ID,Org_or_Min,Veg_Comm,Land_Class,Burn_Severity_Index,Burned_Unburned,Years_Since_Fire,WtMeanCopyNum)
df.sm$Years_Since_Fire=as.factor(df.sm$Years_Since_Fire)
df.sm = dcast(data = df.sm,formula = Site_ID+Org_or_Min+Veg_Comm+Land_Class+Burn_Severity_Index+Burned_Unburned~Years_Since_Fire,value.var=c("WtMeanCopyNum"))
colnames(df.sm)[7:8]=c("WtMeanCopyNum.1","WtMeanCopyNum.5")

# Plot the results
p = ggplot(df.sm,aes(x=WtMeanCopyNum.1,y=WtMeanCopyNum.5,
                     color=Burn_Severity_Index-1,shape=Org_or_Min))
p = p + theme_bw()
p = p + geom_point(size=3,alpha=0.8)
p = p + geom_abline(slope = 1,intercept=0,linetype="dashed")
p = p + scale_colour_gradient(low = "lightgrey",high="black")
p = p + xlab("Weighted mean predicted\n16S rRNA gene copy number\none year post-fire")
p = p + ylab("Weighted mean predicted\n16S rRNA gene copy number\nfive years post-fire")
p = p + ylim(c(1,4.5))+xlim(c(1,4.5))
p = p + guides(shape=guide_legend(title="Organic or\nMineral\nHorizon"),color=guide_colourbar(title="Burn\nSeverity\nIndex"))
p

# Is this apparent decline significant?
# First, are they normally distributed?
shapiro.test(df.sm$WtMeanCopyNum.1)
shapiro.test(df.sm$WtMeanCopyNum.5)
wilcox.test(df.sm$WtMeanCopyNum.1, df.sm$WtMeanCopyNum.5, paired = TRUE)

kruskal.test(WtMeanCopyNum.5~Burn_Severity_Index,data=df.sm)
summary(aov(WtMeanCopyNum.5~Burn_Severity_Index,data=df.sm))

############## Our next questions will be whether we have the same responders year-to-year, and 
### how their burn response differs.

# We will basically just analyze them separately, and then join by OTUs.

ps.1 = prune_samples(sample_data(ps.merged)$Years_Since_Fire=="1", ps.merged)
ps.5 = prune_samples(sample_data(ps.merged)$Years_Since_Fire=="5", ps.merged)

# Remove no-count taxa
ps.1 = prune_taxa(taxa_sums(ps.1)>0,ps.1)
ps.5 = prune_taxa(taxa_sums(ps.5)>0,ps.5)

# Set phyloseq to non-normalized one, with only non-zero taxa
ps.biom = ps.1

# Copy phyloseq object (note, not using normalized ps)
otutab = data.frame(as.matrix(t(otu_table(ps.biom))))
taxtab = data.frame(as(tax_table(ps.biom),"matrix"))
samdat = data.frame(sample_data(ps.biom))

biom = make_biom(data = otutab, observation_metadata=taxtab, sample_metadata=samdat)
biom
# turn phyloseq into a biom table with OTU table, taxonomy, and sample data

biom.MRexp = biom2MRexperiment(biom)
# turns our biom file into the type of file needed for this analysis (an MRexperiment object)
biom.MRexp
# Create .MRexp object

# Remove na for key factors
## Only need this for ps.1 ...
#biom.MRexp = biom.MRexp[,-which(is.na(pData(biom.MRexp)$pH))]

#Preparing object for metegenomeseq analysis
MRexp = biom.MRexp
class(MRexp)
#Stat = cumNormStatFast(MRexp, qFlag=TRUE, pFlag = FALSE)
MRexp = cumNorm(MRexp, cumNormStat(MRexp, qFlag=TRUE, pFlag = FALSE))
ModelData = pData(MRexp)

dim(ModelData)
dim(MRexp)

# Setting factors to be proper values
ModelData$pH = as.numeric(ModelData$pH)
ModelData$TC_pct = as.numeric(ModelData$TC_pct)
ModelData$Veg_Comm = make.names(ModelData$Veg_Comm)
ModelData$Burned_Unburned = as.factor(ModelData$Burned_Unburned)
ModelData$Burned_Unburned = ordered(ModelData$Burned_Unburned, levels = c('Unburned','Burned'))
# Establishing model formula
model = model.matrix(~Veg_Comm+TC_pct+pH+Burned_Unburned, data=ModelData)

# Assigning settings
?zigControl
settings = zigControl(tol = 1e-04, maxit = 30, verbose = TRUE, dfMethod = "default", pvalMethod = "default")

# Running fit
fit = fitZig(obj=MRexp, mod=model, control = settings, useCSSoffset = TRUE, zeroMod = NULL, useMixedModel = FALSE)

# Effective sample size is calculated, and the average values of this is determined
EffSamp = calculateEffectiveSamples(fit)
MeanEffSamp = mean(EffSamp[!is.na(EffSamp)])

# Find rare features
rareFeatures = which(rowSums(MRcounts(MRexp) > 0) < MeanEffSamp)
# These are the taxa that that had less than the average number of effective samples
# As recommended in the vignette: https://www.bioconductor.org/packages/devel/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
MRexp = MRexp[-rareFeatures, ]
# Take the data object and remove the rareFeatures (taxa)
# Re-run the analyses below.
MRexp = cumNorm(MRexp, cumNormStat(MRexp, qFlag=TRUE, pFlag = FALSE))
ModelData = pData(MRexp)

# Resetting the model parameters
ModelData$pH = as.numeric(ModelData$pH)
ModelData$TC_pct = as.numeric(ModelData$TC_pct)
ModelData$Veg_Comm = make.names(ModelData$Veg_Comm)
ModelData$Burned_Unburned = as.factor(ModelData$Burned_Unburned)
ModelData$Burned_Unburned = ordered(ModelData$Burned_Unburned, levels = c('Unburned','Burned'))

model = model.matrix(~Veg_Comm+pH+TC_pct+Burned_Unburned, data=ModelData)
# Creating the model that we will use to analyze our data.
### Need to check on total C results - O horizon sample data zeros
fit = fitZig(obj=MRexp, mod=model, control = settings, useCSSoffset = TRUE, zeroMod = NULL, useMixedModel = FALSE)
#Re-running the fit with the EffSamp - normalized data

modeldesign = fit@fit$design
modelfit = fit@fit
modelfit.treat = treat(modelfit, lfc=0)
dim(modelfit.treat)
resultsBurn=topTreat(modelfit.treat, coef=7, number=dim(modelfit.treat)[1])
results=merge(data.frame(),resultsBurn,by=0, all=TRUE, suffixes=c(".Burn"))
# Extracting the results of interest

row.names(results)=results$Row.names
# Carrying the proper OTU IDs along
results= merge(results,fData(biom.MRexp),by=0,all=TRUE)

results = results[,2:dim(results)[2]]
colnames(results)[1]="OTU"
results$sigBurn = ifelse(results$adj.P.Val<0.05,1,0.5)
# Creating a variable that sets what will be the alpha values for plotting the results

# Summarizing and saving results
#Results.1 = results
#Results.5 = results

#saveRDS(Results.5,"Results.5.17.01.2022")
#saveRDS(Results.5,"Results.5")
#saveRDS(Results.1,"Results.1.17.01.2022")
#saveRDS(Results.1,"Results.1")
Results.1 = readRDS("Results.1.17.01.2022")
Results.5 = readRDS("Results.5.17.01.2022")

# Now we want to combine the two dataframes by OTU
# Each row should be an OTU, and the columns in keep are duplicated.
# Filling NAs for responses where not present
dim(Results.1)
dim(Results.5)

# There are so many taxa, it's a bit hard to merge
# But it might be possible

# We can get rid of the taxa that weren't even tested. That'll be helpful.
Results.1 = Results.1 %>%
    filter(!is.na(OTU))
dim(Results.1)
Results.5 = Results.5 %>%
  filter(!is.na(OTU))
dim(Results.5)

merged = merge(Results.1,Results.5,by="OTU", all=TRUE, suffixes=c(".OneYear",".FiveYears"))
head(merged)
colnames(merged)
dim(merged)

# Make a significance column
merged$Sig = ifelse(merged$sigBurn.OneYear==1 & merged$sigBurn.FiveYears==1,"Both years",
                    ifelse(merged$sigBurn.OneYear==1 & merged$sigBurn.FiveYears!=1,"One year",
                           ifelse(merged$sigBurn.OneYear!=1 & merged$sigBurn.FiveYears==1,"Five years","Not sig")))
merged$MeanAveExpr = (merged$AveExpr.OneYear+merged$AveExpr.FiveYears)/2

# Now, to plot it!
AbundPhyla = c("Acidobacteria","Actinobacteria","Proteobacteria","Chloroflexi","Planctomycetes","Verrucomicrobia")
merged.plot = merged%>%
    filter(Sig != "Not sig")%>%
    filter(Phylum.OneYear %in% AbundPhyla)

# Identifying specific taxa of interest
# Year-one responsive Arthrobacters
sq1 = "AAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTCGCGTCTGCTGTGAAAGACCGGGGCTCAACTCCGGTTCTGCAGTGGGTACGGGCAGACTAGAGTGCAGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCTGTAACTGACGCTGAGGAGCGAAAGCATGG"
sq7 = "AAGCGTTATCCGGAATTATTGGGCGTAAAGAGCTCGTAGGCGGTTTGTCGCGTCTGCTGTGAAAGCCCGGGGCTCAACCCCGGGTCTGCAGTGGGTACGGGCAGACTAGAGTGCAGTAGGGGAGACTGGAATTCCTGGTGTAGCGGTGAAATGCGCAGATATCAGGAGGAACACCGATGGCGAAGGCAGGTCTCTGGGCTGTAACTGACGCTGAGGAGCGAAAGCATGG"

Blasto = merged.plot %>%
  filter(Genus.OneYear=="Blastococcus")%>%
  arrange(-MeanAveExpr)

Blasto

Massilia = merged.plot %>%
  filter(Genus.OneYear=="Massilia")%>%
  arrange(-MeanAveExpr)

Massilia

Micrococc = merged.plot %>%
  filter(Family.OneYear=="Micrococcaceae")%>%
  arrange(-MeanAveExpr)

TaxaOfInterest = c(Micrococc$OTU,Blasto$OTU,sq1,sq7,Massilia$OTU[1])
#TaxaOfInterest = c(Blasto$OTU,sq1,sq7,Massilia$OTU[1])
#TaxaOfInterest = c(sq1,sq7)
#TaxaOfInterest = c(Massilia$OTU[2])
#TaxaOfInterest = c(Massilia$OTU)
#TaxaOfInterest = c(Blasto$OTU)

merged.plot$Interest = 0
merged.plot$Interest[merged.plot$OTU %in% TaxaOfInterest]=1
merged.plot[merged.plot$Interest==1,]

p = ggplot(merged.plot) + theme_bw()
p = p + geom_hline(yintercept=0,color="black")
p = p + geom_vline(xintercept=0,color="black")
p = p + facet_wrap(~Phylum.OneYear)
p = p + geom_point(aes(x=logFC.OneYear,y=logFC.FiveYears,fill=Sig,size=AveExpr.FiveYears),
                   shape=21,alpha=0.5)
p = p + geom_point(data=merged.plot[merged.plot$Interest==1,],aes(x=logFC.OneYear,y=logFC.FiveYears,size=AveExpr.FiveYears),shape=21,alpha=1,stroke=2)
p = p + scale_fill_manual(values=c("gold","darkblue","lightblue"))
p = p + geom_abline(slope=1,intercept=0,linetype="dashed")
p = p + ylab(expression(Five~year~response~to~burn~(log[2]-fold~change)))
p = p + xlab(expression(One~year~response~to~burn~(log[2]-fold~change)))
p = p + scale_size(guide="none")
p = p + guides(fill=guide_legend(title=element_text("Significant\nResponse")))
p = p + theme(strip.text.x = element_text(size = 12, face = "italic"))
p

# How many responders in each? (Replace Results.5 or Results.1)
Responders16S = Results.1 %>%
  filter(!is.na(OTU))%>%
  filter(adj.P.Val<0.05)%>%
  arrange(-logFC)%>%
  mutate(PosNeg = ifelse(logFC>0,"Pos",
                         ifelse(logFC<0,"Neg","Other")))%>%
  group_by(PosNeg)%>%
  summarize(Total=n())
head(Responders16S)  

### If we want to pare down the responders list
cutoff = 2
merged.cat = merged.plot %>%
      dplyr::mutate(Response_Category = ifelse(logFC.FiveYears>cutoff & logFC.OneYear>cutoff, "Both pos",
                  ifelse(logFC.FiveYears>cutoff, "Five pos",
                         ifelse(logFC.OneYear>cutoff,"One pos",
                                ifelse(logFC.FiveYears<(-cutoff) & logFC.OneYear<(-cutoff), "Both neg",
                                       ifelse(logFC.FiveYears<(-cutoff),"Five neg",
                                              ifelse(logFC.OneYear<(-cutoff),"One neg","None")))))))%>%
  filter(Response_Category != "None")

# This gives us the strong responders. We might want to filter by a cutoff for max abundance?
# I.e., must be at least this much
OTU.maxes = data.frame(OTU=taxa_names(ps.merged.norm),
                       MaxAbund = apply(otu_table(ps.merged.norm),2,function(x) max(x)),
                       MeanAbund = apply(otu_table(ps.merged.norm),2,function(x) mean(x)))
merged.cat = plyr::join(merged.cat,OTU.maxes,by="OTU",type="left")

# Ok, got them in, now another filter
cutoff = 0.02
merged.cat = merged.cat %>%
  filter(MaxAbund>cutoff)%>%
  arrange(Response_Category)
merged.cat

#write.csv(merged.cat,"../data/StrongAbundantResponders.17.01.2022.csv")

# It would be cool to look at estimated 16S copy number for the one, both, and five-year responders
# I would predict the one would be highest, both the next highest, then five years.
# However, if both is highest, it might suggest that really fast growth gives taxa an advantage
# that sets them up to dominate over time.

# We already imported the rrnDB info above.
# Get all our responsive taxa
merged.rdp = merged%>%
  filter(Sig != "Not sig")

# Add the rrnDB data to this merged file
merged.rdp = merge(merged.rdp,RDP,by="OTU")
merged.rdp$CopyNum
merged.rdp$CopyNum = as.numeric(paste(merged.rdp$CopyNum))

# Getting taxa that were positive responders in either year
merged.rdp.pos = merged.rdp%>%
    filter(!is.na(CopyNum))%>%
    filter(!is.na(Sig))%>%
    filter(Sig!="Not sig") %>%
    filter(logFC.FiveYears > 0 | logFC.OneYear>0)

# Test our hypothesis
merged.rdp.pos$Sig = ordered(merged.rdp.pos$Sig,levels=c("Both years","One year","Five years"))
p = ggplot(merged.rdp.pos) + theme_bw()
p = p + geom_boxplot(aes(x=Sig,y=CopyNum,fill=Sig),alpha=0.8)
p = p + scale_fill_manual(values=c("gold","lightblue","darkblue"),guide=FALSE)
p = p + ylab("Predicted 16S rRNA gene copy number")
p = p + xlab("Significant Enrichment in Burns")
p

# Copy number clearly not normally distributed
hist(merged.rdp.pos$CopyNum)
shapiro.test(merged.rdp.pos$CopyNum)
# There is a significant difference - not normal.
# Test using kruksal
kruskal.test(CopyNum~Sig,data=merged.rdp.pos)

# Is there a sig diff between both and one or five years?
kruskal.test(CopyNum~Sig,data=merged.rdp.pos[merged.rdp.pos$Sig != "Five years",])
kruskal.test(CopyNum~Sig,data=merged.rdp.pos[merged.rdp.pos$Sig != "Both years",])
kruskal.test(CopyNum~Sig,data=merged.rdp.pos[merged.rdp.pos$Sig != "One year",])
# Nope.
# Similarly:
wilcox.test(merged.rdp.pos[merged.rdp.pos$Sig == "Both years",]$CopyNum,merged.rdp.pos[merged.rdp.pos$Sig == "Five years",]$CopyNum)
wilcox.test(merged.rdp.pos[merged.rdp.pos$Sig == "One year",]$CopyNum,merged.rdp.pos[merged.rdp.pos$Sig == "Five years",]$CopyNum)
wilcox.test(merged.rdp.pos[merged.rdp.pos$Sig == "Both years",]$CopyNum,merged.rdp.pos[merged.rdp.pos$Sig == "One year",]$CopyNum)

### We might want to test just one year and just five year, without pulling out the "both" taxa

merged.rdp.5years = merged.rdp%>%
  filter(!is.na(CopyNum))%>%
  filter(!is.na(Sig))%>%
  filter(adj.P.Val.FiveYears < 0.05) %>%
  filter(logFC.FiveYears > 0)

merged.rdp.1year = merged.rdp%>%
  filter(!is.na(CopyNum))%>%
  filter(!is.na(Sig))%>%
  filter(adj.P.Val.OneYear < 0.05) %>%
  filter(logFC.OneYear > 0)

# Testing
wilcox.test(merged.rdp.1year$CopyNum,merged.rdp.5years$CopyNum)
## No, if we include all the taxa that are in one year, and all those in five years,
# including overlapping taxa, there isn't a significant difference.

# What about taxa that were negative responders in either year?
merged.rdp.neg = merged.rdp%>%
  filter(!is.na(CopyNum))%>%
  filter(!is.na(Sig))%>%
  filter(Sig!="Not sig") %>%
  filter(logFC.FiveYears < 0 | logFC.OneYear<0)

# Test our hypothesis
merged.rdp.neg$Sig = ordered(merged.rdp.neg$Sig,levels=c("Both years","One year","Five years"))
p = ggplot(merged.rdp.neg) + theme_bw()
p = p + geom_boxplot(aes(x=Sig,y=CopyNum,fill=Sig),alpha=0.8)
p = p + scale_fill_manual(values=c("gold","lightblue","darkblue"),guide=FALSE)
p = p + ylab("Predicted 16S rRNA gene copy number")
p = p + xlab("Significant Depletion in Burned Sites")
p

# Copy number clearly not normally distributed
hist(merged.rdp.neg$CopyNum)
shapiro.test(merged.rdp.neg$CopyNum)
# Test using kruksal
kruskal.test(CopyNum~Sig,data=merged.rdp.neg)
# There is not a significant difference for negative responders (kinda makes sense)

# Just another way to see that the positive-responding 
# five years-only taxa do not tend to have lower copy numbers
p = ggplot(merged.rdp.pos) + theme_bw()
p = p + geom_density(aes(x=CopyNum,fill=Sig),alpha=0.5)
p = p + scale_fill_manual(values=wes_palette("Zissou1")[c(5,4,2,3)])
p

# What about degree of response?
# Is the mean l2FC lower 5 years than 1 year after fire?
# Considering only taxa with consistent (same direction) response
merged.pos = merged%>%
    filter(logFC.FiveYears>0 & logFC.OneYear>0)

merged.neg = merged%>%
  filter(logFC.FiveYears<0 & logFC.OneYear<0)

# All non-normal
shapiro.test(merged.pos$logFC.OneYear)
shapiro.test(merged.pos$logFC.FiveYears)
shapiro.test(merged.neg$logFC.OneYear)
shapiro.test(merged.neg$logFC.FiveYears)

# Testing with Mann-Whitney
wilcox.test(merged.pos$logFC.OneYear,merged.pos$logFC.FiveYears)
wilcox.test(merged.neg$logFC.OneYear,merged.neg$logFC.FiveYears)
# For taxa that had significant responses in the same direction both years,
# There is a significant difference, (p=0.0004 pos, p < 0.0001 neg)

# Plotting these general responses
Results.5$Year = "5"
Results.1$Year = "1"
merged.stack = rbind(Results.5,Results.1)
merged.stack$PosNeg = ifelse(merged.stack$logFC>0,"Positive",
                             ifelse(merged.stack$logFC<0,"Negative",NA))

p = ggplot(merged.stack)+theme_bw()
p = p + geom_boxplot(aes(x=PosNeg,y=logFC,fill=Year),alpha=0.8)
p = p + scale_fill_manual(values=c("lightblue","darkblue"))
p = p + guides(fill=guide_legend(title=element_text("Years\nSince\nFire")))
p = p + ylab(expression(Response~to~burn~(log[2]-fold~change)))
p = p + xlab("Direction of Burn Response")
p

# Last, let's test Arthrobacter, Massilia, and Blastococcus
head(merged)
dim(merged)
taxaofinterest = c("Arthrobacter","Massilia","Blastococcus")
familyofinterest = c("Micrococcaceae")
merged.taxa = merged%>%
    filter(Genus.FiveYears %in% taxaofinterest | Genus.OneYear %in% taxaofinterest
           | Family.FiveYears %in% familyofinterest | Family.OneYear %in% familyofinterest)

dim(merged.taxa)
p = ggplot(merged.taxa) + theme_bw()
p = p + geom_point(aes(x=logFC.OneYear,y=logFC.FiveYears,fill=Genus.FiveYears, size=AveExpr.FiveYears),
                   shape=21,alpha=0.5)
p = p + geom_abline(slope=1,intercept=0,linetype="dashed")
p










############ Looking at the 5-year responders more ##################
head(Results.5)

colnames(Results.5)
Results.5$sigBurn = as.factor(Results.5$sigBurn)

p = ggplot(Results.5,aes(x=Phylum, y=logFC, alpha=sigBurn, fill=Phylum))
p = p + geom_jitter(shape=21, stroke=1, aes(size=AveExpr), width=0.2)
p = p + scale_size(guide=FALSE)
p = p + scale_alpha_manual(guide=FALSE, values=c(0.2,0.8))
#p = p + scale_fill_manual(values=c(brewer.pal(11,"Spectral"),"white"))
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(angle=45, hjust=1, size=14))
p = p + theme(axis.title.y = element_text(size=14))
p = p + theme(legend.position="none")
p = p + geom_hline(yintercept=0) 
p = p + ylab(expression(Response~to~burn~(log[2]-fold~change)))
p = p + theme(axis.text.x = element_text(face="italic"),
              legend.text = element_text(face="italic"),
              axis.title.x = element_blank())
p

## Actinobacteria
Results.5.Actinos = Results.5%>%
    filter(Phylum == "Actinobacteria")

p = ggplot(Results.5.Actinos,aes(x=Order, y=logFC, alpha=sigBurn, fill=Order))
p = p + geom_jitter(shape=21, stroke=1, aes(size=AveExpr), width=0.2)
p = p + scale_size(guide=FALSE)
p = p + scale_alpha_manual(guide=FALSE, values=c(0.2,0.8))
p = p + scale_fill_manual(values=c(brewer.pal(11,"Spectral"),"white","black","grey","darkgrey"))
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(angle=45, hjust=1, size=14))
p = p + theme(axis.title.y = element_text(size=14))
p = p + theme(legend.position="none")
p = p + geom_hline(yintercept=0) 
p = p + ylab(expression(Response~to~burn~(log[2]-fold~change)))
p = p + theme(axis.text.x = element_text(face="italic"),
              legend.text = element_text(face="italic"),
              axis.title.x = element_blank())
p

## Proteobacteria
Results.5.Actinos = Results.5%>%
  #filter(Phylum == "Proteobacteria")
  filter(Order %in% c("Rhodospirillales", "Rhizobiales"))

p = ggplot(Results.5.Actinos,aes(x=Family, y=logFC, alpha=sigBurn, fill=Family))
p = p + geom_jitter(shape=21, stroke=1, aes(size=AveExpr), width=0.2)
p = p + scale_size(guide=FALSE)
p = p + scale_alpha_manual(guide=FALSE, values=c(0.2,0.8))
#p = p + scale_fill_manual(values=c(brewer.pal(11,"Spectral"),"white","black","grey","darkgrey"))
p = p + theme_bw()
p = p + theme(axis.text.x = element_text(angle=45, hjust=1, size=14))
p = p + theme(axis.title.y = element_text(size=14))
p = p + theme(legend.position="none")
p = p + geom_hline(yintercept=0) 
p = p + ylab(expression(Response~to~burn~(log[2]-fold~change)))
p = p + theme(axis.text.x = element_text(face="italic"),
              legend.text = element_text(face="italic"),
              axis.title.x = element_blank())
p










################## Testing significant predictors of community composition #####################
# Creating the function to test different models of microbial community composition prediction

ps.1.norm = transform_sample_counts(ps.1, function(x) x/sum(x))
ps.1.hell = transform_sample_counts(ps.1, function(x) (x/sum(x))^0.5)
ps.5.norm = transform_sample_counts(ps.5, function(x) x/sum(x))
ps.5.hell = transform_sample_counts(ps.5, function(x) (x/sum(x))^0.5)


AdonisFunction = function(physeq, method="bray", Org_or_Min=c("O","M"), Land_Class=c("Upland","Wetland")){
  physeq = prune_samples(sample_data(physeq)$Land_Class %in% Land_Class, physeq)
  physeq = prune_samples(!is.na(sample_data(physeq)$pH), physeq)
  physeq = prune_samples(!is.na(sample_data(physeq)$TC_pct), physeq)
  physeq = prune_samples(!is.na(sample_data(physeq)$Sand_pct), physeq)
  print(physeq)
  df = as(sample_data(physeq), "data.frame")
  d = phyloseq::distance(physeq, method = method, weighted=TRUE)  
  d.adonis = adonis(d ~ 
                      + sample_data(physeq)$Veg_Comm
                    + as.numeric(sample_data(physeq)$Moisture_Regime)
                    + sample_data(physeq)$pH
                    + sample_data(physeq)$TC_pct
                    + sample_data(physeq)$Sand_pct
                    + sample_data(physeq)$Burned_Unburned                      
                    #+ sample_data(physeq)$RBR
                    #+ sample_data(physeq)$CFSI
                    #+ sample_data(physeq)$CBI
                    #+ sample_data(physeq)$Understory_CBI
                    #+ sample_data(physeq)$Overstory_CBI
                    #+ sample_data(physeq)$Burn_Severity_Index
                    #+ sample_data(physeq)$Pct_Exposed_Mineral
                    #+ sample_data(physeq)$Mean_Duff_Depth_cm
                    , df)
  d.adonis
}

AdonisFunction(physeq=ps.1.norm)
AdonisFunction(physeq=ps.5.norm)








############################# If we were gonna propose taxon response categories

# Decide required level of abundance
# Work from responders - three categories
# Might be better to place them in 2-D space (actually, that's the L2FC figure)




################ Copynumber figures ################


# We already ran rrnDB classifier RDP Classifier version 2.12 https://rrndb.umms.med.umich.edu/estimate/run_classifier
# Using 0.8 cutoff.
# Importing the results of that run

# Sequence counts
# This file has what they classified each read as
RDP = read.csv("../data/Seqs-processing/rrnDB/mergedResponders/mergedResponders.tsv",header=FALSE,sep=";")

# We want to extract the genus they assigned for each OTU.
# Create function to extract genus, if present
GenusGenerator <- function(taxon) {
  Genus = strsplit(gsub(".*family\\s*|genus*", "", taxon),split="\t")[[1]][2]
  return(Genus)
}


# Extract the genus
RDP$GenusRDP = sapply(RDP$V1,GenusGenerator)

head(RDP)

# Might as well pull out OTU ID to be certain
OTUGenerator <- function(taxon) {
  OTU = strsplit(paste(taxon),split="\t")[[1]][1]
  return(OTU)
}
RDP$OTU = sapply(RDP$V1,OTUGenerator)
head(RDP)

# Ok, we've got what we need.
# Trim it down
RDP = RDP[,c(2,3)]

# Can now pull data from RRNDB.
# Reading in the rrnDB v5.5 file
rrnDB = read.csv("../../../../../../../../../../Volumes/SeagateBackupPlusDrive/Databases/rrnDB-5.5/rrnDB-5.5_pantaxa_stats_RDP.tsv",sep="\t")
head(rrnDB)

# Creating a list of genera in the DB
rrnDBGenera = as.character(rrnDB[rrnDB$rank=="genus",]$name)

# Matching up genus name with predicted copy number
for (i in 1:length(RDP$GenusRDP)){
  GenusRDP = paste(RDP$GenusRDP[i])
  CopyNum = ifelse(GenusRDP %in% rrnDBGenera, rrnDB[rrnDB$name==GenusRDP,9],"")
  RDP$CopyNum[i] = CopyNum
}
head(RDP)

#  Melting phyloseq
mdf.merged.norm = psmelt(ps.merged.norm)

mdf.merged.norm = plyr::join(mdf.merged.norm,RDP,by="OTU")
mdf.merged.norm$CopyNum = as.numeric(mdf.merged.norm$CopyNum)
mdf.merged.norm$Abundance = as.numeric(mdf.merged.norm$Abundance)

meanCopyNum = RDP%>%
  group_by(GenusRDP,CopyNum)%>%
  summarize(N=n())
meanCopyNum = mean(as.numeric(meanCopyNum$CopyNum),na.rm=TRUE)

# Calculating weighted mean copy numbers:
df = mdf.merged.norm %>%
  dplyr::mutate(CopyNum = ifelse(is.na(CopyNum) | CopyNum =="NA",meanCopyNum,CopyNum))%>%
  dplyr::mutate(WtAbund = Abundance/CopyNum)%>%
  dplyr::group_by(Sample)%>%
  dplyr::mutate(AdjAbund = WtAbund/sum(WtAbund))%>%
  dplyr::mutate(WtCopyNum = CopyNum*AdjAbund)%>%
  dplyr::group_by(Years_Since_Fire,Site_ID,Sample,Org_or_Min,Veg_Comm,Land_Class,pH,TC_pct,CBI,Understory_CBI,CFSI,Burn_Severity_Index,RBR,Burned_Unburned)%>%
  dplyr::summarize(WtMeanCopyNum = sum(WtCopyNum,na.rm=TRUE))


# Plotting
p = ggplot(df,aes(x=Burn_Severity_Index-1,y=WtMeanCopyNum,fill=Veg_Comm,color=Years_Since_Fire, shape=Org_or_Min))
p = p + geom_point(alpha=0.8,size=2)
p = p + ylab("Weighted mean predicted 16S\nrRNA gene copy number")
p = p + xlab("Burn Severity Index")
p = p + theme_bw()
p = p + scale_fill_manual(values=wes_palette("Darjeeling1")[c(1,4,3,5,2)])
p = p + scale_color_manual(values=wes_palette("Darjeeling1")[c(1,4,3,5,2)])
p = p + scale_shape_manual(values=c(16,17,23,25,15))
p = p + guides(fill=guide_legend(title=element_blank()),color=guide_legend(title=element_blank()),shape=guide_legend(title=element_blank()))
p
