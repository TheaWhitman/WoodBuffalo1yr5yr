### This is the R code associated with figures and analysis from the paper
### Resilience not yet apparent in soil fungal communities of the 
### boreal forest from one to five years after wildfire across a severity gradient###

# Paper authors: Thea Whitman, Jamie Woolet, Miranda C. Sikora, Dana B. Johnson, Denyse A. Dawe, and Ellen Whitman 

# © 2025. This work is openly licensed via CC BY NC SA 4.0.
# https://creativecommons.org/licenses/by-nc-sa/4.0/

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
library(corrplot)
library(pals)
library(TeachingDemos)

# Import ps object
ps.merged = readRDS("ps.merged.aligned.0pct.glommed.spp.wNAs.noITSx.cutadapt")

# Add FunGuild assignments
FG = read.csv("../Seq-processing/TaxTabFunGuild.taxa.guilds.txt",sep="\t")
FG$trophicModeHighProb = ifelse(FG$confidenceRanking=="Highly Probable",FG$trophicMode,"")
FG$guildHighProb = ifelse(FG$confidenceRanking=="Highly Probable",FG$guild,"")
FG2 = tax_table(FG)
FG2 = FG2[,2:dim(FG2)[2]]
taxa_names(FG2)=FG[,1]

# Had wrong taxonomic categories for some reason
colnames(FG)[4:5]=c("Class","Order")
colnames(FG)
colnames(FG2)=colnames(FG[2:dim(FG)[2]])
head(FG2)
tax_table(ps.merged) = FG2

# Re-extract 1 and 5 year datasets
ps.1 = subset_samples(ps.merged,Years_Since_Fire==1)
ps.5 = subset_samples(ps.merged,Years_Since_Fire==5)
ps.1 = prune_taxa(taxa_sums(ps.1)>0, ps.1)
ps.5 = prune_taxa(taxa_sums(ps.5)>0, ps.5)

mean(sample_sums(ps.1))
mean(sample_sums(ps.5))

# Get shared taxa
shared = taxa_names(ps.1)[taxa_names(ps.1) %in% taxa_names(ps.5)]
ps.shared = prune_taxa(taxa_names(ps.merged) %in% shared, ps.merged)

# Transform to relative abundance
ps.norm = transform_sample_counts(ps.merged, function(x) x/sum(x))
ps.hell = transform_sample_counts(ps.merged, function(x) (x/sum(x))^0.5)

# Create melted ps object
mdf = psmelt(ps.norm)

### TABLE 1 ###

################## Testing significant predictors of community composition #####################
# Creating the function to test different models of microbial community composition prediction

ps.1.norm = transform_sample_counts(ps.1, function(x) x/sum(x))
ps.1.hell = transform_sample_counts(ps.1, function(x) (x/sum(x))^0.5)
ps.5.norm = transform_sample_counts(ps.5, function(x) x/sum(x))
ps.5.hell = transform_sample_counts(ps.5, function(x) (x/sum(x))^0.5)

# Set this according to 1 or 5
physeq = ps.1.hell
physeq = ps.5.hell

physeq = prune_samples(!is.na(sample_data(physeq)$pH), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$TC_pct), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$Sand_pct), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$Moisture_Regime), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$Burned_Unburned), physeq)  

df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray", weighted=TRUE)
d.adonis = adonis2(d ~ 
                     + sample_data(physeq)$Veg_Comm
                   + as.numeric(sample_data(physeq)$Moisture_Regime)
                   + sample_data(physeq)$pH
                   + sample_data(physeq)$TC_pct
                   + sample_data(physeq)$Sand_pct
                   + sample_data(physeq)$Burned_Unburned
                   , df)
d.adonis

# And looking at year effect on full dataset

physeq = ps.hell
physeq = prune_samples(!is.na(sample_data(physeq)$pH), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$TC_pct), physeq)
physeq = prune_samples(!is.na(sample_data(physeq)$Sand_pct), physeq)

df = as(sample_data(physeq), "data.frame")
d = phyloseq::distance(physeq, method = "bray", weighted=TRUE)  
d.adonis = adonis2(d ~   
                     
                     + sample_data(physeq)$Veg_Comm
                   + as.numeric(sample_data(physeq)$Moisture_Regime)
                   + sample_data(physeq)$pH
                   + sample_data(physeq)$TC_pct
                   + sample_data(physeq)$Sand_pct
                   + sample_data(physeq)$Burned_Unburned 
                   + as.factor(sample_data(physeq)$Years_Since_Fire)
                   , df)
d.adonis


### Is veg_comm from five years post fire better predictor than veg_comm before fire,
# for burned upland sites where there was vegetation transitions?

# First, need to add overstory data for all years.
OverstoryDom = read.csv("../data/Dawe_Veg_Transition_Data/AlluvialPlot_VegSpeciesDominance.csv",row.names=1)
head(OverstoryDom)
colnames(OverstoryDom)[5:7] = c("DomVeg2015","DomVeg2017","DomVeg2019")

# Are site IDs all present?
sample_data(ps.merged)$Site_ID %in% OverstoryDom$Site_ID
sample_data(ps.merged)[!(sample_data(ps.merged)$Site_ID %in% OverstoryDom$Site_ID),]$Burned_Unburned
# Only ones not are the unburned sites. We can assume no transition

# Are all tree data plots in phyloseq?
OverstoryDom$Site_ID %in% sample_data(ps.merged)$Site_ID
# One veg data plot is missing from fungal dataset - drop it
OverstoryDom = OverstoryDom[OverstoryDom$Site_ID %in% sample_data(ps.merged)$Site_ID,]

d = data.frame(sample_data(ps.merged))
d$PsSampleName = sample_names(ps.merged)
d = merge(d,OverstoryDom[,c(1,2,5,6,7)],by="Site_ID",all.x=TRUE)
dim(d)
dim(sample_data(ps.merged))
colnames(sample_data(d))
head(d)
d[1:80,c(1,3,42,76:79)]

# The merging changed the order of the samples
sample_names(ps.merged)[1:10]
d$Sample_ID[1:10]

# Turn d back into sample data and set sample names to sample ID, add back to ps.merged
d=sample_data(d)
sample_names(d) = d$PsSampleName
sample_data(ps.merged) = d

# First check only unburned missing (2019 is missing 4 samples)
sample_data(ps.merged)$Burned_Unburned[is.na(sample_data(ps.merged)$Prefire)]
sample_data(ps.merged)$Burned_Unburned[is.na(sample_data(ps.merged)$DomVeg2015)]
sample_data(ps.merged)$Burned_Unburned[is.na(sample_data(ps.merged)$DomVeg2017)]
sample_data(ps.merged)$Burned_Unburned[is.na(sample_data(ps.merged)$DomVeg2019)]

# Look at those samples
sample_data(ps.merged)[is.na(sample_data(ps.merged)$DomVeg2019),][c(1,2,18,19),]

# With 2015 and 2017 POPUTRE, likely POPUTRE for 2019
sample_data(ps.merged)[is.na(sample_data(ps.merged)$DomVeg2019),][c(1,2,18,19),]$DomVeg2019 = sample_data(ps.merged)[is.na(sample_data(ps.merged)$DomVeg2019),][c(1,2,18,19),]$DomVeg2017

# Add Veg Transition flag
sample_data(ps.merged)$VegTransition = ifelse(sample_data(ps.merged)$Prefire==sample_data(ps.merged)$DomVeg2019,"No","Yes")

# Let's consider the upland burned veg transition sites
ps.norm = transform_sample_counts(ps.merged, function(x) x/sum(x))
ps.hell = transform_sample_counts(ps.merged, function(x) (x/sum(x))^0.5)
physeq = ps.hell

# Test on burned upland samples with veg transitions
physeq = prune_samples(sample_data(physeq)$Land_Class =="Upland", physeq)
physeq = prune_samples(sample_data(physeq)$Burned_Unburned =="Burned", physeq)
#physeq = prune_samples(!is.na(sample_data(physeq)$VegTransition), physeq)
physeq = prune_samples(sample_data(physeq)$VegTransition =="Yes", physeq)

distmat = phyloseq::distance(physeq, method = "bray", weighted=TRUE)
df = data.frame(sample_data(physeq))

d.adonis.pre = adonis2(distmat ~ sample_data(physeq)$Prefire, df)
d.adonis.2019 = adonis2(distmat ~ sample_data(physeq)$DomVeg2019, df)
d.adonis.pre
d.adonis.2019


#### Five years post-fire, are fungal communities in burned sites
# significantly more similar to those of unburned sites?

## If we wanted to compare burned sites to unburned,
# controlling for horizon and veg comm
# in 1 vs 5 years post-fire
# We can see if dissim to unburned decreases between 1 and 5 years post-fire.
ps.1.hell = subset_samples(ps.hell,Years_Since_Fire==1)
ps.5.hell = subset_samples(ps.hell,Years_Since_Fire==5)

ps.1.hell = prune_taxa(taxa_sums(ps.1.hell)>0, ps.1.hell)
ps.5.hell = prune_taxa(taxa_sums(ps.5.hell)>0, ps.5.hell)

# Adjust this for running for 1 or 5 years post fire
#physeq = ps.1.hell
#physeq = ps.5.hell

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
df = df[!(df$Burned_Unburned_1=="Unburned" & df$Burned_Unburned_2=="Unburned"),]
dim(df)


# We want to compare the 1 year and 5 year datasets (run above code for each)
#df.1=df
#df.5=df

# Check normality first
shapiro.test(df.1$Mb_dist)
shapiro.test(df.5$Mb_dist)

# df.5 not normally distributed
# So, use Mann-Whitney U test
wilcox.test(df.5$Mb_dist,df.1$Mb_dist)

d = data.frame(Five=c(df.5$Mb_dist),One=c(df.1$Mb_dist),BSI=df.1$Burn_Severity_Index_Diff)
p = ggplot(d)+geom_point(aes(x=One,y=Five,color=BSI))+geom_abline(slope = 1, intercept=0)
p
# There is not a significant difference between the two years.
# I.e., sites did not become more similar to unburned 5 years post-fire

# If we wanted to test only the burned sites, we can do that, too
df.5.burned = df.5[!(df.5$Burned_Unburned_1=="Unburned" & df.5$Burned_Unburned_2=="Unburned"),]
df.1.burned = df.1[!(df.1$Burned_Unburned_1=="Unburned" & df.1$Burned_Unburned_2=="Unburned"),]
# Check normality first
shapiro.test(df.5.burned$Mb_dist)
shapiro.test(df.1.burned$Mb_dist)

#Mann-Whitney U test
wilcox.test(df.5.burned$Mb_dist,df.1.burned$Mb_dist)
# Still not significant.



### FIGURE 2 ###
# Create ordinations
ord.hell.nmds = ordinate(ps.hell,method="NMDS",distance="bray",k=3,maxit=1000)
ord.norm.nmds = ordinate(ps.norm,method="NMDS",distance="bray",k=3,maxit=1000)

ord = ord.hell.nmds
ps = ps.hell


# I want to connect points that are the same site/sample

# For NMDS
x = data.frame(ord$points)$MDS1
y = data.frame(ord$points)$MDS2
z = data.frame(ord$points)$MDS3

sample_data(ps)$x = x
sample_data(ps)$y = y
sample_data(ps)$z = z

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
sample_names(ps)

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


############################################
# Is burn severity significantly correlated with the degree to which 
# fungal communities changed between one and five years post-fire?
ps = ps.hell
sample_data(ps)$Pairs = paste(sample_data(ps)$Site_ID,sample_data(ps)$Org_or_Min,sep="_")

Dist.mb = as.matrix(vegdist(t(otu_table(ps)), method="bray", type="samples"))
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

# Selecting only comparisons between the same site and horizon, diff years
df = df[df$Pairs != "Different",]
dim(df)

# Want to add back in site info
siteinfo = data.frame(SamDat)%>%
  dplyr::group_by(Pairs,Burned_Unburned,Veg_Comm,Burn_Severity_Index,Org_or_Min)%>%
  summarize(N=n())
head(siteinfo)
tail(siteinfo)
dim(siteinfo)

# Add the data and clean it up
df = merge(df,siteinfo,by="Pairs")%>%
  select(Pairs,Burned_Unburned,Veg_Comm,Burn_Severity_Index,Org_or_Min,Mb_dist)
dim(df)
head(df)
tail(df)

# Plot it
p = ggplot(df)+geom_point(aes(x=Burn_Severity_Index,y=Mb_dist,color=Burned_Unburned))
p = p + facet_wrap(~Org_or_Min)
p


# Check whether relationship is significant for burn
# I.e., did burned vs unburned sites change more between 1 and 5 years?

# May want to logit-transform as in Taylor et al. (2010)
df$Mb_dist_logit = log(df$Mb_dist/(1-df$Mb_dist),base=10)

# Is burn severity significantly correlated with the degree to which 
# fungal communities changed between one and five years post-fire?
df.burnedonly = df %>%
  filter(Burned_Unburned=="Burned")
# Kruskal Wallis test - one way anova by ranks
kruskal.test(Mb_dist~Burn_Severity_Index,data=df.burnedonly)
# BSI is not significant for burned sites
# p = 0.07
kruskal.test(Mb_dist_logit~Burn_Severity_Index,data=df.burnedonly)

# None significant - more severely burned sites did not shift more (or less)
# than lower severity sites between 1 and 5 years post-fire.


############## Are fungal communities more or less related to understory veg comms 1 vs 5 years post-fire?



# Subset by horizon
physeq.O.1 = prune_samples(sample_data(ps.1.norm)$Org_or_Min=="O",ps.1.norm)
physeq.M.1 = prune_samples(sample_data(ps.1.norm)$Org_or_Min=="M",ps.1.norm)

Dist.org.mb.1 = as.matrix(distance(physeq.O.1, method="bray", type="samples"))
Dist.min.mb.1 = as.matrix(distance(physeq.M.1, method="bray", type="samples"))

colnames(Dist.org.mb.1) = sample_data(physeq.O.1)$Site_ID
rownames(Dist.org.mb.1) = sample_data(physeq.O.1)$Site_ID
colnames(Dist.min.mb.1) = sample_data(physeq.M.1)$Site_ID
rownames(Dist.min.mb.1) = sample_data(physeq.M.1)$Site_ID

Dist.veg = read.csv("../../../WB2015/data/Veg_properties/WBNPNWT_Vegetation_Braydistance_2015.csv")
row.names(Dist.veg)=Dist.veg$X
Dist.veg=Dist.veg[,2:dim(Dist.veg)[2]]
colnames(Dist.veg)=row.names(Dist.veg)

OrgSamples.1 = colnames(as.matrix(Dist.org.mb.1))[colnames(as.matrix(Dist.org.mb.1)) %in% colnames(as.matrix(Dist.veg))]
MinSamples.1 = colnames(as.matrix(Dist.min.mb.1))[colnames(as.matrix(Dist.min.mb.1)) %in% colnames(as.matrix(Dist.veg))]

Dist.org.veg.1 = Dist.veg[row.names(Dist.veg) %in% OrgSamples.1,colnames(Dist.veg) %in% OrgSamples.1]
Dist.min.veg.1 = Dist.veg[row.names(Dist.veg) %in% MinSamples.1,colnames(Dist.veg) %in% MinSamples.1]
Dist.org.mb.1 = Dist.org.mb.1[row.names(Dist.org.mb.1) %in% OrgSamples.1,colnames(Dist.org.mb.1) %in% OrgSamples.1]
Dist.min.mb.1 = Dist.min.mb.1[row.names(Dist.min.mb.1) %in% MinSamples.1,colnames(Dist.min.mb.1) %in% MinSamples.1]

plot.1.O=plot(as.dist(Dist.org.mb.1),as.dist(Dist.org.veg.1),xlim=c(0,1), ylim=c(0,1))
plot.1.M=plot(as.dist(Dist.min.mb.1),as.dist(Dist.min.veg.1),xlim=c(0,1), ylim=c(0,1))

Mantel.1.O = vegan::mantel(as.dist(Dist.org.mb.1),as.dist(Dist.org.veg.1), method="spearman", permutations = 999)
Mantel.1.M = vegan::mantel(as.dist(Dist.min.mb.1),as.dist(Dist.min.veg.1), method="spearman", permutations = 999)

Mantel.1.O
Mantel.1.M


### Veg dist 2019 - 5 years post-burn
# Import the dataset
UnderstoryVegComm.2019 = read.csv("../../data/Veg_data/WBNPNWT_Understory-Veg-Cover_2015_2017_2019_11172021.csv")
# Subset to only the 2019 year
UnderstoryVegComm.2019 = UnderstoryVegComm.2019[UnderstoryVegComm.2019$Year==2019,]
dim(UnderstoryVegComm.2019)
UnderstoryVegComm.2019$Site_ID

# S in Veg data is not capitalized; need to replace
# This should fix it
gsub("15s","15S",UnderstoryVegComm.2019$Site_ID)
UnderstoryVegComm.2019$Site_ID=gsub("15s","15S",UnderstoryVegComm.2019$Site_ID)

# Create 2019 veg dist matrix
# Get rid of ID row
row.names(UnderstoryVegComm.2019)=UnderstoryVegComm.2019$Site_ID
UnderstoryVegComm.2019 = UnderstoryVegComm.2019[,3:dim(UnderstoryVegComm.2019)[2]]
rowSums(UnderstoryVegComm.2019)

Dist.veg.2019 = vegdist(as.matrix(UnderstoryVegComm.2019),method="bray")
Dist.veg.2019 = as.matrix(Dist.veg.2019)
row.names(Dist.veg.2019)=row.names(UnderstoryVegComm.2019)
colnames(Dist.veg.2019)=row.names(UnderstoryVegComm.2019)

# Set up the 2019 MB data

# Subset by horizon
physeq.O.5 = prune_samples(sample_data(ps.5.norm)$Org_or_Min=="O",ps.5.norm)
physeq.M.5 = prune_samples(sample_data(ps.5.norm)$Org_or_Min=="M",ps.5.norm)
# Generate dist matrices
Dist.org.mb.5 = as.matrix(distance(physeq.O.5, method="bray", type="samples"))
Dist.min.mb.5 = as.matrix(distance(physeq.M.5, method="bray", type="samples"))
# Retain sample IDs
colnames(Dist.org.mb.5) = sample_data(physeq.O.5)$Site_ID
rownames(Dist.org.mb.5) = sample_data(physeq.O.5)$Site_ID
colnames(Dist.min.mb.5) = sample_data(physeq.M.5)$Site_ID
rownames(Dist.min.mb.5) = sample_data(physeq.M.5)$Site_ID

OrgSamples.5 = colnames(as.matrix(Dist.org.mb.5))[colnames(as.matrix(Dist.org.mb.5)) %in% colnames(as.matrix(Dist.veg.2019))]
MinSamples.5 = colnames(as.matrix(Dist.min.mb.5))[colnames(as.matrix(Dist.min.mb.5)) %in% colnames(as.matrix(Dist.veg.2019))]

Dist.org.veg.5 = Dist.veg.2019[row.names(Dist.veg.2019) %in% OrgSamples.5,colnames(Dist.veg.2019) %in% OrgSamples.5]
Dist.min.veg.5 = Dist.veg.2019[row.names(Dist.veg.2019) %in% MinSamples.5,colnames(Dist.veg.2019) %in% MinSamples.5]
Dist.org.mb.5 = Dist.org.mb.5[row.names(Dist.org.mb.5) %in% OrgSamples.5,colnames(Dist.org.mb.5) %in% OrgSamples.5]
Dist.min.mb.5 = Dist.min.mb.5[row.names(Dist.min.mb.5) %in% MinSamples.5,colnames(Dist.min.mb.5) %in% MinSamples.5]

plot.5.O=plot(as.dist(Dist.org.mb.5),as.dist(Dist.org.veg.5),xlim=c(0,1), ylim=c(0,1))
plot.5.M=plot(as.dist(Dist.min.mb.5),as.dist(Dist.min.veg.5),xlim=c(0,1), ylim=c(0,1))
Mantel.5.O = vegan::mantel(as.dist(Dist.org.mb.5),as.dist(Dist.org.veg.5), method="spearman", permutations = 999)
Mantel.5.M = vegan::mantel(as.dist(Dist.min.mb.5),as.dist(Dist.min.veg.5), method="spearman", permutations = 999)

Mantel.5.O
Mantel.5.M

Mantel.1.O
Mantel.1.M

# Seems like the O horizon is more strongly correlated 5 years post-fire,
# whereas the M horizon is similarly to less well-correlated 5 years later.
# maybe because the O horizon is most affected?


#### FIGURE 3 - Comparing Asco vs Basidios ###
df.phyla = mdf %>%
  dplyr::group_by(Phylum,Sample,Site_ID,Org_or_Min,Years_Since_Fire,Burned_Unburned,Burn_Severity_Index,Severity_Class,Veg_Comm)%>%
  summarize(Relabund = sum(Abundance))

p = ggplot(df.phyla)
p = p + geom_boxplot(aes(x=Severity_Class, y=Relabund,colour=Years_Since_Fire))
p = p + facet_wrap(~Phylum,scales="free")
p


df.AscoBasido.Ratio = df.phyla%>%
  filter(Phylum %in% c("Ascomycota","Basidiomycota"))%>%
  group_by(Sample,Site_ID,Org_or_Min,Years_Since_Fire,Burned_Unburned,Burn_Severity_Index,Severity_Class,Veg_Comm)%>%
  select(Sample,Site_ID,Org_or_Min,Years_Since_Fire,Burned_Unburned,Burn_Severity_Index,Severity_Class,Veg_Comm,Phylum,Relabund)%>%
  arrange(Sample)

df.AscoBasido.Ratio = reshape2::dcast(df.AscoBasido.Ratio, Sample + Site_ID + Org_or_Min+Years_Since_Fire+Burned_Unburned+Burn_Severity_Index+Severity_Class+Veg_Comm ~ Phylum)
df.AscoBasido.Ratio = df.AscoBasido.Ratio %>%
  mutate(ABratio=Ascomycota/Basidiomycota)

# Plotting ratio of ascos to basidios one vs. five years post-fire, for all severity classes
p = ggplot(df.AscoBasido.Ratio[df.AscoBasido.Ratio$ABratio<100,])
p = p + geom_boxplot(aes(x=Severity_Class, y=ABratio,colour=Years_Since_Fire))
p

p = ggplot(df.AscoBasido.Ratio)
p = p + geom_boxplot(aes(x=Severity_Class, y=log(ABratio,10),fill=Years_Since_Fire))
p = p + theme_bw()
p = p + scale_fill_manual(values=c("lightblue","darkblue"))
p = p + ylab("log(Asco:Basidio Ratio)") + xlab("Severity Class")
p = p + guides(fill=guide_legend(title="Years\nSince\nFire"))
#p = p + facet_grid(~Org_or_Min~Years_Since_Fire)
p

lm = lm(log(ABratio,10)~Burned_Unburned+Years_Since_Fire,data=df.AscoBasido.Ratio)
summary(lm)

lm = lm(log(ABratio,10)~Burn_Severity_Index+as.factor(Years_Since_Fire),data=df.AscoBasido.Ratio)
summary(lm)


############## Our next questions will be whether we have the same responders year-to-year, and 
### how their burn response differs.

# We will basically just analyze them separately, and then join by OTUs.

ps.1 = prune_samples(sample_data(ps.merged)$Years_Since_Fire=="1", ps.merged)
ps.5 = prune_samples(sample_data(ps.merged)$Years_Since_Fire=="5", ps.merged)

# Remove no-count taxa
ps.1 = prune_taxa(taxa_sums(ps.1)>0,ps.1)
ps.5 = prune_taxa(taxa_sums(ps.5)>0,ps.5)

# Set phyloseq to non-normalized one, with only non-zero taxa
#ps.biom = ps.1
#ps.biom = ps.5

# Copy phyloseq object (note, not using normalized ps)
otutab = data.frame(as.matrix((otu_table(ps.biom))))
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
results = results[,c(1,3:dim(results)[2])]
colnames(results)[1]="OTU"
results$sigBurn = ifelse(results$adj.P.Val<0.05,1,0.5)
# Creating a variable that sets what will be the alpha values for plotting the results

# Summarizing and saving results
#Results.1 = results
#Results.5 = results

#saveRDS(Results.1,"Results.1.06.08.2024")
#saveRDS(Results.5,"Results.5.06.08.2024")
Results.1 = readRDS("Results.1.06.08.2024")
Results.5 = readRDS("Results.5.06.08.2024")

# Now we want to combine the two dataframes by OTU
# Each row should be an OTU, and the columns in keep are duplicated.
# Filling NAs for responses where not present
dim(Results.1)
dim(Results.5)

Results.5 = Results.5 %>%
  arrange(-logFC)
Results.1 = Results.1 %>%
  arrange(-logFC)

# Merge both years' response results
merged = merge(Results.1,Results.5,by="OTU", all=TRUE, suffixes=c(".OneYear",".FiveYears"))

# Set NAs that weren't tested to zero
# Need to find l2FC for taxa that weren't tested. Should be possible.
logFC.OneYear.calc = mdf %>%
  filter(Years_Since_Fire==1)%>%
  group_by(OTU,Burned_Unburned)%>%
  summarize(MeanRelabund = mean(Abundance))%>%
  mutate(Ratio.OneYear.calc = MeanRelabund[Burned_Unburned=="Burned"]/MeanRelabund[Burned_Unburned=="Unburned"])%>%
  mutate(Ratio.OneYear.calc.inverse = MeanRelabund[Burned_Unburned=="Unburned"]/MeanRelabund[Burned_Unburned=="Burned"])%>%
  filter(Burned_Unburned=="Burned")%>%
  select(OTU,Ratio.OneYear.calc,Ratio.OneYear.calc.inverse)

logFC.FiveYears.calc = mdf %>%
  filter(Years_Since_Fire==5)%>%
  group_by(OTU,Burned_Unburned)%>%
  summarize(MeanRelabund = mean(Abundance))%>%
  mutate(Ratio.FiveYears.calc = MeanRelabund[Burned_Unburned=="Burned"]/MeanRelabund[Burned_Unburned=="Unburned"])%>%
  mutate(Ratio.FiveYears.calc.inverse = MeanRelabund[Burned_Unburned=="Unburned"]/MeanRelabund[Burned_Unburned=="Burned"])%>%
  filter(Burned_Unburned=="Burned")%>%
  select(OTU,Ratio.FiveYears.calc,Ratio.FiveYears.calc.inverse)

# If .calc is infinite, then it means unburned was zero and burned was a number
# If .calc.inverse is infinite, then it means burned was zero and unburned was a number
# If there is NaN, that means that the organism wasn't detected in that year, I believe.

# Could replace with dummy values so long as this is clear
UpperDummy = 6
LowerDummy = -6
logFC.OneYear.calc$logFC.OneYear.calc = 
  ifelse(is.na(logFC.OneYear.calc$Ratio.OneYear.calc),0,
         ifelse(!is.finite(logFC.OneYear.calc$Ratio.OneYear.calc),UpperDummy,
                ifelse(!is.finite(logFC.OneYear.calc$Ratio.OneYear.calc.inverse),LowerDummy,log(logFC.OneYear.calc$Ratio.OneYear.calc,base=2))))

logFC.FiveYears.calc$logFC.FiveYears.calc = 
  ifelse(is.na(logFC.FiveYears.calc$Ratio.FiveYears.calc),0,
         ifelse(!is.finite(logFC.FiveYears.calc$Ratio.FiveYears.calc),UpperDummy,
                ifelse(!is.finite(logFC.FiveYears.calc$Ratio.FiveYears.calc.inverse),LowerDummy,log(logFC.FiveYears.calc$Ratio.FiveYears.calc,base=2))))

# Then want to add these values to our merged wherever we have zeroes
merged = merge(merged,logFC.OneYear.calc[,c("OTU","logFC.OneYear.calc")],by="OTU",all=TRUE)
merged = merge(merged,logFC.FiveYears.calc[,c("OTU","logFC.FiveYears.calc")],by="OTU",all=TRUE)

merged$logFC.OneYear[is.na(merged$logFC.OneYear)]=merged$logFC.OneYear.calc[is.na(merged$logFC.OneYear)]
merged$logFC.FiveYears[is.na(merged$logFC.FiveYears)]=merged$logFC.FiveYears.calc[is.na(merged$logFC.FiveYears)]

# Add not significant alpha values for NAs
merged$sigBurn.OneYear[is.na(merged$sigBurn.OneYear)]=0.5
merged$sigBurn.FiveYears[is.na(merged$sigBurn.FiveYears)]=0.5

# Actually bring mean relative abundance across full dataset
# Could obviously do in burned only but this is just rough anyway
MeanRelabund = mdf %>%
  group_by(OTU)%>%
  summarize(MeanRelabund = mean(Abundance))

merged = merge(merged,MeanRelabund,by="OTU", all.x=TRUE)

# Make a significance column
merged$Sig = ifelse(merged$sigBurn.OneYear==1 & merged$sigBurn.FiveYears==1,"Both years",
                    ifelse(merged$sigBurn.OneYear==1,"One year",
                           ifelse(merged$sigBurn.FiveYears==1,"Five years","Not sig")))

merged$Phylum = ifelse(is.na(merged$Phylum.FiveYears),merged$Phylum.OneYear,merged$Phylum.FiveYears)
merged$Family = ifelse(is.na(merged$Family.FiveYears),merged$Family.OneYear,merged$Family.FiveYears)
merged$Genus = ifelse(is.na(merged$Genus.FiveYears),merged$Genus.OneYear,merged$Genus.FiveYears)
merged = merged %>%
  arrange(-logFC.OneYear)

# Now, to plot significant values, dropping outlier estimates of l2FC
merged.plot = merged%>%
  filter(Sig != "Not sig")

##### FIGURE 4 ######

p = ggplot(merged.plot) + theme_bw()
p = p + geom_hline(yintercept=0,color="black")
p = p + geom_vline(xintercept=0,color="black")
p = p + geom_point(aes(x=logFC.OneYear,y=logFC.FiveYears,fill=Sig,size=MeanRelabund),
                   shape=21,alpha=0.5)
p = p + geom_point(data=merged.plot[merged.plot$Interest==1,],aes(x=logFC.OneYear,y=logFC.FiveYears,size=AveExpr.FiveYears),shape=21,alpha=1,stroke=2)
p = p + scale_fill_manual(values=c("gold","darkblue","lightblue"))
p = p + geom_abline(slope=1,intercept=0,linetype="dashed")
p = p + ylab(expression(Five~year~response~to~burn~(log[2]-fold~change)))
p = p + xlab(expression(One~year~response~to~burn~(log[2]-fold~change)))
p = p + scale_size(guide="none")
p = p + geom_point(data=merged.plot[merged.plot$Genus=="Penicillium",],aes(x=logFC.OneYear,y=logFC.FiveYears,size=MeanRelabund),
                   shape=21,color="black",stroke=1)
p = p + geom_point(data=merged.plot[merged.plot$Genus=="Calyptrozyma",],aes(x=logFC.OneYear,y=logFC.FiveYears,size=MeanRelabund),
                   shape=21,color="red3",stroke=1)
p = p + geom_point(data=merged.plot[merged.plot$Genus=="Coniochaeta",],aes(x=logFC.OneYear,y=logFC.FiveYears,size=MeanRelabund),
                   shape=21,color="orange",stroke=1)
p = p + geom_point(data=merged.plot[merged.plot$Family=="Venturiaceae",],aes(x=logFC.OneYear,y=logFC.FiveYears,size=MeanRelabund),
                   shape=21,color="pink2",stroke=1)
p = p + guides(fill=guide_legend(title=element_text("Significant\nResponse")))
p = p + theme(strip.text.x = element_text(size = 12, face = "italic"))
p

# How many responders in each? (Replace Results.5 or Results.1)
RespondersITS = Results.5 %>%
  filter(!is.na(OTU))%>%
  filter(adj.P.Val<0.05)%>%
  arrange(-logFC)%>%
  mutate(PosNeg = ifelse(logFC>0,"Pos",
                         ifelse(logFC<0,"Neg","Other")))%>%
  group_by(PosNeg)%>%
  summarize(Total=n())
head(RespondersITS)  


# How many Both Positive responders?
BothPos = merged.plot %>%
  filter(logFC.OneYear>0 &
           logFC.FiveYears>0 &
              Sig == "Both years")
dim(BothPos)

BothNeg = merged.plot %>%
  filter(logFC.OneYear<0 &
           logFC.FiveYears<0 &
           Sig == "Both years")
dim(BothNeg)


### Penicillium responders and FIGURE 5 ####
Gen.plot = merged%>%
  filter(Sig != "Not sig")%>%
  filter(Genus.OneYear == "Penicillium")

p = ggplot(Gen.plot) + theme_bw()
p = p + geom_hline(yintercept=0,color="black")
p = p + geom_vline(xintercept=0,color="black")
p = p + facet_wrap(~Genus.OneYear)
p = p + geom_point(aes(x=logFC.OneYear,y=logFC.FiveYears,fill=Sig,size=MeanRelabund),
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

# Wondering if some of these are really the same OTU, though.
Penns = Biostrings::DNAStringSet(Gen.plot$OTU)
aln = DECIPHER::AlignSeqs(Penns, processors = 18)
DECIPHER::BrowseSeqs(aln, highlight=0)

d = DECIPHER::DistanceMatrix(aln, processors = NULL, includeTerminalGaps = FALSE,
                             penalizeGapLetterMatches = TRUE,
                             penalizeGapGapMatches = FALSE)

# Cluster at 3% similarity
clusters = DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03,
  processors = NULL)
clusters

# Some OTUS cluster at 97%
d<0.03

Gen.plot$cluster = clusters$cluster

p = ggplot(Gen.plot) + theme_bw()
p = p + geom_hline(yintercept=0,color="black")
p = p + geom_vline(xintercept=0,color="black")
p = p + facet_wrap(~Genus.OneYear)
p = p + geom_point(aes(x=logFC.OneYear,y=logFC.FiveYears,fill=as.factor(cluster),size=MeanRelabund),
                   shape=21,alpha=0.5)
p = p + geom_point(data=merged.plot[merged.plot$Interest==1,],aes(x=logFC.OneYear,y=logFC.FiveYears,size=AveExpr.FiveYears),shape=21,alpha=1,stroke=2)
#p = p + scale_fill_manual(values=c("gold","darkblue","lightblue"))
p = p + geom_abline(slope=1,intercept=0,linetype="dashed")
p = p + ylab(expression(Five~year~response~to~burn~(log[2]-fold~change)))
p = p + xlab(expression(One~year~response~to~burn~(log[2]-fold~change)))
p = p + scale_size(guide="none")
p = p + guides(fill=guide_legend(title=element_text("Significant\nResponse")))
p = p + theme(strip.text.x = element_text(size = 12, face = "italic"))
p

# The 5 year responders do not seem to be mis-clustered 1-year responders, as they are different colours here

# Let's plot a subset of the Penicillium OTUs of interest from this study.

# Make a list of the Penicillium OTUs in order of mean relative abundance in burned samples
AbundPenn = mdf%>%
  filter(Genus=="Penicillium")%>%
  group_by(OTU)%>%
  summarize(MeanRelabundBurned = mean(Abundance[Burned_Unburned=="Burned"]),
            MeanRelabundUnburned = mean(Abundance[Burned_Unburned=="Unburned"]))%>%
  arrange(-MeanRelabundBurned)

# Add an OTU number to help organize them later
AbundPenn$OTUNumber = 1:dim(AbundPenn)[1]
AbundPenn

# Make a list of the Penicillium burn responders
PennResponders = merged %>%
  filter(Genus=="Penicillium")%>%
  filter(Sig != "Not sig")%>%
  select(OTU,Sig)


#### Make a list of OTUs within 99% ID of a Day taxon (as determined externally using BLAST)
DayTaxa = c("TCAGCGGGTAGTCCCGCCTGATTTGAGGTCAAAATGTCAAGAGTTGTCCAAGTTAATGGACGGTTATGAGCTGAACCCCATGTAAAGCTGCTTCACGGCAATGGCGTAGATAATTATCACACCAGAGACGGTCCACAAAGGTTCCGCTAATGCATTTAAGGAGAGCCGACCTCTGAAGAAGCCGGCAACCCCCACATCCAAGCCTACGCCGACTCGTAAAAGCTGGCAAGGTTGAGAATTTAATGACACTCAAACAGGCGTGCCCCTCGGAATACCAAGGGGCGCAAGGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCGAGAGCCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGAAAAGGATGATGGGTTGTCGGCTGGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACTCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCGGGGGGACGGAGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTCTGCAATTCACATTAGTTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGAAAAGGATGATTGGTTGTCGGCTGGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACTCGGTGCCGCCGCTGCCTTTCGGACCCGTCCCCGGGGGGACGGAGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTCTGCAATTCACATTAGTTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGATAAAATTTTGGGTTGATCGGCAAGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACGCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCCGGAATCGGAGGACGGGGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTCTGCAATTCACATTACGTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGATAAAAATTTGGGTTGATCGGCAAGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACGCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCCGGAATCGGAGGACGGGGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTCTGCAATTCACATTACGTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGATAAAAATTTGGGTTGATCGGCAAGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACGCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCCGGAATCGGAGGACGGGGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTTGCAATTCACATTACGTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGATAAAAATTTGGGTTGATCGGCAAGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACGCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCCGGAATCGGAGGACGGGGCCCAACACACAAGCCGGGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTTGCAATTCACATTACGTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGATAAAAATTTGGGTTGATCGGCAAGCGCCGGCCGGGCCTACAGAGCGGGTGACAAAGCCCCATACGCTCGAGGACCGGACGCGGTGCCGCCGCTGCCTTTCGGGCCCGTCCCCCGAAATCGGAGGACGGGGCCCAACACACAAGCCGTGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTTGCAATTCACATTACGTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTACTCTTACCTGATCCGAGGTCAACCTTGTGTTAAAGGGGGTCTTTAACGGCTGGAACTCGCTGTAACTCCCTAGCGATGATAGTGTACTACTACGCTCGGAGTTGTAGCGAGCCCGCCACTTCTTTTCAGGGCCTACGGCAGCCGTAGGGCCCCAACACCAAGCAGGGCTTGAGGGTTGAAATGACGCTCGAACAGGCATGCCCGCTAGAATACTAGCGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAGCCA",
            
            "TCAGCGGGTACTCTTACCTGATCCGAGGTCAACCTTGTGTTAAAGGGGGTCTTTAACGGCTGGAACTCGCTGTAACTCCCTAGCGACGATAGTGTGCTACTACGCTCGGAGTTGTAGCGAGCCCGCCACTTCTTTTCAGGGCCTACGGCAGCCGTAGGGCCCCAACACCAAGCAGGGCTTGAGGGTTGAAATGACGCTCGAACAGGCATGCCCGCTAGAATACTAGCGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAGCCA",
            
            "TCAGCGGGTATCCCTACCTGATCCGAGGTCAACCTGGAAAAAAAAGGTTGGTAAGTCGGCAGACTATCGGCCGTTCCTACTAGAGCGGAGTGACGAAGCCCCATACGCTCGAGGACCGGACGCGAAGTCGCCGCTGCCTTTCGGGCCCGTCCCCCCCGGGAGGGGAGGACGGCGACCCAACACACAAGCCGGGCTTGAGGGCAGCAATGACGCTCGGACAGGCATGCCCCCCGGAATACCAGGGGGCGCAATGTGCGTTCAAAGACTCGATGATTCACTGAATTCTGCAATTCACATTAGTTATCGCATTTCGCTGCGTTCTTCATCGATGCCGGAACCA",
            
            "TCAGCGGGTATCCCTACCCGATCCGAGGTCAACCATAGAAATTTAGGGGTTGATGGCAAGCATCCACCAGGACCCTGTAGCGAGAAGTATTACTACGCTTAGAGCCAGATGGCACCGCCACTGATTTTAAGGGCTGCCGGTAACAGCAGGCCCCAACACCAAGCTGGGCTTGAGGGGTTATAATGACGCTCGAACGGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAACCA",
            
            "TCAGCGGGTATCCCTACCCGATCCGAGGTCAACCATAGAAATTTAGGGGTTGATGGCAAGCATCCACCGGGACCCTGTAGCGAGAAGTATTACTACGCTTAGAGCCAGATGGCACCGCCACTGATTTTAAGGGCTGCCGGTAACAGCAGGCCCCAACACCAAGCTGGGCTTGAGGGGTTATAATGACGCTCGAACGGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAACCA",
            
            "TCAGCGGGTATCCCTACCCGATCCGAGGTCAACCATAGAAATTTAGGGGTTGATGGCAAGTATCCACCAGGACCCTGTAGCGAGAAGTATTACTACGCTTAGAGCCAGATGGCACCGCCACTGATTTTAAGGGCTGCCGGTAACAGCAGGCCCCAACACCAAGCTAGGCTTGAGGGGTTATAATGACGCTCGAACGGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAACCA",
            
            "TCAGCGGGTATCCCTACCCGATCCGAGGTCAACCATAGAAGTTTAGGGGTTGATGGCAAGCATCCACCAGGACCCTGTAGCGAGAAGTATTACTACGCTTAGAGCCAGATGGCACCGCCACTGATTTTAAGGGCTGCCGGCAACAGCAGGCCCCAACACCAAGCTGGGCTTGAGGGTTATAATGACGCTCGAACGGGCATGCCCCTCGGAATACCAAGGGGCGCAATGTGCGTTCAAAGATTCGATGATTCACTGAATTCTGCAATTCACATTACTTATCGCATTTCGCTGCGTTCTTCATCGATGCCAGAACCA")


#### Get Penicillium taxa

df.penn = mdf %>%
  filter(Genus=="Penicillium")

# Add their Abundance-ordered OTU number
df.penn = merge(df.penn,AbundPenn[,c("OTU","OTUNumber")],by="OTU",all=TRUE)

# Annotate 
df.penn = df.penn%>%
  mutate(Pairs=paste(Site_ID,Org_or_Min,sep="_"))%>%
  group_by(Sample,Pairs,Burned_Unburned,Severity_Class,pH,TC_pct,Burn_Severity_Index,Org_or_Min,Veg_Comm,Genus,Species,OTU,guildHighProb,Years_Since_Fire,OTUNumber)%>%
  summarize(Relabund = sum(Abundance))%>%
  group_by(OTU,Org_or_Min)%>%
  mutate(MeanRelabund = mean(Relabund))%>%
  mutate(MaxRelabund = max(Relabund))%>%
  filter(OTU %in% AbundPenn$OTU)%>%
  mutate(YSFbox = ifelse(Years_Since_Fire==1,0.5,5.5))%>%
  mutate(OTULab = ifelse(is.na(Species),"Penicillium_sp.",Species))%>%
  mutate(OTULab = ifelse(OTULab=="Penicillium_sp.",paste("Penicillium_sp.",OTUNumber),OTULab))%>%
  mutate(OTULab = ifelse(OTU %in% PennResponders$OTU,paste(OTULab,"*"),OTULab))%>%
  mutate(OTULab = ifelse(OTU %in% DayTaxa,paste(OTULab,"+"),OTULab))

#Find Day Penicillium and check they were annotated correctly above
DayPenn = df.penn %>%
  filter(OTU %in% DayTaxa)%>%
  group_by(OTULab)%>%
  summarize(n(),)
DayPenn

####

# Gen.plot had our Penicillium responders (either year, pos or neg) from above
# Which of those were in the Day paper?
Gen.plot$OTU %in% DayTaxa

# Here we're choosing a few OTUs to illustrate different responses to fire
PennToPlot = c(Gen.plot$OTU[c(1,3,9)],AbundPenn$OTU[1:2])

df.penn.plot = df.penn%>%
  filter(OTU %in% PennToPlot)

df.penn.plot$Response = ifelse(df.penn.plot$OTU %in% Gen.plot$OTU,"Fire responder","Non-responder")
df.penn.plot$DayTaxa = ifelse(df.penn.plot$OTU %in% DayTaxa,"Day OTU","Not Day OTU")


p = ggplot(df.penn.plot)
p = p + geom_boxplot(aes(x=Severity_Class,y=(Relabund),fill=Severity_Class)) + theme_bw()
p = p + facet_wrap(~Years_Since_Fire~Response*DayTaxa*OTU*Species,scales="free",ncol=5)
p = p + xlab("Severity Class")
p = p + ylab("Relative Abundance")
p = p + scale_fill_manual(values=c("darkred","red3","orange","gold"))
p

### Saving csv of strong responders abundant OTUs ###
cutoff = 1
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
OTU.maxes = data.frame(OTU=taxa_names(ps.norm),
                       MaxAbund = apply(t(otu_table(ps.norm)),2,function(x) max(x)),
                       MeanAbund = apply(t(otu_table(ps.norm)),2,function(x) mean(x)))
merged.cat = plyr::join(merged.cat,OTU.maxes,by="OTU",type="left")

# Ok, got them in, now another filter
cutoff = 0.005
merged.cat = merged.cat %>%
  filter(MaxAbund>=cutoff)%>%
  arrange(Response_Category)%>%
  arrange(-logFC.OneYear)
merged.cat

#write.csv(merged.cat,"StrongAbundantResponders.ITS.28.03.2025.csv")


### FIGURE 6 ###

#### Looking at FunGuilds info

#### Looking at FunGuilds info

df = mdf %>%
  group_by(Sample,Burned_Unburned,Burn_Severity_Index,Org_or_Min,Veg_Comm,trophicMode,trophicModeHighProb,Years_Since_Fire)%>%
  summarize(Relabund = sum(Abundance))


# Plotting all samples
p = ggplot(df)
p = p + geom_bar(aes(x=Sample,y=Relabund,fill=trophicModeHighProb),stat="identity", position="stack")
#p = p + geom_bar(aes(x=Sample,y=Relabund,fill=trophicMode),stat="identity", position="stack")
p = p + facet_wrap(~Org_or_Min~Veg_Comm~Years_Since_Fire~Burned_Unburned,scales="free",ncol=8)
p

# Plotting by trophic mode
p = ggplot(df)
p = p + geom_boxplot(aes(x=Org_or_Min,y=Relabund,colour=Burned_Unburned))
p = p + facet_wrap(~trophicModeHighProb~paste(Org_or_Min,Years_Since_Fire),scales="free",ncol=4)
p

# Figures for Sapros, Pathos, and Symbios
OM.labs = c("Mineral", "Organic")
names(OM.labs) = c("M", "O")

p = ggplot(df[df$trophicModeHighProb=="Saprotroph",])
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_grid(~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of Saprotrophs")
p

aov.sapro.O = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Saprotroph" & df$Org_or_Min=="O",])
summary(aov.sapro.O)
TukeyHSD(aov.sapro.O, conf.level=.95)

aov.sapro.M = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Saprotroph" & df$Org_or_Min=="M",])
summary(aov.sapro.M)
TukeyHSD(aov.sapro.M, conf.level=.95)


p = ggplot(df[df$trophicModeHighProb=="Pathotroph",])
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_grid(~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of Pathotrophs")
p

aov.patho.O = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Pathotroph" & df$Org_or_Min=="O",])
summary(aov.patho.O)
TukeyHSD(aov.patho.O, conf.level=.95)

aov.patho.M = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Pathotroph" & df$Org_or_Min=="M",])
summary(aov.patho.M)
TukeyHSD(aov.patho.M, conf.level=.95)

p = ggplot(df[df$trophicModeHighProb=="Symbiotroph",])
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_grid(~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of Symbiotrophs")
p

aov.symb.O = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Symbiotroph" & df$Org_or_Min=="O",])
summary(aov.symb.O)
TukeyHSD(aov.symb.O, conf.level=.95)

aov.symb.M = aov(Relabund~Years_Since_Fire*Burned_Unburned,data=df[df$trophicModeHighProb=="Symbiotroph" & df$Org_or_Min=="M",])
summary(aov.symb.M)
TukeyHSD(aov.symb.M, conf.level=.95)

### Illustrating which taxa are driving responses ###


# See which taxa are driving response in Saprotrophs
cutoff = 0.0005
levels(as.factor(mdf$guildHighProb))
df.sapro = mdf %>%
  group_by(Sample,Burned_Unburned,Burn_Severity_Index,Org_or_Min,Veg_Comm,Genus,guildHighProb,trophicModeHighProb,Years_Since_Fire)%>%
  summarize(Relabund = sum(Abundance))%>%
  filter(trophicModeHighProb=="Saprotroph")%>%
  group_by(Org_or_Min,Genus)%>%
  mutate(MeanRelabund = mean(Relabund))%>%
  mutate(MaxRelabund = max(Relabund))%>%
  filter(MeanRelabund > cutoff)

p = ggplot(df.sapro)
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_wrap(~Genus~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of abundant Saprotrophs")
p

# See which ectos are driving response in Symbiotrophs
cutoff = 0.01

df.ectos = mdf %>%
  group_by(Sample,Burned_Unburned,Burn_Severity_Index,Org_or_Min,Veg_Comm,Genus,guildHighProb,Years_Since_Fire)%>%
  summarize(Relabund = sum(Abundance))%>%
  filter(guildHighProb=="Ectomycorrhizal")%>%
  group_by(Org_or_Min,Genus)%>%
  mutate(MeanRelabund = mean(Relabund))%>%
  mutate(MaxRelabund = max(Relabund))%>%
  filter(MeanRelabund > cutoff)

p = ggplot(df.ectos)
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_wrap(~Genus~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of Ectomycorrhizae")
p

# See which endophytes are driving response in Symbiotrophs

cutoff = 0.005

df.endos = mdf %>%
  group_by(Sample,Burned_Unburned,Burn_Severity_Index,Org_or_Min,Veg_Comm,Genus,guildHighProb,Years_Since_Fire)%>%
  summarize(Relabund = sum(Abundance))%>%
  filter(guildHighProb=="Endophyte")%>%
  group_by(Org_or_Min,Genus)%>%
  mutate(MeanRelabund = mean(Relabund))%>%
  mutate(MaxRelabund = max(Relabund))%>%
  filter(MeanRelabund > cutoff)

p = ggplot(df.endos)
p = p + theme_bw()
p = p + geom_boxplot(aes(x=Years_Since_Fire,y=Relabund,fill=Burned_Unburned))
p = p + scale_fill_manual(values=c("red3","gold"))
p = p + facet_wrap(~Genus~Org_or_Min,scales="free",labeller=labeller(Org_or_Min=OM.labs))
p = p + xlab("Years Since Fire") + ylab("Relative Abundance of Endophytes")
p


