# Creating co-occurrence networks

library(phyloseq)
library(ggplot2)
library(reshape)
library(plyr)
library(dplyr)
library(vegan)
library(ade4)
library(wesanderson)
library(igraph)
library(tidyr)
library(threejs)
library(htmlwidgets)
library(RColorBrewer)

# Set working directory
setwd("~/Box/WhitmanLabMaster/WhitmanLab/Projects/WoodBuffalo/NWT-Fire-Microbes/WB2019/code")

# Load ps object
ps.merged = readRDS("ps.merged")
ps.merged

# Add pH and pct C data from 2019
# Importing new pH values for 2019 soils
pH.2019 = read.csv("../data/Soils_data/Summaries/WB2019 pH means_04-30-21.csv")

# Make sure same sample IDs are present
as.factor(pH.2019$Sample_ID) %in% sample_data(ps.merged)$Sample_ID
sample_data(ps.merged)$Sample_ID %in% as.factor(pH.2019$Sample_ID)
# Yes, good, all sample IDs in our phyloseq object are found in the pH table
# we have a few measured pH values that aren't in our dataset - that's fine


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


# Normalize counts
ps.merged.norm = transform_sample_counts(ps.merged, function(x) x/sum(x))


# We want to remove the least abundant taxa, to reduce computational load and to not bother with low-abundance taxa
# We could interpret this as being, across all samples, at least X abundant
# That way, something low abundance but widely present, or very high abundance, will get included

# The cutoff is, of course, somewhat arbitrary - present across all samples at 0.005 seems quite low, though.
# You could also likely argue that it should be different across fungi/bact/plants, but how to justify the precise number?
cutoff = 0.005

ps.mini = prune_taxa((taxa_sums(ps.merged.norm) > cutoff), ps.merged.norm)
ps.mini

# Start with just 2019 data, normalized
ps.2019 = subset_samples(ps.mini,Project_ID=="WB19")
ps.2019 = prune_taxa((taxa_sums(ps.2019) > 0), ps.2019)
ps.2019

# Can compare to 2015 data from same sites
ps.2015 = subset_samples(ps.mini,is.na(Project_ID))
ps.2015 = prune_taxa((taxa_sums(ps.2015) > 0), ps.2015)
ps.2015

# Okay. After trimming the full dataset to a cutoff, and then removing the zero-count taxa,
# That's pretty close total numbers between the two datasets - 2866 and 2926. Diff of ~50 (2%)
# Ok to proceed.

# Set phyloseq object to one of these.
ps=ps.2019
#ps=ps.2015

#### Matrix with no noise added ####

# Creating a cutoff function to collect the Spearman rho cutoff and the fraction of OTUs that are included
# in largest cluster
cutoff_function = function(cutoff,ps){
  adjacency_matrix = as.matrix(ps)
  adjacency_matrix[abs(adjacency_matrix)<cutoff] = 0
  adjacency_matrix[abs(adjacency_matrix)>cutoff] = 1
  am = graph.adjacency(adjacency_matrix)
  c = clusters(am)
  MaxSize = c$csize[1]
  return(c(cutoff,MaxSize/dim(adjacency_matrix)[1]))
}

# Test it against this range of values
# Started out with wide range from 0 to 1, narrowed in on ranges of interest
# to pinpoint optimal rho cutoff value
inputs = seq(0.35,0.45,0.01)

# Setting up the data frame
Rep=c()
cutoff=c()
PercentIncl=c()
df.no=data.frame(Rep,cutoff,PercentIncl)

network_function = function(ps){
  for (i in 1:1){
    
    # Record iteration
    Rep = i
    
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Calculate spearman correlations
    
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs, cutoff_function, ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.no)
    df.no=rbind(df.no,df)
  }
  colnames(df.no) = c("Rep","Cutoff","PctIncl")
  return(df.no)
}

df.n = network_function(ps)


# Plot percent of OTUs included in largest cluster with increasing Spearman cutoffs
p = ggplot(df.n)
p = p + geom_point(aes(x=Cutoff,y=PctIncl, color=Rep))
p
# Dropoff happens around 0.42

#### Matrix with just noise added ####

## To calculate minimum distance between taxa across the whole dataset
# Basically, it will be 1/[the most total sequences across all samples]
mindist = min(1/max(sample_sums(ps)))
# Across all datasets, the minimum distance between taxa within a sample
delta = mindist*10^-2
# Just to be generous

# Setting up the data frame
Rep=c()
cutoff=c()
PercentIncl=c()
df.list=data.frame(Rep,cutoff,PercentIncl)

noise_function = function(ps){
  # Running it multiple times
  for (i in 1:10){
    
    # Record iteration
    Rep = i
    
    # Making matrix of random values that are less than smallest difference between two samples
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Add the noise to the matrix
    
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Calculate spearman correlations
    
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs,cutoff_function,ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.list)
    df.list=rbind(df.list,df)
  }
  colnames(df.list) = c("Rep","Cutoff","PctIncl")
  return(df.list)
}

df.e = noise_function(ps)

p = ggplot(df.e,aes(x=Cutoff,y=PctIncl, color=Rep))
p = p + geom_point()
p
# As we can see, including random error to break ties
# reduces the number of taxa included in the network.
# Now the dropoff starts a bit earlier

#### Matrix with noise and permuations added ####

# Setting up the data frame
Rep=c()
Cutoff=c()
PctIncl=c()
df.list.permute=data.frame(Rep,Cutoff,PctIncl)


permute_noise_function = function(ps){
  # Running it multiple times
  for (i in 1:30){
    
    # Record iteration
    Rep = i
    
    # Making matrix of random values as above
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Add the noise to the matrix
    
    M = matrix(otu_table(ps))
    # Get the OTU table as matrix
    
    n = dim(M)[1]*dim(M)[2]
    # number of elements in matrix
    
    M.new = matrix(base::sample(M, n, replace=TRUE),nrow=dim(otu_table(ps))[1])
    # Sample the matrix with replacement; make new random matrix
    
    ps.dist = cor(M.new, use="everything", method="spearman") 
    # Calculate spearman correlations
    
    # Apply the function to the range of values and turn it into a dataframe
    df = t(sapply(inputs, cutoff_function, ps=ps.dist))
    df = data.frame(df)
    df = data.frame(i,df)
    colnames(df) = colnames(df.list.permute)
    df.list.permute=rbind(df.list.permute,df)
  }
  colnames(df.list.permute) = c("Rep","Cutoff","PctIncl")
  return(df.list.permute)
}

df.p = permute_noise_function(ps)

p = ggplot(df.p,aes(x=Cutoff,y=PctIncl, color=Rep))
p = p + geom_point()
p
# As we can see, permuting the matrix randomly still makes networks.
# We must choose a rho cutoff that is above the random threshold.

# Join the three types of runs together
df.n$Set = "Nothing"
df.p$Set = "Error+Permute"
df.e$Set = "Error"
df.full = rbind(df.n,df.e,df.p)

df.full.2019=df.full
#df.full.2015=df.full
df.full = df.full.2015

max(df.full[df.full$Cutoff==0.43  & df.full$Set=="Error+Permute",]$PctIncl)
# Choosing a cutoff value that keeps the max below 0.01 (1%) in the error+permuted dataset
# For 2019, 0.43 cutoff gives <0.01 max
# For 2015, 0.42 cutoff gives <0.01 max, but 0.43 is closer to 2019, so can go with 0.43 for both of them.
# Will use those cutoffs for comparability.

Cutoff=0.43

p = ggplot(df.full,aes(x=Cutoff,y=PctIncl, color=Set)) + geom_point() #+ ylim(c(0,0.1))
p

#The above code basically does what we want - calculates how large the largest cluster is 
#for each of the rho cutoff values, with no alterations, adding random error, and adding 
#random error plus rearranging the matrix randomly to result in no expected correlations. 
#Thus, we can see where the rho cutoff values fall for all options. We should choose a cutoff 
#that is higher than where the permuted one reaches ~1% of the OTUs being included in the 
#largest cluster. We can then look at that cutoff in the actual (with added random error) data.

# To be really confident in this, any chosen cutoff should be examined to make sure small variations in its value do not affect ecological conclusions

# Now we should be ready to use that cutoff to calculate the various parameters of the system.
# Now, to test our confidence in the overall network developed, we will caculate a series of 
# parameters about the network, and determine their null distributions.

# First, we will use the derived network, based on the cutoff that we chose above.

#### Matrix with just noise added ####

# Chosen based on above calculations
rho.2015 = 0.43
rho.2019 = 0.43
rho = rho.2015

matrix_parameters = function(ps,rho){
  df.reps=data.frame(m=c(),n=c(),k=c(),apl=c(),c=c())
  #Starting out running it 100 times
  for (i in 1:100){
    
    # Record iteration
    Rep = i
    
    # Making matrix of random values
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Add the noise to the matrix
    
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Calculate spearman correlations
    
    adjacency_matrix = as.matrix(ps.dist)
    # Turns Spearman into matrix
    
    adjacency_matrix[abs(adjacency_matrix)<rho] = 0
    adjacency_matrix[abs(adjacency_matrix)>rho] = 1
    adjacency_matrix = adjacency_matrix[rowSums(adjacency_matrix)>0,]
    # Set correlations below rho cutoff to 0, and above to 1
    # Remove all taxa that have no correlations above the cutoff.
    
    adjacency_matrix = graph.adjacency(adjacency_matrix, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Create adjacency matrix from correlation matrix
    
    edge_list = get.edgelist(adjacency_matrix)
    # Gets edge list from adjacency matrix
    
    N = graph_from_data_frame(edge_list,directed=FALSE)
    # Creates the network
    
    m = ecount(N)
    # Number of edges
    n = vcount(N)
    # Number of nodes
    k = 2*m/n
    # Average degree
    apl = mean_distance(N,directed=FALSE)
    # Average path length
    c = transitivity(N)
    # Clustering coefficient of whole graph
    
    # Apply the function to the range of values and turn it into a dataframe
    df = data.frame(m,n,k,apl,c)
    colnames(df) = colnames(df.reps)
    df.reps=rbind(df.reps,df)
  }
  colnames(df.reps)=c("m","n","k","apl","c")
  return(df.reps)
}

df.reps = matrix_parameters(ps,rho)


k.ave = mean(df.reps$k)
n.ave = mean(df.reps$n)
edges = k.ave * n.ave / 2
p = edges * 2 / (n.ave*n.ave)
# Probability that any pair of verticies is connected
# Thus we can generate a random matrix that preserved the degree of the network

ER_function = function(k.ave,n.ave,edges,p){
  df.ER=data.frame(Rep=c(),m=c(),n=c(),k=c(),apl=c(),c=c())
  
  for (i in 1:100){
    Rep=i
    ER.thresh = matrix(runif(round(n.ave)*round(n.ave)),ncol=round(n.ave))
    # Makes a random uniform matrix with the average number of nodes
    ER.thresh[upper.tri(ER.thresh,diag=FALSE) & ER.thresh < (1-p)] = 0
    ER.thresh[upper.tri(ER.thresh,diag=FALSE) & ER.thresh > (1-p)] = 1
    # Set anything in the upper triangle that is below the threshold (1-p) to zero and anything above to 1
    ER.thresh[lower.tri(ER.thresh,diag=FALSE)] = 0
    # Reflect the upper triangle to the lower triangle
    diag(ER.thresh) = 0
    # Set the diagonal to 0 (no self-connections)
    # This is the ER random thresholded matrix
    
    adjacency_matrix = graph.adjacency(ER.thresh, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Create adjacency matrix from correlation matrix
    
    edge_list = get.edgelist(adjacency_matrix)
    # Gets edge list from adjacency matrix
    
    N = graph_from_data_frame(edge_list,directed=FALSE)
    # Creates the network
    
    m = ecount(N)
    # Number of edges
    n = vcount(N)
    # Number of nodes
    k = 2*m/n
    # Average degree
    apl = mean_distance(N,directed=FALSE)
    # Average path length
    c = transitivity(N)
    # Clustering coefficient of whole graph
    
    # Apply the function to the range of values and turn it into a dataframe
    df = data.frame(Rep,m,n,k,apl,c)
    colnames(df) = colnames(df.ER)
    df.ER=rbind(df.ER,df)
  }
  colnames(df.ER)=c("Rep","m","n","k","apl","c")
  return(df.ER)
}

ER = ER_function(k.ave,n.ave,edges,p)

head(ER)
head(df.reps)

plot(density(ER$apl),xlim=c(2.5,4.5),col="red")
lines(density(df.reps$apl))

plot(density(ER$c),xlim=c(0,0.5),col="red")
lines(density(df.reps$c))

# The red traces represent the distribution of values for a random network with similar properties
# Since our true values lie outside of these, generally, this is good - if the network is randomly all connected,
# Average path length would be low - everything is connected. More structure increases APL.

# Since our network statistics seem to lie well outside of the null distributions 
# (more clustering and longer path length), it seems like there is a good chance that 
# our network is not just representing random noise. The clustering coefficient is around 0.3 - 0.34; 
# the Connor paper was 0.38.

# Now we are ready to determine a consensus network. 
# Connor et al. ran 2000 network simulations with added noise and 
# using their permutationally-verified rho cutoff, and included edges 
# that were present in 90% of their simulations.
# We will do 1000 at 95%.

#### Running consensus network ####

consensus_network = function(ps,rho){
  EL = data.frame(Rep=c(),X1=c(),X2=c())
  
  # Running it 1000 times
  for (i in 1:1000){
    
    # Record iteration
    Rep = i
    
    # Making matrix of random values
    b = delta/1000
    E = replicate(dim(otu_table(ps))[1], rnorm(dim(otu_table(ps))[2]))
    E = 2*b*E
    E = -b + E
    
    otu_table(ps) = otu_table(as.matrix(otu_table(ps))+abs(t(E)),taxa_are_rows=FALSE)
    # Add the noise to the matrix
    
    ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
    # Calculate spearman correlations
    
    adjacency_matrix = as.matrix(ps.dist)
    # Turns Spearman into matrix
    
    adjacency_matrix[abs(adjacency_matrix)<rho] = 0
    adjacency_matrix[abs(adjacency_matrix)>rho] = 1
    adjacency_matrix = adjacency_matrix[rowSums(adjacency_matrix)>0,]
    # Set correlations below rho cutoff to 0, and above to 1
    # Remove all taxa that have no correlations above the cutoff.
    
    adjacency_matrix = graph.adjacency(adjacency_matrix, diag=FALSE, mode=c("undirected"), weighted=TRUE)
    # Create adjacency matrix from correlation matrix
    
    edge_list = get.edgelist(adjacency_matrix)
    # Gets edge list from adjacency matrix
    
    edge_list = data.frame(Rep,edge_list)
    
    EL = rbind(EL,edge_list)
    
  }
  EL$X3=""
  EL$X3 = apply(EL,1,function(x) paste(sort(c(paste(x[2]),paste(x[3])))[1],sort(c(paste(x[2]),paste(x[3])))[2]))
  # The issue is that the pairing "OTU1 OTU2" can also be written "OTU2 OTU1"
  # If we want to match the pairings, we need to make sure they are always written the same way
  # To do this, we take X1 and X2 (the two nodes) and put them in a consistent order (~alphabetical)
  # (The above code is quite a bit (1-2x) faster than a paralellized mapply or a for loop to do this.)
  return(EL)
}

EL = consensus_network(ps,rho)
#saveRDS(EL,"EL.201X.x-x-x")
#saveRDS(EL,"EL.2015.06-04-2020")
EL = readRDS("EL.2019.06-04-2020")
#EL = readRDS("EL.2015.06-04-2020")


# Finds how many instances of each pairing there are

Common = EL %>%
  dplyr::group_by(X3)%>%
  dplyr::summarize(n())

# We want those that are present 95% of the time
# For 1000 runs, that's n>=950 or higher.
# For 100 runs, that's n>=95 or higher.
#Consensus = Common[Common[,2]>=1900,]
Consensus = Common[Common[,2]>=950,]

dim(EL)
dim(Common)
dim(Consensus)


# Divide the 95%-present paired responders back into two columns
# This is our consensus edge list
edge_list_consensus = separate(Consensus,1,into=c("X1","X2"),sep=" ")
edge_list_consensus = edge_list_consensus[,1:2]

edge_list = edge_list_consensus
# Setting edge list to the consensus value
# Then running through all previous graphing code

colnames(edge_list)=c("A","B")
edge_list=data.frame(edge_list)
edge_list$CorVal=""
edge_list$CorSign=""

# Want to add correlation values and colors for the consensus edges
# Because there are many edges, this takes a while.

Cor_sign_estimate = function(ps,edge_list){
  ps.dist = cor(otu_table(ps), use="everything", method="spearman") 
  for (i in 1:dim(edge_list)[1]){
    CorVal = ps.dist[which(rownames(ps.dist)==data.frame(edge_list)$A[i]),
                     which(colnames(ps.dist)==data.frame(edge_list)$B[i])]
    CorSign = ifelse(CorVal>0,"Positive","Negative")
    edge_list[i,3]=CorVal
    edge_list[i,4]=CorSign
  }
  edge_list$EdgeColor[edge_list$CorSign=="Positive"] = "black"
  edge_list$EdgeColor[edge_list$CorSign=="Negative"] = "red"
  
  return(edge_list)
}

edge_list = Cor_sign_estimate(ps,edge_list)

igraph = graph_from_data_frame(edge_list,directed=FALSE)

# Get the taxon info for nodes
node_list = data.frame(V(igraph)$name)

# Add the OTU ID column
colnames(node_list)[1] = "OTU"

# Add fire-responsive data to the node list
l2FC = readRDS("Results.5") # 2019
#l2FC = readRDS("Results.1") # 2015

m = l2FC%>%
  dplyr::group_by(logFC,sigBurn,OTU,Phylum,Class,Order,Family,Genus)%>%
  dplyr::summarize(MeanAbund=mean(AveExpr))%>%
  dplyr::group_by(logFC,MeanAbund,sigBurn,OTU,Phylum,Class,Order,Family,Genus)%>%
  dplyr::summarize()%>%
  dplyr::filter(!is.na(logFC))%>%
  dplyr::filter(sigBurn==1)%>%
  dplyr::filter(logFC>0)%>%
  dplyr::arrange(OTU)
m.fire.pos = levels(as.factor(m$OTU))

m = l2FC%>%
  dplyr::group_by(logFC,sigBurn,OTU,Phylum,Class,Order,Family,Genus)%>%
  dplyr::summarize(MeanAbund=mean(AveExpr))%>%
  dplyr::group_by(logFC,MeanAbund,sigBurn,OTU,Phylum,Class,Order,Family,Genus)%>%
  dplyr::summarize()%>%
  dplyr::filter(!is.na(logFC))%>%
  dplyr::filter(sigBurn==1)%>%
  dplyr::arrange(OTU)%>%
  dplyr::filter(logFC<0)
m.fire.neg = levels(as.factor(m$OTU))

# Assigning fire response variables from l2FC calculation results
Assign_fireresp = function(node_list){
  node_list$FireResponseColor = ""
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    FireResponseColor = ifelse((OTU %in% m.fire.pos),"red",
                               ifelse((OTU %in% m.fire.neg),"lightskyblue3","white"))
    node_list$FireResponseColor[i] = FireResponseColor
  }
  return(node_list)
}

node_list = Assign_fireresp(node_list)

# Adding taxonomy to the nodes list
taxonomy_labeller = function(ps,node_list){
  node_list$Species = c()
  node_list$Genus = c()
  node_list$Family = c()
  node_list$Order = c()
  node_list$Class = c()
  node_list$Phylum = c()
  node_list$Kingdom = c()
  node_list$Guild = c()
  node_list$Trophic.Mode = c()
  node_list$Growth.Morphology=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    Species = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Species"]),"")
    Genus = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Genus"]),"")
    Family = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Family"]),"")
    Order = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Order"]),"")
    Class = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Class"]),"")
    Phylum = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Phylum"]),"")
    Kingdom = ifelse(OTU %in% row.names(tax_table(ps)),paste(tax_table(ps)[paste(OTU),"Kingdom"]),"")
    #Guild = ifelse(OTU %in% row.names(tax_table(F.guilds)),paste(tax_table(F.guilds)[paste(OTU),"Guild"]),"")
    #Trophic.Mode = ifelse(OTU %in% row.names(tax_table(F.guilds)),paste(tax_table(F.guilds)[paste(OTU),"Trophic.Mode"]),"")
    #GuildGrowth.Morphology= ifelse(OTU %in% row.names(tax_table(F.guilds)),paste(tax_table(F.guilds)[paste(OTU),"Growth.Morphology"]),"")
    node_list$Species[i] = Species
    node_list$Genus[i] = Genus
    node_list$Family[i] = Family
    node_list$Order[i] = Order
    node_list$Class[i] = Class
    node_list$Phylum[i] = Phylum
    node_list$Kingdom[i] = Kingdom
    #node_list$Guild[i] = Guild
    #node_list$Trophic.Mode[i] = Trophic.Mode
    #node_list$Growth.Morphology[i] = GuildGrowth.Morphology
    node_list$SpecialName = paste("My name is",node_list$OTU)
  }
  return(node_list)
}
node_list = taxonomy_labeller(ps,node_list)


# Calculate standard properties of networks

network_properties = function(edge_list){
  N = graph_from_data_frame(edge_list,directed=FALSE)
  # Creates the network
  
  m = ecount(N)
  # Number of edges
  n = vcount(N)
  # Number of nodes
  k = 2*m/n
  # Average degree
  apl = mean_distance(N,directed=FALSE)
  # Average path length
  c = transitivity(N, type="global")
  cAve = transitivity(N, type = "average")
  # Clustering coefficient of whole graph - 
  # Transitivity measures the probability that the adjacent vertices of a vertex are connected.
  
  #cl.mean = mean(closeness(N))
  #cl.sd = sd(closeness(N))
  # closeness of graph
  ed = edge_density(N)
  # edge density of graph
  d = diameter(N)
  # diameter of graph
  
  # Turn it into a dataframe
  #df = data.frame(m,n,k,apl,c,cAve,cl.mean,cl.sd,ed,d)
  df = data.frame(m,n,k,apl,c,cAve,ed,d)
  return(df)
}

net.prop = network_properties(edge_list)
net.prop


# From https://chengjunwang.com/web_data_analysis/demo2_simulate_networks/
# plot and fit the power law distribution
fit_power_law = function(edge_list) {
  N = graph_from_data_frame(edge_list,directed=FALSE)
  # calculate degree
  d = degree(N, mode = "all")
  dd = degree.distribution(N, mode = "all", cumulative = FALSE)
  degree = 1:max(d)
  probability = dd[-1]
  # delete blank values
  nonzero.position = which(probability != 0)
  probability = probability[nonzero.position]
  degree = degree[nonzero.position]
  reg = lm(log(probability) ~ log(degree))
  cozf = coef(reg)
  power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
  alpha = -cozf[[2]]
  R.square = summary(reg)$r.squared
  print(paste("Alpha =", round(alpha, 3)))
  print(paste("R square =", round(R.square, 3)))
  # plot
  plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
       col = 1, main = "Degree Distribution")
  curve(power.law.fit, col = "red", add = T, n = length(d))
}

fit_power_law(edge_list)

# Modules can be detected using the cluster walktrap method - similar to greedy, somewhat smaller modules
# Pascal Pons, Matthieu Latapy: Computing communities in large networks using random walks
# http://arxiv.org/abs/physics/0512106

Modules = cluster_walktrap(graph=igraph, merges = TRUE, modularity = TRUE, membership = TRUE)



# Modularity (M) is an index measuring the extent to which a network is divided into modules, 
# and we used M > 0.4 as the threshold to define modular structures (Newman 2006)
# All our networks have meaningful modularity.
modularity(Modules)

# Connectivity of each node can be determined based on its within-module connectivity (Zi)
# and among-module connectivity (Pi) (Guimera & Amaral 2005)

# Zi and Pi can be used to classify the nodes based on the topological roles they play in the network
# Node topologies are organised into four categories: module hubs (highly connected nodes within modules, Zi > 2.5)
# network hubs (highly connected nodes within entire network, Zi > 2.5 and Pi > 0.62)
# connectors (nodes that connect modules, Pi > 0.62)
# and peripherals (nodes connected in modules with few outside connections, Zi < 2.5 and Pi < 0.62) (Olesen et al. 2007; Zhou et al. 2010; Deng et al. 2012).

# To calculate the Zi and Pi of each node:
# First, find out which module it is in
# Make a list of all the other nodes in that module
# Calculate the connectivity of that node to all those other nodes
# Do this for each node
# Then, Zi is calculated as:
# (number of links from a given node to other nodes in the module - the average number for nodes in this module)
# Divided by the standard deviation of this value for nodes in this module.

# Then, do the same, but make the nodes list all the nodes in other modules

# First attempt at calculating Zi

adding_Zi = function(node_list,Modules,igraph){
  
  Zi=data.frame(Name=node_list$OTU,ModuleNumber=rep(0,length(node_list$OTU)),CON=rep(0,length(node_list$OTU)))
  # Establish empty dataframe
  
  for (i in 1:length(node_list$OTU)){
    node = paste(node_list$OTU[i])
    ModuleNumber = Position(function(x) node %in% x, Modules[], nomatch = 0)
    if(ModuleNumber !=0){
      ModuleNodes = Modules[[ModuleNumber]]
      modgraph = induced_subgraph(graph=igraph, v=c(ModuleNodes,node))
      CON=try(as.numeric(lengths(adjacent_vertices(modgraph,node))),TRUE)
      if(isTRUE(class(CON)=="try-error")) { CON=NA } else {CON = as.numeric(lengths(adjacent_vertices(modgraph,node)))}
    }
    if(ModuleNumber ==0){CON=NA}
    Zi$Name[i]=as.factor(node_list$OTU[i])
    Zi$ModuleNumber[i]=ModuleNumber
    Zi$CON[i]=CON
  }
  
  # (number of links from a given node to other nodes in the module - the average number for nodes in this module)
  # Divided by the standard deviation of this value for nodes in this module.
  Zi = Zi %>%
    filter(!is.na(CON))%>%
    group_by(ModuleNumber)%>%
    mutate(MeanCON=mean(CON))%>%
    mutate(SdCON=sd(CON))%>%
    mutate(Zi=((CON-MeanCON)/SdCON))
  return(Zi)
}
Zi = adding_Zi(node_list,Modules,igraph)

# Next, we add Pi
adding_Pi = function(node_list,Modules,igraph){
  
  Pi=data.frame(Name=rep(node_list$OTU,dim(Modules[])),HomeModuleNumber=rep(0,length(node_list$OTU)),OtherModuleNumber=rep(0,length(node_list$OTU)),TotalCON=rep(0,length(node_list$OTU)),CON=rep(0,length(node_list$OTU)))
  # Establish empty dataframe
  
  for (i in 1:length(node_list$OTU)){
    node = paste(node_list$OTU[i])
    HomeModuleNumber = Position(function(x) node %in% x, Modules[], nomatch = 0)
    ModuleNumbers = 1:dim(Modules[])
    n = length(ModuleNumbers)
    TotalCON = as.numeric(lengths(adjacent_vertices(igraph,node)))
    lowend = (i-1)*n+1
    highend = n*i
    Pi$Name[lowend:highend]=node_list$OTU[i]
    Pi$HomeModuleNumber[lowend:highend]=HomeModuleNumber
    Pi$OtherModuleNumber[lowend:highend]=ModuleNumbers
    Pi$TotalCON[lowend:highend]=TotalCON
    if(HomeModuleNumber !=0){    
      for (j in ModuleNumbers){
        OtherModuleNumber = j
        NodesInOtherModule = Modules[[OtherModuleNumber]]
        modgraph = induced_subgraph(graph=igraph, v=c(node,NodesInOtherModule))
        CON=as.numeric(lengths(adjacent_vertices(modgraph,node)))
        Pi$CON[Pi$HomeModuleNumber==HomeModuleNumber & Pi$OtherModuleNumber==OtherModuleNumber & Pi$Name==node]=CON
      }
    }
  }
  Pi$kk2 = (Pi$CON/Pi$TotalCON)^2
  return(Pi)
}

Pi = adding_Pi(node_list,Modules,igraph)

Pifinal = Pi %>%
  dplyr::group_by(Name,HomeModuleNumber,TotalCON)%>%
  dplyr::summarize(Sum=sum(kk2))%>%
  dplyr::mutate(Pi=1-Sum)

# Calculating Zi from Pi data
Zinew = Pi %>% filter(HomeModuleNumber==OtherModuleNumber) %>% mutate(MeanCON=mean(CON),SdCON=sd(CON),Zi=((CON-MeanCON)/SdCON))

# Bringing module data together
# Thresholds based on paper
Making_module_data = function(Pifinal,Zinew){
  Pthresh = 0.62
  Zthresh = 2.5
  ModuleData=data.frame(Name=Pifinal$Name,Module=Pifinal$HomeModuleNumber,TotalCON=Pifinal$TotalCON,ModuleCON=Zinew$MeanCON,Pi=Pifinal$Pi,Zi=Zinew$Zi)
  ModuleData$Class = ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi>Pthresh,"Network Hub",
                            ifelse(ModuleData$Zi>Zthresh & ModuleData$Pi<Pthresh,"Module Hub",
                                   ifelse(ModuleData$Zi<Zthresh & ModuleData$Pi>Pthresh,"Connector", "Peripheral")))
  return(ModuleData)
}

ModuleData = Making_module_data(Pifinal,Zinew)

p = ggplot(ModuleData)
p = p + geom_point(aes(x=Pi,y=Zi,color=Class))
p


## Add this info to the nodes list

add_modInfo = function(node_list,ModuleData){
  node_list$Pi=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Pi!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Pi, NA))
    node_list$Pi[i] = x
  }
  
  node_list$Zi=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Zi!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Zi, NA))
    node_list$Zi[i] = x
  }
  
  node_list$NetworkRole=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Class!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Class, NA))
    node_list$NetworkRole[i] = x
  }
  
  node_list$Module=c()
  for (i in 1:dim(node_list)[1]){
    OTU = node_list$OTU[i]
    x = ifelse(OTU %in% ModuleData$Name, 
               ifelse(ModuleData[ModuleData$Name==OTU,]$Module!=(-Inf),ModuleData[ModuleData$Name==OTU,]$Module, NA))
    node_list$Module[i] = x
  }

    node_list$NetworkRoleColour = ifelse(node_list$NetworkRole=="Connector","red",
                                       ifelse(node_list$NetworkRole=="Module Hub","navy","white"))
  
  return(node_list)
}

node_list = add_modInfo(node_list,ModuleData)

dim(node_list)

# Check most abundant modules
ModProps = node_list %>% dplyr::group_by(Module) %>% dplyr::summarize(Total=n()) %>% dplyr::arrange(-Total)

head(ModProps)
dim(ModProps)
# Pull out all modules with more than (modcutoff) nodes

modcutoff = 20
AbundMods = ModProps %>% filter(Total>modcutoff) %>% select(Module)
AbundMods = AbundMods$Module

# Reporting properties of each module wrt bacteria (for modules with more than modcutoff taxa)

get_ModProps_bact = function(node_list, AbundMods){
  ModProps = node_list %>%
    dplyr::filter(Kingdom=="Bacteria")%>%
    dplyr::filter(Module %in% AbundMods)%>%
    dplyr::group_by(Module,Kingdom,Phylum)%>%
    dplyr::summarize(Sum = n())%>%
    dplyr::group_by(Module,Kingdom)%>%
    dplyr::mutate(Fract = Sum/sum(Sum))%>%
    dplyr::mutate(Tot = sum(Sum))%>%
    dplyr::arrange(Module,-Fract)%>%
    dplyr::mutate(Phylum=ifelse(Phylum=="Bacteria_unclassified(0)","Unclassified",Phylum))
  ModProps$Module = as.factor(ModProps$Module)
  ModProps$Module = droplevels((ModProps$Module))
  
  ModProps$Phylum = as.factor(ModProps$Phylum)
  ModProps$Phylum = droplevels((ModProps$Phylum))
  PhyOrder = ModProps%>%
    dplyr::group_by(Phylum)%>%
    dplyr::summarize(Total=sum(Tot))%>%
    dplyr::arrange(-Total)
  PhyOrder = PhyOrder$Phylum
  ModProps$Phylum = as.factor(ModProps$Phylum)
  ModProps$Phylum = factor(ModProps$Phylum,levels=PhyOrder)
  
  return(ModProps)
}

head(node_list)



#saveRDS(edge_list,"edge_list.2019.06-04-2020")
#saveRDS(node_list,"node_list.2019.06-04-2020")
#saveRDS(igraph,"igraph.2019.06-04-2020")

#node_list = readRDS("node_list.2019.06-04-2020")
#edge_list = readRDS("edge_list.2019.06-04-2020")
#igraph = readRDS("igraph.2019.06-04-2020")

#edge_list = readRDS("edge_list.2015.19-11-2020")
#node_list = readRDS("node_list.2015.19-11-2020")
#igraph = readRDS("igraph.2015.19-11-2020")


##### Making network figures ######

# Check out basic graph
# Creating consensus graph
set.seed(100)

# Edge colours
E(igraph)$color = c(edge_list$EdgeColor)
E(igraph)$size = 4

# Setting palette for when colouring by phylum
rbPal <- colorRampPalette(brewer.pal(11,"Spectral"))

# Assigning node (taxa) properties to the network
V(igraph)$color=node_list$FireResponseColor
#V(igraph)$color=rbPal(8)[c(5,6,7,3,4,1)][as.numeric(paste(node_list$Module))]
#V(igraph)$color=c("white","white","red","white","white","white")[as.numeric(paste(test$Module.2015))]
#V(igraph)$color=node_list$NetworkRoleColour
#V(igraph)$label=node_list$Genus

# Plot the 3D network graph
graphjs(igraph,
        width=700,height=500,vertex.shape="circle",bg="white",fg="white",vertex.size=0.5,
        edge.alpha=0.75)

# Modules figure

set.seed(9)
# 9 good for 2019
# 4 is good for 2015

# Make a list of all the individual modules as their own graph
Modules = cluster_walktrap(graph=igraph, merges = TRUE, modularity = TRUE, membership = TRUE)
igraph.modules = list()
for (i in 1:length(levels(as.factor(node_list$Module)))){
  Mod = levels(as.factor(node_list$Module))[i]
  igraph.modules[[i]] = induced_subgraph(graph=igraph, v=c(Modules[[Mod]]))
}

# Just get the coordinates of each individual module graph
graphs = igraph.modules
# Put each module's nodes into a circle
layouts = lapply(graphs,layout_in_circle)
# take our module graphs, and their layouts, and put them all in the same space
# They are placed in order of size
# This results in a long matrix with the position of each vertex
lay = merge_coords(graphs, layouts)
# In this step, the graphs are all merged
# It assumes no overlap (or relabels vertices to create this)
g = disjoint_union(graphs)

# Take our set of coordinates
layOrder = data.frame(lay)
# Then we assume they are in the same/correct order as in the graph list
layOrder$OTU = V(g)$name
# Here, we then put them in the same order as the igraph names
layOrder$OTU = factor(layOrder$OTU, levels = V(igraph)$name)
layOrder = layOrder%>%
  dplyr::arrange(OTU)
lay = as.matrix(layOrder[,1:2])

# Go back to the original network, and give it this layout
igraph.Modules = igraph
igraph.Modules$layout = lay
# This tells us where each node should be

plot.igraph(igraph.Modules, vertex.size=3, vertex.label=NA,
            edge.width=c(0.2),edge.alpha=c(0.1)
)

#### If we wanted to plot just the three most abundant modules ####

# 2,6,3 for 2015, 1,2,5 for 2019
AbundMods = c(1,2,5)
igraph.3Modules = induced_subgraph(igraph, node_list$OTU[node_list$Module %in% c(AbundMods)])
node_list.3Modules = node_list[node_list$Module %in% c(AbundMods),]

set.seed(9)
# 9 good for 2019

# Make a list of all the individual modules as their own graph

# Get the same modules as the original graph
Modules = cluster_walktrap(graph=igraph, merges = TRUE, modularity = TRUE, membership = TRUE)
# Make list for only the modules of interest in the 3-module sub-graph
igraph.modules.3Modules = list()
for (i in 1:length(levels(as.factor(node_list.3Modules$Module)))){
  Mod = levels(as.factor(node_list.3Modules$Module))[i]
  igraph.modules.3Modules[[i]] = induced_subgraph(graph=igraph.3Modules, v=c(Modules[[Mod]]))
}

# Just get the coordinates of each individual module graph
graphs = igraph.modules.3Modules

# Put each module's nodes into a circle
layouts = lapply(graphs,layout_in_circle)
# take our module graphs, and their layouts, and put them all in the same space
# They are placed in order of size
# This results in a long matrix with the position of each vertex
lay = merge_coords(graphs, layouts)
# In this step, the graphs are all merged
# It assumes no overlap (or relabels vertices to create this)
g = disjoint_union(graphs)

# Take our set of coordinates
layOrder = data.frame(lay)
# Then we assume they are in the same/correct order as in the graph list
layOrder$OTU = V(g)$name
# Here, we then put them in the same order as the igraph names
layOrder$OTU = factor(layOrder$OTU, levels = V(igraph.3Modules)$name)
layOrder = layOrder%>%
  dplyr::arrange(OTU)
lay = as.matrix(layOrder[,1:2])
# Go back to the original network, and give it this layout
igraph.3Modules$layout = lay
# This tells us where each node should be

V(igraph.3Modules)$color=node_list.3Modules$FireResponseColor
dim(node_list.3Modules)

plot.igraph(igraph.3Modules, vertex.size=3, vertex.label=NA,
            edge.width=c(0.2),edge.alpha=c(0.1)
)





# Cross-module properties
ps = ps.2019
node_list = node_list.19

mdf = psmelt(ps)

mdf2 = mdf %>%
  dplyr::filter(OTU %in% node_list$OTU)

# Get only network OTUs, change naming scheme to match node_list
M = node_list[,c("OTU","Module")]
# Get just the two relevant columns
mdf2 = merge(mdf2,M,by="OTU",all.x=TRUE)
# Merge the two dataframes


OTUbyMod = mdf2 %>%
  dplyr::group_by(Sample,Module,CBI,pH,Org_or_Min,Moisture_Regime,Veg_Comm,Land_Class,Severity_Class,TC_pct,Kingdom)%>%
  dplyr::summarize(Abundance = sum(Abundance))%>%
  dplyr::mutate(ModuleName = paste("Module",Module))%>%
  dplyr::mutate(pHGroup = ifelse(pH>7,">7",
                          ifelse(pH>5,"5-7","<5")))%>%
  dplyr::mutate(OMGroup = ifelse(TC_pct>35,">35",
                          ifelse(TC_pct>5,"5-35",
                                 ifelse(TC_pct>2,"2-5","<2"))))
OTUbyMod$pHGroup = factor(OTUbyMod$pHGroup, levels = c("<5", "5-7", ">7"))
OTUbyMod$OMGroup = factor(OTUbyMod$OMGroup, levels = c("<2", "2-5","5-35",">35"))

# 2,6,3 for 2015, 1,2,3 for 2019
plotmods = c(1,2,5) # 2019
#plotmods = c(2,6,3) # 2015

d = OTUbyMod %>%
  dplyr::filter(Module %in% plotmods)%>%
  dplyr::filter(!is.na(pHGroup))%>%
  dplyr::mutate(StripName = paste(Module))%>%
  dplyr::group_by(Sample,pHGroup,Module)%>%
  dplyr::summarize(Abundance=sum(Abundance))


p = ggplot(data=d, aes(x=pHGroup,y=Abundance,color=pHGroup))
p = p + geom_boxplot()
p = p + theme_bw()
palette = (brewer.pal(9, "YlGnBu")[3:9])[c(4,3,2)]
p = p + scale_color_manual(values=palette)
p = p + facet_wrap(~Module, scales="free_y",ncol=6) + expand_limits(y=0)
p = p + ylab("Relative Abundance") + xlab("pH range") + guides(color = "none")
p = p + theme(axis.text=element_text(size=15),
              axis.title=element_text(size=15),
              strip.text=element_text(size=15))
p

### Questions for networks comparisons ###
# Which co-occurrences persisted?
# Did the modules remain the same?
node_list.15 = readRDS("node_list.2015.06-04-2020")
edge_list.15 = readRDS("edge_list.2015.06-04-2020")

node_list.19 = readRDS("node_list.2019.06-04-2020")
edge_list.19 = readRDS("edge_list.2019.06-04-2020")

# Network dimensions
# Nodes
dim(node_list.15)
dim(node_list.19)
sum(node_list.19$OTU %in% node_list.15$OTU)
# Edges
dim(edge_list.15)
dim(edge_list.19)

# All possible matchups of OTUs
combos.15 = c(paste(edge_list.15$A,edge_list.15$B),paste(edge_list.15$B,edge_list.15$A))
combos.19 = c(paste(edge_list.19$A,edge_list.19$B))
sum(combos.19 %in% combos.15)
length(combos.19)
# Opposite direction
combos.19 = c(paste(edge_list.19$A,edge_list.19$B),paste(edge_list.19$B,edge_list.19$A))
combos.15 = c(paste(edge_list.15$A,edge_list.15$B))
sum(combos.15 %in% combos.19)
length(combos.15)

# Are the co-occurrences the same sign?
combos.15 = c(paste(edge_list.15$A,edge_list.15$B,edge_list.15$CorSign),paste(edge_list.15$B,edge_list.15$A,edge_list.15$CorSign))
combos.19 = c(paste(edge_list.19$A,edge_list.19$B,edge_list.19$CorSign))
sum(combos.19 %in% combos.15)
length(combos.19)
length(combos.15)
# Mostly they are, but a few switched direction of correlation (about 10%)

# Getting a list of those switchers (including both orders of OTUs to catch all)
combos.15 = c(paste(edge_list.15$A,edge_list.15$B),paste(edge_list.15$B,edge_list.15$A))
combos.19 = c(paste(edge_list.19$A,edge_list.19$B))
present.both=(combos.19 %in% combos.15)
combos.15 = c(paste(edge_list.15$A,edge_list.15$B,edge_list.15$CorSign),paste(edge_list.15$B,edge_list.15$A,edge_list.15$CorSign))
combos.19 = c(paste(edge_list.19$A,edge_list.19$B,edge_list.15$CorSign))
same.both = (combos.19 %in% combos.15)

combos.19 = c(paste(edge_list.19$A,edge_list.19$B,edge_list.15$CorSign))
combos.19.same.both = combos.19[present.both & same.both]
combos.19.diff.both = combos.19[present.both & !same.both]
length(combos.19.same.both) + length(combos.19.diff.both)
length(combos.19.same.both)
length(combos.19.diff.both)
# Seems to have worked - have both lists.
# Could look at things like how similar each pairing is?
# Idea would be more closely related taxa more likely to stay same correlation.
# Could also look at simple taxonomy to start.

diff.both = (strsplit(combos.19.diff.both,split=" "))
diff.both = as.data.frame(do.call(rbind, diff.both))  
colnames(diff.both) = c("A","B","CorSign19")
head(diff.both)

taxonomy = data.frame(as(tax_table(ps.2019),"matrix"))
taxonomy$A = row.names(taxonomy)
diff.both = merge(diff.both,taxonomy,by="A")
colnames(taxonomy)[9]="B"
diff.both = merge(diff.both,taxonomy,by="B",suffixes=c(".A",".B"))
head(diff.both)

sum(diff.both$Phylum.A == diff.both$Phylum.B)/dim(diff.both)[1]
# 20% of the pairs that changed sign were from the same phylum
sum(diff.both$Family.A == diff.both$Family.B)/dim(diff.both)[1]
# 6% for family level.

same.both = (strsplit(combos.19.same.both,split=" "))
same.both = as.data.frame(do.call(rbind, same.both))  
colnames(same.both) = c("A","B","CorSign19")
head(same.both)

taxonomy = data.frame(as(tax_table(ps.2019),"matrix"))
taxonomy$A = row.names(taxonomy)
same.both = merge(same.both,taxonomy,by="A")
colnames(taxonomy)[9]="B"
same.both = merge(same.both,taxonomy,by="B",suffixes=c(".A",".B"))
head(same.both)

sum(same.both$Phylum.A == same.both$Phylum.B)/dim(same.both)[1]
# 27% of the pairs that didn't change sign were from the same phylum.
# Somewhat similar to the taxa that did change sign.
sum(same.both$Family.A == same.both$Family.B)/dim(same.both)[1]
# 10% for family level.

# We can also compare how the sign shifted
sum(diff.both$CorSign19 == "Positive")
dim(diff.both)[1]
sum(diff.both$CorSign19 == "Positive")/dim(diff.both)[1]
# About 20% of the pairings went from negative in 2015 to positive in 2019 (7 pairs)
# So, conversely, 80% of the pairings went from positive to negative
sum(same.both$CorSign19 == "Positive")/dim(same.both)[1]
# Almost all of the pairings that stayed the same year to year were positive

sum(node_list.19$OTU %in% node_list.15$OTU)
length(node_list.15$OTU)

### One of my other questions is whether taxa stayed in the same modules, generally.
node_list.19.mod15 = merge(node_list.19,node_list.15[,c("OTU","Module")],by="OTU",suffixes=c(".19",".15"))
ggplot(node_list.19.mod15,aes(x=Module.15,y=Module.19))+geom_jitter(width=0.2,height=0.2)
# Basically, Module 6 from 2015 is Module 5 from 2019;
# Module 1 from 2015 is Module 1 from 2019
# Module 2 from 2015 got broken into multiple modules - 2,3,4,6 - in 2019, and some went over to Module 3.
# Some taxa not in both networks, so not all modules included.
# This is pretty neat.

# Modules X from 2015 were high pH
# Module X from 2015 was low pH

# Overarching findings is general stability of modules between years
# Issue of different number of nodes in the networks.
# Could play with cutoffs to get closer to same expected false positive level (currently the same vals)
# Could be some interesting questions about which taxa switched modules.
# E.g., which taxa were in module 1, but jumped to module 3?

Switch = node_list.19.mod15 %>%
  dplyr::filter(Module.15 ==1 & Module.19 == 3)%>%
  dplyr::group_by(Order)%>%
  dplyr::summarize(TotalOTUs = n())%>%
  dplyr::arrange(-TotalOTUs)
Switch

NoSwitch = node_list.19.mod15 %>%
  dplyr::filter(Module.15 ==1 & Module.19 == 1)%>%
  dplyr::group_by(Order)%>%
  dplyr::summarize(TotalOTUs = n())%>%
  dplyr::arrange(-TotalOTUs)
NoSwitch

# Similar order-level taxa are staying behind vs. making the switch.
# Not totally sure what would be causing the switch.

# Could look at fire response, too.
Resp.15 = readRDS("Results.1")
Resp.19 = readRDS("Results.5")
Resp.19$Response = ifelse(!is.na(Resp.19$logFC) & Resp.19$logFC>0 & Resp.19$sigBurn==1.0,"positive",
                          ifelse(!is.na(Resp.19$logFC) & Resp.19$logFC<0 & Resp.19$sigBurn==1.0,"negative","neutral"))
Resp.15$Response = ifelse(!is.na(Resp.15$logFC) & Resp.15$logFC>0 & Resp.15$sigBurn==1.0,"positive",
                          ifelse(!is.na(Resp.15$logFC) & Resp.15$logFC<0 & Resp.15$sigBurn==1.0,"negative","neutral"))

node_list.19.mod15 = merge(node_list.19.mod15,Resp.19[,c("OTU","Response")],by="OTU")
node_list.19.mod15 = merge(node_list.19.mod15,Resp.15[,c("OTU","Response")],by="OTU",suffixes=c(".19",".15"))
ggplot(node_list.19.mod15,aes(x=Module.15,y=Module.19, fill=Response.15, color = Response.19))+
  geom_jitter(width=0.4,height=0.25,shape=21,stroke=2,size=2)+
  scale_fill_manual(values=c("yellow","darkgoldenrod2","red3"))+
  scale_colour_manual(values=c("lemonchiffon1","darkgoldenrod1","red"))

# Second figure checking modules
ggplot(node_list.19.mod15,aes(x=Module.15,y=Module.19, fill=as.factor(Module.15)))+
  geom_jitter(width=0.3,height=0.3,shape=21,stroke=2,size=2,alpha=0.5)
  #scale_fill_manual(values=c("yellow","darkgoldenrod2","red3"))
  #scale_colour_manual(values=c("lemonchiffon1","darkgoldenrod1","red"))

# This figure is cool, because it indicates the following:
# 1. There are two pretty persistent positive responder modules (2/4) and (3/3)
# 2. The second positive responder module, 3/3, is joined by some negative or neutral responders
# in 2019.
# 3. These responders come from a larger negative/neutral responder module from 2015,
# which, in 2019, is split into three, with one module joining the 3/3 module
# 4. Now we might ask, the 2/4 module seems the most persistent / indicative of positive response
# I might guess this would be higher pH module, as it may indicate sites that had the greatest
# pH shift with fire.
# This could suggest that the effects of fire on the lower pH soils is less severe?
# Will be great to get the pH shift data.

# Let's look at the consistent module 2/4 taxa.
ConsistPosMod = node_list.19.mod15[node_list.19.mod15$Module.15==2 & node_list.19.mod15$Module.19==4,"OTU"]

mdf$ConsistPosMod = ifelse(mdf$OTU %in% ConsistPosMod,"Yes","No")

PosMod = mdf %>%
  dplyr::filter(ConsistPosMod=="Yes")
colnames(PosMod)
ggplot(PosMod,aes(x=pH,y=log(Abundance)))+geom_point()+facet_wrap(~OTU)
ggplot(PosMod,aes(x=Burn_Severity_Index,y=log(Abundance)))+geom_point()+facet_wrap(~OTU)

# Compare to less strong pos mod
ConsistMod3 = node_list.19.mod15[node_list.19.mod15$Module.15==3 & node_list.19.mod15$Module.19==3,"OTU"]

mdf$ConsistMod3 = ifelse(mdf$OTU %in% ConsistMod3,"Yes","No")

Mod3 = mdf %>%
  dplyr::filter(ConsistMod3=="Yes")
ggplot(Mod3,aes(x=pH,y=log(Abundance)))+geom_point()+facet_wrap(~OTU)
ggplot(Mod3,aes(x=Burn_Severity_Index,y=log(Abundance)))+geom_point()+facet_wrap(~OTU)

# I think in general, the consistent 3/3 responders are more neutral or even negative to BSI
# and tend also to non-responsive or negatively associated with high pH sites.

ggplot(Mod3,aes(x=pH,y=Burn_Severity_Index))+geom_point()
# There's a weak positive correlation between burn severity and pH (I think) [in burned sites]

# Possible hypothesis - sites that shift in pH more produce greatest community changes
# and are also then more "locked in" to the shifted community composition.

## Thought of another cool way to illustrate the module changes

# Harmonize the equivalent values with a new variable
d = node_list.19.mod15 %>% 
  #dplyr::filter(Module.15 %in% c(1,2,3,4,5,6)) %>% 
  #dplyr::filter(Module.19 %in% c(1,2,3,4,5,6))
  dplyr::filter(Module.15 %in% c(2,3,6)) %>% 
  dplyr::filter(Module.19 %in% c(1,2,5))
#d[d$Module.15==2,]$Module.15 = 5
#d[d$Module.19==4,]$Module.19 = 5

# Remapping the modules to keep the same numeric IDs
d[d$Module.15==2,]$Module.15 = 1
d[d$Module.15==3,]$Module.15 = 2
d[d$Module.15==6,]$Module.15 = 3
d[d$Module.19==2,]$Module.19 = 3
d[d$Module.19==1,]$Module.19 = 2
d[d$Module.19==5,]$Module.19 = 1


# We need to kind of expand this so the response is in the same column, corresponding to different x-vals (1 and 2)
head(d)
d = melt(d, measure.vars = c("Module.19","Module.15"),variable_name="Module")
d$Colour = ifelse(d$Module == "Module.15",d$Response.15,d$Response.19)
d$xval = ifelse(d$Module=="Module.15",1,5)
d$xval = as.numeric(d$xval)
d$value = as.numeric(d$value)
# Also get wanted colour for each point

# We need a midpoint
midpoints = d %>%
  dplyr::group_by(OTU)%>%
  dplyr::summarize(value=mean(value),Colour="Module.Mid",xval=3)
head(midpoints)

midpoints = merge(node_list.19.mod15,midpoints,by="OTU")
midpoints=midpoints%>% select(!Module.19)
midpoints=midpoints%>% select(!Module.15)
midpoints$Module = "Midpoint"
head(midpoints)
colnames(midpoints)
colnames(d)
# Rearrange column names to match for rbind
midpoints = midpoints[,c(1:16,20,17,18,19)]
colnames(midpoints) == colnames(d)
d = rbind(d,midpoints)
d = d %>%arrange(OTU,xval)

d = d %>%
  dplyr::group_by(OTU)%>%
  dplyr::mutate(NewColour=Colour[xval==5])
d[d$xval==3,]$Colour = d[d$xval==3,]$NewColour 

# First attempt plot
ggplot(d,aes(x=xval, y=value, group = OTU, colour=Response.15, linetype=Response.19))+
  geom_line(position=position_jitter(w=0, h=0.2))+
  geom_jitter(height=0.1,width=0)+
  scale_colour_manual(values=c("lemonchiffon1","darkgoldenrod1","red"))+
  scale_linetype_manual(values=c(3,2,1))+
  ylab("Module ID") + xlab("Years since Fire")

# Could probably get plot to be coloured by specifc year by plotting each change separately.

d1 = d[d$Response.15=="positive" & d$Response.19=="positive",]
d2 = d[d$Response.15=="positive" & d$Response.19=="neutral",]
d3 = d[d$Response.15=="positive" & d$Response.19=="negative",]
d4 = d[d$Response.15=="neutral" & d$Response.19=="positive",]
d5 = d[d$Response.15=="neutral" & d$Response.19=="neutral",]
d6 = d[d$Response.15=="neutral" & d$Response.19=="negative",]
d7 = d[d$Response.15=="negative" & d$Response.19=="positive",]
d8 = d[d$Response.15=="negative" & d$Response.19=="neutral",]
d9 = d[d$Response.15=="negative" & d$Response.19=="negative",]

p = ggplot() + theme_bw()
p = p + geom_line(data=d1,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d2,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d3,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d4,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d5,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d6,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d7,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d8,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + geom_line(data=d9,position=position_jitter(w=0, h=0.2),aes(x=xval, y=value, group = OTU, colour=Colour),alpha=0.3,size=1)
p = p + ylab("Module ID") + xlab("Years Since Fire")
p = p + geom_jitter(data=d[d$xval != 3,],aes(x=xval,y=value,fill=Colour),shape=21,width=0.1,height=0.2,alpha=0.75)
p = p + scale_colour_manual(values=c("lightskyblue3","lightgrey","red"))
p = p + scale_fill_manual(values=c("lightskyblue3","white","red"))
p = p + guides(fill=guide_legend(title="Fire Response"), color=guide_legend(title="Fire Response"))
p = p + scale_y_continuous(breaks=c(1,2,3))
p = p + scale_x_continuous(breaks=c(1,5))
p
# Could also add taxa to this network that were lost from the next year's network
# e.g. clusters of dots that appear to the left and to the right of the modules here
# Could also represent clusters as stacked 100% bar charts to better represent module breakdown
# Could be interesting to plot total module abundance across pH, etc.


### Quickly checking abundance of Caballeronia sordidicola across veg comms
mdf.2015 = psmelt(ps.2015)
mdf.2019 = psmelt(ps.2019)
mdf = psmelt(ps.merged.norm)
head(d)
d = mdf %>%
  filter(OTU == "AAGCGTTAATCGGAATTACTGGGCGTAAAGCGTGCGCAGGCGGTCCGCTAAGACAGATGTGAAATCCCCGGGCTTAACCTGGGAACTGCATTTGTGACTGGCGGGCTAGAGTATGGCAGAGGGGGGTAGAATTCCACGTGTAGCAGTGAAATGCGTAGAGATGTGGAGGAATACCGATGGCGAAGGCAGCCCCCTGGGCCAATACTGACGCTCATGCACGAAAGCGTGG")
colnames(d)
ggplot(d,aes(x=Severity_Class,y=log(Abundance), color=Years_Since_Fire))+geom_boxplot()+facet_wrap(~Veg_Comm)
