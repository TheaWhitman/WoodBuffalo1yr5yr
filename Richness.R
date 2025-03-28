### This is the R code associated with figures and analysis related to richness from the paper
### Resilience not yet apparent in soil fungal communities of the 
### boreal forest from one to five years after wildfire across a severity gradient###

# Paper authors: Thea Whitman, Jamie Woolet, Miranda C. Sikora, Dana B. Johnson, Denyse A. Dawe, and Ellen Whitman 

# © 2025. This work is openly licensed via CC BY NC SA 4.0.
# https://creativecommons.org/licenses/by-nc-sa/4.0/


library("breakaway")
library("phyloseq")
library("plyr")
library("dplyr")
library("ggplot2")
library("wesanderson")

# Import ps object
ps = readRDS("ps.merged.aligned.0pct.glommed.spp.wNAs.noITSx.cutadapt")

# Generates the frequency table summary needed by breakaway
FreqTableGenerator = function(Sample){
  df = data.frame((otu_table(ps))[,Sample])
  # Grab the OTU table
  colnames(df)="Frequency"
  # Add a column for the frequencies of the OTUs
  df = df %>%
    group_by(Frequency)%>%
    summarize(nOTUs=n())%>%
    arrange(Frequency)
  # Summarize the total OTUs that are present at each frequency
  df = df[df$Frequency>1,]
  # Cut out the 0 and 1 frequencies (because dada2 trimmed singletons)
  colnames(df)=NULL
  # Omit column names
  df = as.matrix(df)
  # Spit out a nice matrix
  df
}

Samples = sample_names(ps)

report = data.frame(name="",Richness_estimate="",Richness_stderr="",Richness_model="")

breakawayrunner = function(SampleName){
  df = FreqTableGenerator(SampleName)
  if(df[1,1]==2 & (sum(df[1:6,1])==sum(c(2:7)))){
    # First, check that the first frequency count is, indeed, 2,
    # and there are at least 6 consecutive counts (could actually just do this test)
    m = breakaway_nof1(df, answers=TRUE, plot=FALSE, print=FALSE)
    # Run breakaway for the no singletons data
    Richness_model = m$name
    Richness_estimate=m$est
    Richness_stderr = m$seest
    name = SampleName
    # Grab the outputs
    report = data.frame(name,Richness_estimate,Richness_stderr,Richness_model)
    # Generate a report
  } else {
    name = SampleName
    report=data.frame(name)
    # If breakaway isn't going to work anyway, just spit out the name (and NAs)
  }
  
  report
}

Reports = plyr::mdply(Samples,breakawayrunner)

Report_Summary = Reports %>%
  group_by(Richness_model)%>%
  summarize(n())
Report_Summary

# Some were estimated with Poisson, some with breakaway

sample_data(ps)[,colnames(Reports)[3:5]]=Reports[,3:5]

mdf = sample_data(ps)

# Plotting
p = ggplot(mdf, aes(Severity_Class, Richness_estimate, colour=Years_Since_Fire))
p = p + geom_boxplot()
#p = p + geom_point(size=3,pch=21, alpha=0.75)
#p = p + geom_errorbar(aes(ymin = Richness_estimate-1.96*Richness_stderr, ymax = Richness_estimate+1.96*Richness_stderr), width = 1)
#p = p + ylim(0,1000)
p = p + facet_wrap(~Veg_Comm~Org_or_Min)
p = p + theme(axis.text.x=element_text(angle=45, hjust=1)
              ,axis.ticks.x=element_blank()
              #,panel.grid=element_blank()
              ,panel.background=element_blank()
)
p


# Running by treatment
BettaRunner = function(Years_Since_Fire,model,OM){
  df=data.frame(mdf)
  df = df %>% filter(Richness_model==model) %>% filter(Org_or_Min == OM)
  d = df[df$Years_Since_Fire==Years_Since_Fire,]
  CoVars = data.frame(Years_Since_Fire = d$Years_Since_Fire, Severity_Class=d$Severity_Class)
  CoVars$Severity = c(rep("",dim(CoVars)[1]))
  SeverityOptions = levels(CoVars$Severity_Class)
  N = length(SeverityOptions)
  for (i in 1:dim(CoVars)[1]){
    Severity_Class = paste(CoVars$Severity_Class[i])    
    if("Unburned" %in% SeverityOptions){CoVars$Unburned[i] = ifelse(Severity_Class=="Unburned",1,0)}
    if("Low" %in% SeverityOptions){CoVars$Low[i] = ifelse(Severity_Class=="Low",1,0)}
    if("Moderate" %in% SeverityOptions){CoVars$Moderate[i] = ifelse(Severity_Class=="Moderate",1,0)}
    if("High" %in% SeverityOptions){CoVars$High[i] = ifelse(Severity_Class=="High",1,0)}
  }
  nmax = 3+N
  CoVars = as.matrix(CoVars[,4:nmax])
  Richness_estimates = as.vector(d$Richness_estimate)
  Richness_stderr = as.vector(d$Richness_stderr)
  BettaEst = betta(Richness_estimates,Richness_stderr,CoVars)
  X=data.frame(BettaEst$table)
  X$Years_Since_Fire = Years_Since_Fire
  X$Org_or_Min = OM
  X$Model = model
  X$Severity_Class = row.names(X)
  return(X)
}

x = rbind(BettaRunner(Years_Since_Fire="1",model="PoissonModel",OM="O"),BettaRunner(Years_Since_Fire="1",model="PoissonModel",OM="M"),
          BettaRunner(Years_Since_Fire="1",model="breakaway_nof1",OM="O"),BettaRunner(Years_Since_Fire="1",model="breakaway_nof1",OM="M"),
          BettaRunner(Years_Since_Fire="5",model="PoissonModel",OM="O"), BettaRunner(Years_Since_Fire="5",model="PoissonModel",OM="M"),
          BettaRunner(Years_Since_Fire="5",model="breakaway_nof1",OM="O"), BettaRunner(Years_Since_Fire="5",model="breakaway_nof1",OM="M"))

x$p.adj = p.adjust(x$p.values,method="bonferroni")
x$ymin = x$Estimates-1.96*x$Standard.Errors
x$ymax = x$Estimates+1.96*x$Standard.Errors
x$comb = paste(x$Years_Since_Fire,x$Severity_Class)
x$Severity_Class = factor(x$Severity_Class,levels=c("Unburned","Low","Moderate","High"))
x

# Plot
p = ggplot(x)
p = p + geom_hline(yintercept=0,colour="grey")
p = p + geom_errorbar(aes(x=Severity_Class,min = ymin, max = ymax,color=Years_Since_Fire),width=0.5)
p = p + geom_point(aes(x=Severity_Class,y=Estimates,color=Years_Since_Fire,shape=Years_Since_Fire),size=3)
palette = wes_palette("Darjeeling1")[c(2,3,4,1,5)]
p = p + scale_colour_manual(values=palette)
p = p + facet_grid(~Model~Org_or_Min)
p = p + theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,strip.background = element_rect(colour="white", fill="white"))
p = p + theme(axis.text.x = element_text(angle=55,hjust=1))
p = p + ylab("Mean Richness Estimate")
p = p + xlab("Burn Severity Class")
p = p + ylim(-175,550)
p = p + geom_hline(yintercept=0,colour="grey")
p = p + guides(colour=guide_legend(title="Years Since Fire"),shape=guide_legend(title="Years Since Fire"))
p

