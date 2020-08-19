#analyze the probabiliy of being killed for each plant growth form versus controls
#I first need to create a new dataframe with a separate row for each Epi_fernepiphyte

#coxph runs a mixed effect survival analysis that can handle linear predictors
#if I run a mixed effects cox regression, then I might be interested in the "strata" term
####"strata" calculates a separate baseline hazard for each group
#I am working with this example for the coxph models:https://stats.idre.ucla.edu/r/dae/mixed-effects-cox-regression/
####a better alternative might be: https://cran.r-project.org/web/packages/survival/vignettes/survival.pdf



#open relevant packages
library(dplyr)
library(lme4)
library(visreg)
library(rcompanion)
library(ggplot2)
library(DescTools)
#we are plotting, so I'll load the fonts as well
library(extrafont)
loadfonts(device="win")

#remove old data
rm(list=ls())

#import data
alldata_allyears <- read.csv("C:/Users/Evan Gora/Google Drive/Lightning/lightning_ecology_paper/data/cleaned_2015-2019_lightningcensuses.csv")
str(alldata_allyears)

##################################################
#we need to create an appropriate dataframe

###################################################
#first, select only the paired strikes and controls

#drop the 2019 controls
alldata_no2019controls<-alldata_allyears[alldata_allyears$strike!="19003N",]
alldata_no2019controls<-alldata_no2019controls[alldata_no2019controls$strike!="19004N",]
#confirm row count makes sense
nrow(alldata_no2019controls)-nrow(alldata_allyears)#looks good

#first subset the data to only include the paired strikes from 2018
control_levels<-as.character(unique(alldata_no2019controls$strike[grepl("N",alldata_no2019controls$strike)]))
control_levels#this now only has the levels of the control strikes
paired_strike_levels <-gsub("N","",control_levels)
paired_strike_levels

#use old code to grab all the controls and strikes
strike_only<-alldata_allyears[alldata_allyears$strike==18001|alldata_allyears$strike==18002|alldata_allyears$strike==18003|alldata_allyears$strike==18004|alldata_allyears$strike==18005|alldata_allyears$strike==18006|alldata_allyears$strike==18007|alldata_allyears$strike==18008|alldata_allyears$strike==18009|alldata_allyears$strike==18010|alldata_allyears$strike==18011|alldata_allyears$strike==18015,]
str(strike_only)
as.character(unique(strike_only$strike))==paired_strike_levels#we have all of our strikes

#use grepl to select all strike names that include "N"
control_only<-alldata_no2019controls[grepl("N",alldata_no2019controls$strike),]
str(control_only)
as.character(unique(control_only$strike))==control_levels#we have all of our controls

#combine the two datasets
paired_strikes<-bind_rows(strike_only,control_only)
str(paired_strikes)

#######################################################
#create a column that references the time (in days) of the previous census
paired_strikes<-paired_strikes %>%
  group_by(tree_strike) %>%
  arrange(days.post.strike) %>%
  mutate(prior_census_days = lag(days.post.strike,default = NA))

#View(paired_strikes[,c("prior_census_days","days.post.strike")])
#this works great - it might be better to employ this after re-organizing the data

#######################################################
#now create a dataframe with a separate row for each value of epiphyte, liana, and vine

#create seed df to set structure for the for loop below
seeddf<-data.frame("A","A",1,1,1,1,1,1,"A","A",1)
colnames(seeddf)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")


for(i in 1:length(paired_strikes$tree_strike)){
  #create df from healthy Epi_fernepiphtyes for that row
  #first create a value that equals
  healthy_Epi_fern_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_fern_healthy[i])))
  healthy_Epi_fern_response<-rep("healthy",paired_strikes$Epi_fern_healthy[i])
  healthy_Epi_fern_ID<-if(paired_strikes$Epi_fern_healthy[i]==0) {rep("healthy",paired_strikes$Epi_fern_healthy[i])} else {paste("fern",seq(1,paired_strikes$Epi_fern_healthy[i],1),sep="")}
  
  healthy_df<-data.frame(as.character(healthy_Epi_fern_strike),as.character(healthy_Epi_fern_tree),healthy_Epi_fern_census,healthy_Epi_fern_distance,healthy_Epi_fern_DBH,healthy_Epi_fern_days.post,healthy_Epi_fern_days.prior,healthy_Epi_fern_Dieback,as.character(healthy_Epi_fern_treatment),healthy_Epi_fern_response,healthy_Epi_fern_ID)
  colnames(healthy_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  damaged_Epi_fern_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_fern_damaged[i])))
  damaged_Epi_fern_response<-rep("damaged",paired_strikes$Epi_fern_damaged[i])
  damaged_Epi_fern_ID<-if(paired_strikes$Epi_fern_damaged[i]==0) {rep("damaged",paired_strikes$Epi_fern_damaged[i])} else {paste("fern",seq(1,paired_strikes$Epi_fern_damaged[i],1),sep="")}
  
  damaged_df<-data.frame(as.character(damaged_Epi_fern_strike),as.character(damaged_Epi_fern_tree),damaged_Epi_fern_census,damaged_Epi_fern_distance,damaged_Epi_fern_DBH,damaged_Epi_fern_days.post,damaged_Epi_fern_days.prior,damaged_Epi_fern_Dieback,as.character(damaged_Epi_fern_treatment),damaged_Epi_fern_response,damaged_Epi_fern_ID)
  colnames(damaged_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  dead_Epi_fern_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_fern_dead[i])))
  dead_Epi_fern_response<-rep("dead",paired_strikes$Epi_fern_dead[i])
  dead_Epi_fern_ID<-if(paired_strikes$Epi_fern_dead[i]==0) {rep("dead",paired_strikes$Epi_fern_dead[i])} else {paste("fern",seq(1,paired_strikes$Epi_fern_dead[i],1),sep="")}
  
  dead_df<-data.frame(as.character(dead_Epi_fern_strike),as.character(dead_Epi_fern_tree),dead_Epi_fern_census,dead_Epi_fern_distance,dead_Epi_fern_DBH,dead_Epi_fern_days.post,dead_Epi_fern_days.prior,dead_Epi_fern_Dieback,as.character(dead_Epi_fern_treatment),dead_Epi_fern_response,dead_Epi_fern_ID)
  colnames(dead_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  seeddf<-rbind.data.frame(seeddf,healthy_df,damaged_df,dead_df)
  
}
head(seeddf)
fern_data<-seeddf[-1,]
str(fern_data)
table(fern_data$response)

#confirm that all rows have data
summary(is.na(fern_data))
#all rows have data

#double check that there are the same total number of epiphytes
sum(paired_strikes$Epi_fern_healthy)+sum(paired_strikes$Epi_fern_dead)+sum(paired_strikes$Epi_fern_damaged)

#for the damage "survival" analysis - convert all "healthy" to 0 and all "damaged"/"dead" to 1

#for the dead survival analysis - convert all "healthy" and "damaged" to 0 and "dead" to 1



#######################################################
#now create a dataframe with a separate row for each value of epiphyte, liana, and vine

#create seed df to set structure for the for loop below
seeddf<-data.frame("A","A",1,1,1,1,1,1,"A","A",1)
colnames(seeddf)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")


for(i in 1:length(paired_strikes$tree_strike)){
  #create df from healthy Epi_othrepiphtyes for that row
  #first create a value that equals
  healthy_Epi_othr_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_othr_healthy[i])))
  healthy_Epi_othr_response<-rep("healthy",paired_strikes$Epi_othr_healthy[i])
  healthy_Epi_othr_ID<-if(paired_strikes$Epi_othr_healthy[i]==0) {rep("healthy",paired_strikes$Epi_othr_healthy[i])} else {paste("othr",seq(1,paired_strikes$Epi_othr_healthy[i],1),sep="")}
  
  healthy_df<-data.frame(as.character(healthy_Epi_othr_strike),as.character(healthy_Epi_othr_tree),healthy_Epi_othr_census,healthy_Epi_othr_distance,healthy_Epi_othr_DBH,healthy_Epi_othr_days.post,healthy_Epi_othr_days.prior,healthy_Epi_othr_Dieback,as.character(healthy_Epi_othr_treatment),healthy_Epi_othr_response,healthy_Epi_othr_ID)
  colnames(healthy_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  damaged_Epi_othr_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_othr_damaged[i])))
  damaged_Epi_othr_response<-rep("damaged",paired_strikes$Epi_othr_damaged[i])
  damaged_Epi_othr_ID<-if(paired_strikes$Epi_othr_damaged[i]==0) {rep("damaged",paired_strikes$Epi_othr_damaged[i])} else {paste("othr",seq(1,paired_strikes$Epi_othr_damaged[i],1),sep = "")}
  
  damaged_df<-data.frame(as.character(damaged_Epi_othr_strike),as.character(damaged_Epi_othr_tree),damaged_Epi_othr_census,damaged_Epi_othr_distance,damaged_Epi_othr_DBH,damaged_Epi_othr_days.post,damaged_Epi_othr_days.prior,damaged_Epi_othr_Dieback,as.character(damaged_Epi_othr_treatment),damaged_Epi_othr_response,damaged_Epi_othr_ID)
  colnames(damaged_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  dead_Epi_othr_strike<-rep(paired_strikes$strike[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_tree<-rep(paired_strikes$tree_strike[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_census<-rep(paired_strikes$census[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_distance<-rep(paired_strikes$distance.from.center[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_DBH<-rep(paired_strikes$DBH.mm[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_days.post<-rep(paired_strikes$days.post.strike[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_days.prior<-rep(paired_strikes$prior_census_days[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_Dieback<-rep(paired_strikes$Dieback[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_treatment<-rep(paired_strikes$control_strike[i],with(paired_strikes, sum(Epi_othr_dead[i])))
  dead_Epi_othr_response<-rep("dead",paired_strikes$Epi_othr_dead[i])
  dead_Epi_othr_ID<-if(paired_strikes$Epi_othr_dead[i]==0) {rep("dead",paired_strikes$Epi_othr_dead[i])} else {paste("othr",seq(1,paired_strikes$Epi_othr_dead[i],1),sep="")}
  
  dead_df<-data.frame(as.character(dead_Epi_othr_strike),as.character(dead_Epi_othr_tree),dead_Epi_othr_census,dead_Epi_othr_distance,dead_Epi_othr_DBH,dead_Epi_othr_days.post,dead_Epi_othr_days.prior,dead_Epi_othr_Dieback,as.character(dead_Epi_othr_treatment),dead_Epi_othr_response,dead_Epi_othr_ID)
  colnames(dead_df)<-c("strike","tree_strike","census","distance.from.center","DBH.mm","days.post.strike","prior_census_days","Dieback","control_strike","response","epi_ID")
  
  seeddf<-rbind.data.frame(seeddf,healthy_df,damaged_df,dead_df)
  
}
head(seeddf)
othr_data<-seeddf[-1,]
str(othr_data)
table(othr_data$response)

#confirm that all rows have data
summary(is.na(othr_data))
#all rows have data

#double check that there are the same total number of epiphytes
sum(paired_strikes$Epi_othr_healthy)+sum(paired_strikes$Epi_othr_dead)+sum(paired_strikes$Epi_othr_damaged)
#the counts match up

################################################
#combine the dataframes for analysis of epiphyte damage and death patterns

#combine the other epiphyte and fern epiphyte data
epiphyte_data<-bind_rows(fern_data,othr_data)
str(epiphyte_data)

###################################################
#create vectors with response variables for damaged and dead epiphytes
#for the damage "survival" analysis - convert all "healthy" to 0 and all "damaged"/"dead" to 1
#for the dead survival analysis - convert all "healthy" and "damaged" to 0 and "dead" to 1

#create damaged epiphyte binary response variable
inter_dmg1<-gsub("healthy",0,epiphyte_data$response)
inter_dmg2<-gsub("damaged",1,inter_dmg1)
epiphyte_data$epi_dmg_binary<-as.numeric(gsub("dead",1,inter_dmg2))

#create dead epiphyte binary response variable
inter_dead1<-gsub("healthy",0,epiphyte_data$response)
inter_dead2<-gsub("damaged",0,inter_dead1)
epiphyte_data$epi_dead_binary<-as.numeric(gsub("dead",1,inter_dead2))

##################################################
#look at counts of lianas - report counts from the first census as the sample size
epiphyte_data %>%
  group_by(control_strike,census) %>%
  summarise(total_lianas=length(response))

###################################################
#try out putting together a logistic regression
#####There are major issues with singularity here

#I can use different modeling approaches to handle singularity
#####package MCMCglmm might be the right approach to handle singularity (this is full bayesian)
#####package blme can use partially bayesian methods to do something similar
library(blme)

#create strike factor that is shared between strike and control trees
epiphyte_data$paired_factor<-as.numeric(gsub("N","",epiphyte_data$strike))

#re-scale continous variables
epiphyte_data$scale.distance.from.center<-scale(epiphyte_data$distance.from.center, center = FALSE)
epiphyte_data$scale.days.post.strike<-scale(epiphyte_data$days.post.strike, center = FALSE)

####################################################

#save a dataset of liana death and damage data
#"strike","census","paired_factor","distance.from.center","days.post.strike","epi_ID","liana_dmg_binary","liana_dead_binary","control_strike"
epiphyte_df<-epiphyte_data[,c(1:4,6,9,11:13)]

#convert the "treetag" name to a more general ID name and convert liana dead and dmg binary to just "dead" and "dmg
colnames(epiphyte_df)[7:9]<-c("plantID","damaged","dead")

#create plant type column
epiphyte_df$plant_type <-rep("Epiphytes",nrow(epiphyte_df))

#save df for plotting or comparative analyses
#write.csv(epiphyte_df,"epiphyte_df.csv")

####################################################
#what if we use blme (function glmer) to run logistic regression? This handles complete separation
######this sets a weak prior on the variance to push it away from complete separation for trees with only living or dead epiphytes
glm2_dead_mod1<-glmer(epi_dead_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center +
                        control_strike:scale.days.post.strike + control_strike:scale.distance.from.center + 
                        (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod1)#check convergence issues

#check for singularity
tt <- getME(glm2_dead_mod1,"theta")
ll <- getME(glm2_dead_mod1,"lower")
min(tt[ll==0])#not an issue

#compute gradient and Hessian with more intensity
library(numDeriv)
derivs1 <- glm2_dead_mod1@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))#basically the same - still very close to threshold

#redo the calculations with numDeriv
dd <- update(glm2_dead_mod1,devFunOnly=TRUE)
pars <- unlist(getME(glm2_dead_mod1,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))#even closer to the threshold - try more iterations

#try restarting model and try it with more iterations
ss <- getME(glm2_dead_mod1,c("theta","fixef"))
m2 <- update(glm2_dead_mod1,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))
#no problems now

glm2_dead_mod2<-glmer(epi_dead_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center +
                        control_strike:scale.distance.from.center + 
                        (1|paired_factor), 
                        data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod2)
anova(m2,glm2_dead_mod2)#no effect of the control/strike:time interaction

glm2_dead_mod3<-glmer(epi_dead_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center + 
                        (1|paired_factor), 
                        data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod3)
anova(glm2_dead_mod3,glm2_dead_mod2)#no effect of control/strike:distance interaction
####THIS IS THE BEST FIT MODEL

glm2_dead_mod4<-glmer(epi_dead_binary ~  scale.days.post.strike + scale.distance.from.center + 
                          (1|paired_factor), 
                        data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod4)
anova(glm2_dead_mod4,glm2_dead_mod3)#strong effect of control/strike

glm2_dead_mod5<-glmer(epi_dead_binary ~ control_strike + scale.days.post.strike + 
                        (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod5)
anova(glm2_dead_mod5,glm2_dead_mod3)#strong effect of distance from center

glm2_dead_mod6<-glmer(epi_dead_binary ~ control_strike + scale.distance.from.center + 
                        (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod6)
anova(glm2_dead_mod6,glm2_dead_mod3)#strong effect of time as well

#look at model 
visreg(glm2_dead_mod3, xvar = "scale.days.post.strike", by = "control_strike",overlay = T)
visreg(glm2_dead_mod3, xvar = "scale.distance.from.center", by = "control_strike",overlay = T)
#strong effect exactly as we would expect

##############################################################
#contrast model coefficients with the bglmer approach
bglm2_dead_mod3<-bglmer(epi_dead_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center + 
                          (1|paired_factor), 
                        data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dead_mod3)
summary(bglm2_dead_mod3)
#the fixed effects coefficients and standard errors are very similar between the models
#the random effects do differ a little bit, but not by too much


##############################################################
#plot expected probabilities in quantiles versus real proportions
#create a dataframe with the information needed for plotting this relationship only

#predict log-odds and probability for both good models
dead_logodds<-predict(glm2_dead_mod3)
dead_probs<-exp(dead_logodds)/(1+exp(dead_logodds))

#create separate vectors for control and strike data
control_probs<-dead_probs[epiphyte_data$control_strike=="Control"]
strike_probs<-dead_probs[epiphyte_data$control_strike=="Strike"]

#now combine them all into 1 dataframe
plot_quantiles_control<-data.frame("control_strike" = epiphyte_data$control_strike[epiphyte_data$control_strike=="Control"],
                                 "dead.binary" = epiphyte_data$epi_dead_binary[epiphyte_data$control_strike=="Control"],control_probs)
plot_quantiles_strike<-data.frame("control_strike" = epiphyte_data$control_strike[epiphyte_data$control_strike=="Strike"],
                                   "dead.binary" = epiphyte_data$epi_dead_binary[epiphyte_data$control_strike=="Strike"],strike_probs)

#use "cut()" to create 20 breaks based on quantiles of each model (2 different vectors - one for each control and strike)
control_quantiles<-quantile(control_probs, probs = c(0,seq(.1,.9,.1),1))
strike_quantiles<-quantile(strike_probs, probs = c(0,seq(.1,.9,.1),1))

bin_labels <- seq(1,10,1)

#look at the distributions of the quantiles
plot(c(bin_labels,11),control_quantiles)
plot(c(bin_labels,11),strike_quantiles)

#cut the data by quantile
plot_quantiles_control$cut<-cut(control_probs,control_quantiles,bin_labels)
plot_quantiles_strike$cut<-cut(strike_probs,strike_quantiles,bin_labels)

#confirm that the data are pretty evenly distributed among bins
table(plot_quantiles_control$cut)
table(plot_quantiles_strike$cut)

########################################
#now use those quantiles to calculate proporitons of real data with confidence intervals
#############################################
#now create confidence intervals based on the binomial distribution

#aggregate the number of damaged and killed trees by 5m distance class
control_dead<-aggregate(dead.binary~cut,data=plot_quantiles_control,FUN="sum")
control_total<-aggregate(dead.binary~cut,data=plot_quantiles_control,FUN="length")

#create a complete dataframe from all three variables
control_dead_prop<-data.frame(control_total,"Dead"=control_dead$dead.binary)
colnames(control_dead_prop)<-c("Quantile","Total","Dead")
str(control_dead_prop)

#use DescTools to calculate the binomial CI
CI_dataframe_control_dead<-BinomCI(control_dead_prop$Dead,control_dead_prop$Total,
                               conf.level=0.95, method = "clopper-pearson")
plottingCIs_control_dead<-data.frame(as.data.frame(CI_dataframe_control_dead),"plotting.numeric" = seq(.025,.975,.05),
                                    "plotting.label" = bin_labels, "data.type" = rep("Real proportions",10), "control_strike" = rep("Control",10))
str(plottingCIs_control_dead)

#aggregate the number of damaged and killed trees by 5m distance class
strike_dead<-aggregate(dead.binary~cut,data=plot_quantiles_strike,FUN="sum")
strike_total<-aggregate(dead.binary~cut,data=plot_quantiles_strike,FUN="length")

#create a complete dataframe from all three variables
strike_dead_prop<-data.frame(strike_total,"Dead"=strike_dead$dead.binary)
colnames(strike_dead_prop)<-c("Quantile","Total","Dead")
str(strike_dead_prop)

#use DescTools to calculate the binomial CI
CI_dataframe_strike_dead<-BinomCI(strike_dead_prop$Dead,strike_dead_prop$Total,
                                   conf.level=0.95, method = "clopper-pearson")
plottingCIs_strike_dead<-data.frame(as.data.frame(CI_dataframe_strike_dead),"plotting.numeric" = seq(.025,.975,.05),
                                    "plotting.label" = bin_labels, "data.type" = rep("Real proportions",10), "control_strike" = rep("Strike",10))
str(plottingCIs_strike_dead)


###########################################
#calculate confidence intervals of probability with "groupwiseMean"
control_CIsinter<-groupwiseMean(control_probs~cut, data = plot_quantiles_control, conf = .95)
control_CIs<-data.frame(control_CIsinter[-11,c(1,3,5,6)],data.type=rep("Predicted probabilities",10), "control_strike" = rep("Control",10))

strike_CIsinter<-groupwiseMean(strike_probs~cut, data = plot_quantiles_strike, conf = .95)
strike_CIs<-data.frame(strike_CIsinter[-11,c(1,3,5,6)],data.type=rep("Predicted probabilities",10), "control_strike" = rep("Strike",10))

#################################################


##########################
#now create 2 plotting dataframe
#rename the columns to facilitate combination
colnames(control_CIs)[c(1,3,4)]<-c("plotting.label","upr.ci","lwr.ci")
colnames(strike_CIs)[c(1,3,4)]<-c("plotting.label","upr.ci","lwr.ci")
colnames(plottingCIs_control_dead)[1]<-"Mean"
colnames(plottingCIs_strike_dead)[1]<-"Mean"

#add plotting label variable and change label factor to numeric
control_CIs$plotting.label<-as.numeric(as.character(control_CIs$plotting.label))
strike_CIs$plotting.label<-as.numeric(as.character(strike_CIs$plotting.label))
control_CIs$plotting.numeric<-plottingCIs_control_dead$plotting.numeric
strike_CIs$plotting.numeric<-plottingCIs_strike_dead$plotting.numeric

#bind rows to combine data frames
control_plotting<-bind_rows(plottingCIs_control_dead,control_CIs)
str(control_plotting)
strike_plotting<-bind_rows(plottingCIs_strike_dead,strike_CIs)
str(strike_plotting)

#create a combined plotting variable
combined_plotting<-bind_rows(strike_plotting,control_plotting)

################################
#create ggplot theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=16))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 14))

#now plot data
epideath_glmfit<-ggplot(data = combined_plotting, aes(y = Mean, x = plotting.label, shape = data.type,
                                                        color = control_strike))+
  geom_errorbar(aes(ymax = upr.ci, ymin = lwr.ci),width = .2, position = position_dodge(width = .6))+
  geom_point(size = 2, position = position_dodge(width = .6))+facet_wrap(~control_strike, nrow = 2,scales = "free_x")+
  scale_shape_manual(values = c(1,16)) + theme_basis +
  theme(strip.text.x = element_text(family = "Arial",colour="black", face ="bold",size=14),
        strip.background = element_blank())+
  guides(color = FALSE)+
  theme(legend.title=element_blank(),legend.background = element_blank(),
        legend.text=element_text(family = "Arial",colour="black", face ="bold",size=12), 
        legend.position=c(.7,.9),   legend.spacing.y = unit(.001, 'cm'),
        legend.key = element_blank(),legend.key.height=unit(.9,"line"))+
  scale_x_continuous(name = "Quantiles of Epiphyte Death", breaks = seq(0,10,1),labels = c("0.0","","0.2","","0.4","","0.6","","0.8","","1.0"))+
  scale_y_continuous(name = "Probability of Death",breaks = seq(0,.6,.1),labels = c("0.0","","0.2","","0.4","","0.6"))

epideath_glmfit#this looks pretty good (not perfect, but good enough)

ggsave("epideath_glmfit.tiff",epideath_glmfit,dpi = 600,width = 3.3, height = 3, scale = 1.6, compression = "lzw")
getwd()


####################################################
#what if we use blme (function glmer) to run logistic regression? This handles complete separation
######this sets a weak prior on the variance to push it away from complete separation for trees with only living or dmg epiphytes
glm2_dmg_mod1<-glmer(epi_dmg_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center +
                        control_strike:scale.days.post.strike + control_strike:scale.distance.from.center + 
                        (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dmg_mod1)

glm2_dmg_mod2<-glmer(epi_dmg_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center +
                        control_strike:scale.distance.from.center + 
                        (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dmg_mod2)
anova(glm2_dmg_mod1,glm2_dmg_mod2)#no interaction between control/strike:time interaction

glm2_dmg_mod3<-glmer(epi_dmg_binary ~ control_strike + scale.days.post.strike + scale.distance.from.center + 
                      (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dmg_mod3)
anova(glm2_dmg_mod1,glm2_dmg_mod3)#significant effect of control/strike:distance interaction

glm2_dmg_mod4<-glmer(epi_dmg_binary ~  scale.days.post.strike + scale.distance.from.center + 
                       control_strike:scale.distance.from.center +  (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dmg_mod4)
anova(glm2_dmg_mod4,glm2_dmg_mod2)#strong effect of control/strike

glm2_dmg_mod6<-glmer(epi_dmg_binary ~ control_strike + scale.distance.from.center + 
                       control_strike:scale.distance.from.center + (1|paired_factor), 
                      data = epiphyte_data, family = binomial(link = "logit"))
summary(glm2_dmg_mod6)
anova(glm2_dmg_mod6,glm2_dmg_mod2)#non-significant effect of time as well
#BEST FIT MODEL

#look at model 
visreg(glm2_dmg_mod6, xvar = "scale.distance.from.center", by = "control_strike",overlay = T)
#strong effect exactly as we would expect

######################################################################
##############################################################
#plot expected probabilities in quantiles versus real proportions
#create a dataframe with the information needed for plotting this relationship only

#predict log-odds and probability for both good models
dmg_logodds<-predict(glm2_dmg_mod6)
dmg_probs<-exp(dmg_logodds)/(1+exp(dmg_logodds))

#create separate vectors for control and strike data
control_probs<-dmg_probs[epiphyte_data$control_strike=="Control"]
strike_probs<-dmg_probs[epiphyte_data$control_strike=="Strike"]

#now combine them all into 1 dataframe
plot_quantiles_control<-data.frame("control_strike" = epiphyte_data$control_strike[epiphyte_data$control_strike=="Control"],
                                   "dmg.binary" = epiphyte_data$epi_dmg_binary[epiphyte_data$control_strike=="Control"],control_probs)
plot_quantiles_strike<-data.frame("control_strike" = epiphyte_data$control_strike[epiphyte_data$control_strike=="Strike"],
                                  "dmg.binary" = epiphyte_data$epi_dmg_binary[epiphyte_data$control_strike=="Strike"],strike_probs)

#use "cut()" to create 20 breaks based on quantiles of each model (2 different vectors - one for each control and strike)
control_quantiles<-quantile(control_probs, probs = c(0,seq(.1,.9,.1),1))
strike_quantiles<-quantile(strike_probs, probs = c(0,seq(.1,.9,.1),1))

bin_labels <- seq(1,10,1)

#look at the distributions of the quantiles
plot(c(bin_labels,11),control_quantiles)
plot(c(bin_labels,11),strike_quantiles)

#cut the data by quantile
plot_quantiles_control$cut<-cut(control_probs,control_quantiles,bin_labels)
plot_quantiles_strike$cut<-cut(strike_probs,strike_quantiles,bin_labels)

#confirm that the data are pretty evenly distributed among bins
table(plot_quantiles_control$cut)
table(plot_quantiles_strike$cut)

########################################
#now use those quantiles to calculate proporitons of real data with confidence intervals
#############################################
#now create confidence intervals based on the binomial distribution

#aggregate the number of damaged and killed trees by 5m distance class
control_dmg<-aggregate(dmg.binary~cut,data=plot_quantiles_control,FUN="sum")
control_total<-aggregate(dmg.binary~cut,data=plot_quantiles_control,FUN="length")

#create a complete dataframe from all three variables
control_dmg_prop<-data.frame(control_total,"dmg"=control_dmg$dmg.binary)
colnames(control_dmg_prop)<-c("Quantile","Total","dmg")
str(control_dmg_prop)

#use DescTools to calculate the binomial CI
CI_dataframe_control_dmg<-BinomCI(control_dmg_prop$dmg,control_dmg_prop$Total,
                                   conf.level=0.95, method = "clopper-pearson")
plottingCIs_control_dmg<-data.frame(as.data.frame(CI_dataframe_control_dmg),"plotting.numeric" = seq(.025,.975,.05),
                                     "plotting.label" = bin_labels, "data.type" = rep("Real proportions",10), "control_strike" = rep("Control",10))
str(plottingCIs_control_dmg)

#aggregate the number of damaged and killed trees by 5m distance class
strike_dmg<-aggregate(dmg.binary~cut,data=plot_quantiles_strike,FUN="sum")
strike_total<-aggregate(dmg.binary~cut,data=plot_quantiles_strike,FUN="length")

#create a complete dataframe from all three variables
strike_dmg_prop<-data.frame(strike_total,"dmg"=strike_dmg$dmg.binary)
colnames(strike_dmg_prop)<-c("Quantile","Total","dmg")
str(strike_dmg_prop)

#use DescTools to calculate the binomial CI
CI_dataframe_strike_dmg<-BinomCI(strike_dmg_prop$dmg,strike_dmg_prop$Total,
                                  conf.level=0.95, method = "clopper-pearson")
plottingCIs_strike_dmg<-data.frame(as.data.frame(CI_dataframe_strike_dmg),"plotting.numeric" = seq(.025,.975,.05),
                                    "plotting.label" = bin_labels, "data.type" = rep("Real proportions",10), "control_strike" = rep("Strike",10))
str(plottingCIs_strike_dmg)


###########################################
#calculate confidence intervals of probability with "groupwiseMean"
control_CIsinter<-groupwiseMean(control_probs~cut, data = plot_quantiles_control, conf = .95)
control_CIs<-data.frame(control_CIsinter[-11,c(1,3,5,6)],data.type=rep("Predicted probabilities",10), "control_strike" = rep("Control",10))

strike_CIsinter<-groupwiseMean(strike_probs~cut, data = plot_quantiles_strike, conf = .95)
strike_CIs<-data.frame(strike_CIsinter[-11,c(1,3,5,6)],data.type=rep("Predicted probabilities",10), "control_strike" = rep("Strike",10))

#################################################


##########################
#now create 2 plotting dataframe
#rename the columns to facilitate combination
colnames(control_CIs)[c(1,3,4)]<-c("plotting.label","upr.ci","lwr.ci")
colnames(strike_CIs)[c(1,3,4)]<-c("plotting.label","upr.ci","lwr.ci")
colnames(plottingCIs_control_dmg)[1]<-"Mean"
colnames(plottingCIs_strike_dmg)[1]<-"Mean"

#add plotting label variable and change label factor to numeric
control_CIs$plotting.label<-as.numeric(as.character(control_CIs$plotting.label))
strike_CIs$plotting.label<-as.numeric(as.character(strike_CIs$plotting.label))
control_CIs$plotting.numeric<-plottingCIs_control_dmg$plotting.numeric
strike_CIs$plotting.numeric<-plottingCIs_strike_dmg$plotting.numeric

#bind rows to combine data frames
control_plotting<-bind_rows(plottingCIs_control_dmg,control_CIs)
str(control_plotting)
strike_plotting<-bind_rows(plottingCIs_strike_dmg,strike_CIs)
str(strike_plotting)

#create a combined plotting variable
combined_plotting<-bind_rows(strike_plotting,control_plotting)

################################
#create ggplot theme
theme_basis<-theme(axis.title = element_text(family = "Arial", color="black", face="bold", size=16))+
  theme(axis.text.y=element_text(family = "Arial", 
                                 colour="black", face ="bold", size = 14)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),axis.ticks = element_line(colour="black"))+
  theme(axis.text.x = element_text(family = "Arial", colour="black", face ="bold", size = 14))

#now plot data
epiidmg_glmfit<-ggplot(data = combined_plotting, aes(y = Mean, x = plotting.label, shape = data.type,
                                                    color = control_strike))+
  geom_errorbar(aes(ymax = upr.ci, ymin = lwr.ci),width = .2, position = position_dodge(width = .6))+
  geom_point(size = 2, position = position_dodge(width = .6))+facet_wrap(~control_strike, nrow = 2, scales = "free_x")+
  scale_shape_manual(values = c(1,16)) + theme_basis +
  theme(strip.text.x = element_text(family = "Arial",colour="black", face ="bold",size=14),
        strip.background = element_blank())+
  guides(color = FALSE)+
  theme(legend.title=element_blank(),legend.background = element_blank(),
        legend.text=element_text(family = "Arial",colour="black", face ="bold",size=12), 
        legend.position=c(.7,.9),   legend.spacing.y = unit(.001, 'cm'),
        legend.key = element_blank(),legend.key.height=unit(.9,"line"))+
  scale_x_continuous(name = "Quantiles of Epiphyte Damage", breaks = seq(0,10,1),labels = c("0.0","","0.2","","0.4","","0.6","","0.8","","1.0"))+
  scale_y_continuous(name = "Probability of Damage",breaks = seq(0,1,.1),labels = c("0.0","","0.2","","0.4","","0.6","","0.8","","1.0"))

epiidmg_glmfit#this looks ok - some variation along the way, but the general patterns are being captured!

ggsave("epiidmg_glmfit.tiff",epiidmg_glmfit,dpi = 600,width = 3.3, height = 3, scale = 1.6, compression = "lzw")
getwd()
