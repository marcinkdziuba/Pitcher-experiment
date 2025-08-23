# loading packages
library(readxl)
library(ggplot2)
library(mgcv)
library(emmeans)
library(dplyr)
library(tidyr)
library(Rmisc)
library(gridExtra)
library(cowplot)
library(patchwork)
library(betareg)
library(dunn.test)
library(nlme)
library(car)

###############################################################################
################ INDIVIDUAL LEVEL EXP #########################################
###############################################################################

########### REPRODUCTION ##############

repMGvM <- read_excel("Individual reprod.xlsx") # Load data
repMGvM<-repMGvM[repMGvM$clutch<15,]
# data curation
repMGvM$NB<-as.numeric(repMGvM$NB)
repMGvM$age<-as.numeric(repMGvM$age)
repMGvM$parasite<-as.factor(repMGvM$parasite)
repMGvM<-na.omit(repMGvM)
# change parasite level names to match the population experiment
repMGvM$parasite <- factor(repMGvM$parasite, 
                                levels = c("C", "MET", "MG"),
                                labels = c("C", "MB", "OP"))
#data summary overview
repMGvMse<-summarySE(data=repMGvM, measurevar = "NB", groupvars = c("parasite","clutch"),na.rm=T)
repMGvMse

# gam model with interaction of parasite treatment and clutch number, with common smoother for all treatments
MGvM_M1<-gam(NB~parasite*clutch+s(clutch),data = repMGvM)
anova(MGvM_M1) #stats report
summary(MGvM_M1)

# model diagnostics
op <- par(mfrow = c(2, 2))
gam.check(MGvM_M1)
par(mfrow = c(1, 1))
vis.gam(MGvM_M1, theta = 150, color = "heat") # visualisation of gam model

# The pattern in response vs Fitted is due to smaller variation in number of 
# offspring produced in first few clutches, when the clutches are small 

# Use emmeans to get the emmeans object for parasite groups across levels of clutch
emmsRep <- emmeans(MGvM_M1, ~ parasite | clutch, at = list(clutch = seq(min(repMGvM$clutch), max(repMGvM$clutch), by = 1)))

# Perform pairwise comparisons
contrast(emmsRep, method = "pairwise", adjust = "bonferroni")
contrastRep<-contrast(emmsRep, method = "pairwise", adjust = "bonferroni")
summary(contrastRep)
emmsRep
# Exporting post-hoc comparisons into .csv files
pairwise_comparisons_Rep <- as.data.frame(summary(contrastRep))
write.csv(pairwise_comparisons_Rep, file = "pairwise_comparisons_Rep.csv", row.names = FALSE)
write.csv(emmsRep, file = "emmsRep.csv", row.names = FALSE)

# GAM generally looks good except for a bit smaller variance for the early clutches,
# so just to confirm the patterns, let's run a non-parametric test as well
# analyzing parasite within each clutch
repMGvM$fclutch<-as.factor(repMGvM$clutch)
RP_KW_c1<-kruskal.test(NB[repMGvM$fclutch=="1"]~parasite[repMGvM$fclutch=="1"],data = repMGvM)
RP_KW_c1
dunn.test(repMGvM$NB[repMGvM$fclutch=="1"], repMGvM$parasite[repMGvM$fclutch=="1"], method = "bonferroni")

RP_KW_c2<-kruskal.test(NB[repMGvM$fclutch=="2"]~parasite[repMGvM$fclutch=="2"],data = repMGvM)
RP_KW_c2
dunn.test(repMGvM$NB[repMGvM$fclutch=="2"], repMGvM$parasite[repMGvM$fclutch=="2"], method = "bonferroni")

RP_KW_c3<-kruskal.test(NB[repMGvM$fclutch=="3"]~parasite[repMGvM$fclutch=="3"],data = repMGvM)
RP_KW_c3
dunn.test(repMGvM$NB[repMGvM$fclutch=="3"], repMGvM$parasite[repMGvM$fclutch=="3"], method = "bonferroni")

RP_KW_c4<-kruskal.test(NB[repMGvM$fclutch=="4"]~parasite[repMGvM$fclutch=="4"],data = repMGvM)
RP_KW_c4
dunn.test(repMGvM$NB[repMGvM$fclutch=="4"], repMGvM$parasite[repMGvM$fclutch=="4"], method = "bonferroni")

RP_KW_c5<-kruskal.test(NB[repMGvM$fclutch=="5"]~parasite[repMGvM$fclutch=="5"],data = repMGvM)
RP_KW_c5

RP_KW_c6<-kruskal.test(NB[repMGvM$fclutch=="6"]~parasite[repMGvM$fclutch=="6"],data = repMGvM)
RP_KW_c6
dunn.test(repMGvM$NB[repMGvM$fclutch=="6"], repMGvM$parasite[repMGvM$fclutch=="6"], method = "bonferroni")

RP_KW_c7<-kruskal.test(NB[repMGvM$fclutch=="7"]~parasite[repMGvM$fclutch=="7"],data = repMGvM)
RP_KW_c7
dunn.test(repMGvM$NB[repMGvM$fclutch=="7"], repMGvM$parasite[repMGvM$fclutch=="7"], method = "bonferroni")
# qualitatively the non-parametric test outcomes looks similar, except it's a bit more
# conservative so the difference between C and OP is visible in later clutches in comparison to GAM

#### FIUGRE #######

# a few lines of code to design new dataframe with values predicted with the GAM model
#create new dataframe
fitREP_MGvM1 <- data.frame()

# Loop through each level of the parasite group
for (p in levels(repMGvM$parasite)) {
  # Subset the data for the current parasite level
  subset_data <- subset(repMGvM, parasite == p)
  
  # Create a sequence of clutch values specific to the current parasite level
  clutch_seq <- seq(min(subset_data$clutch), max(subset_data$clutch), length.out = 1000)
  
  # Create a data frame for the current parasite level
  temp_df <- expand.grid(clutch = clutch_seq, parasite = p)
  
  # Combine with the main data frame
  fitREP_MGvM1 <- rbind(fitREP_MGvM1, temp_df)
}

# Create new column with densities predicted for each timepoint using the GAM model
predNB <- predict(MGvM_M1, newdata = fitREP_MGvM1, se.fit = TRUE)
fitREP_MGvM1$NB <- predNB$fit
fitREP_MGvM1$se <- predNB$se.fit



# plot the fitted values and errorbars
REP1<-
ggplot(fitREP_MGvM1, aes(x=clutch,y=NB,group = parasite,color = parasite))+
  geom_smooth(data=fitREP_MGvM1,aes(ymin=NB-se, ymax=NB+se),stat='identity',alpha = .3)+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange'))+
  scale_y_continuous(limits = c(0, 20), expand = c(0,0), name = "# offspring in clutch")+
  scale_x_continuous(breaks = c(0,3,6,9,12),expand=c(0,0),name="Clutch (reproduction event)")+
  geom_rect(data=fitREP_MGvM1, mapping=aes(xmin=3.8, xmax=14, ymin=0.75, ymax=1.05), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitREP_MGvM1, mapping=aes(xmin=1.8, xmax=4, ymin=0.2, ymax=0.5), fill="red", color="red", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.1,0.85),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs( color = "Parasite")


#### SURVIVAL #################################################################
#loading data and data curation
mortMGvM <- read_excel("Individual mort.xlsx")
mortMGvM$lifespan<-as.numeric(mortMGvM$lifespan)
mortMGvM<-na.omit(mortMGvM)
mortMGvM$parasite<-as.factor(mortMGvM$parasite)
mortMGvM$parasite <- factor(mortMGvM$parasite, 
                            levels = c("C", "MET", "MG"),
                            labels = c("C", "MB", "OP"))

# determines that variance is variable between parasite groups
varIDMort <- varIdent(form=~1 | parasite) 
# generalized least squares with variance estimated independently for each parasite level
LSM3<-gls(lifespan~parasite, data=mortMGvM, weights = varIDMort)
anova(LSM3)
summary(LSM3)
#diagnostics
hist(resid(LSM))
plot(LSM)

#contrasts
emmLSM <- emmeans(LSM, ~ parasite)
contrast(emmLSM, method = "pairwise", adjust = "bonferroni")

# the model has some small overdispersion, so just for safety let's analyze the same data
# with a non-parametric model
# Kruskal-Wallis test
KW_LSM<-kruskal.test(lifespan~parasite, data=mortMGvM)
KW_LSM
# Post-hoc
dunn.test(mortMGvM$lifespan,mortMGvM$parasite, method = "bonferroni")
# very similar resilts

####### FIGURE ############
LS_SE<-summarySE(data = mortMGvM, measurevar = "lifespan", groupvars = "parasite", na.rm=T)

LSplot<-
  ggplot(data = LS_SE, aes(x = parasite, y = lifespan))+
  geom_point(data=mortMGvM,shape=21, aes(x=parasite,y=lifespan, color = parasite, fill=parasite),
             position = position_jitter(width = 0.15, height = 0))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=lifespan-se, ymax=lifespan+se), width=.1)+
  scale_color_manual(values = c('dodgerblue', 'red','darkorange'))+
  scale_fill_manual(values = c('dodgerblue', 'red','darkorange'))+
  scale_y_continuous(limits= c(0,80),expand = c(0, 0), name = "Lifespan [days]")+
  scale_x_discrete(name="Parasite")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(fill="Treatment", color = "Treatment")



combined_plot_LT <- REP1 + LSplot + plot_layout(widths = c(2, 1)) + 
  plot_annotation(tag_levels = 'A')
combined_plot_LT <- combined_plot_LT & theme(plot.tag = element_text(size = 28))
ggsave("Fig1.png",combined_plot_LT, width = 12, height = 6) 

############ AGE AT REPRODUCTION ###############################

repMGvM_4C<-repMGvM[repMGvM$clutch<5,] # create a new database by subsetting first four clutches
                                       # of the individual-level experiment
repageSE<-summarySE(data=repMGvM_4C, measurevar = "age", groupvars = c("fclutch","parasite"))

varId<-varIdent(form = ~1 | parasite)
M_lme <- lme(age ~ parasite ,
             random = ~1 | clutch,         # use your grouping variable here
             data = repMGvM_4C,
             weights = varId)
summary(M_lme)
anova(M_lme)

plot(M_lme)
qqPlot(M_lme$resid)
emm_M_lme <- emmeans(M_lme, ~ parasite)
contrast(emm_M_lme, method = "pairwise", adjust = "bonferroni")
emmip(M_lme, ~ parasite)

### FIGURE ###
AgePlot<-
  ggplot(repageSE,aes(x=fclutch,y=age,color=parasite))+
  geom_point(position = position_dodge(.8))+
  geom_errorbar(aes(ymin=age-se,ymax=age+se),position = position_dodge(.8))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange'))+
  scale_y_continuous(limits = c(5,20),name="Age at reproduction")+
  scale_x_discrete(name="# clutch")+
  theme_classic()+
  theme(legend.position = c(0.15,0.8),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Parasite")
AgePlot
ggsave("FigS1.png",AgePlot, width = 4, height = 4) 

###############################################################################
################ POPULATION EXPERIMENT ########################################
###############################################################################

# loading database
full_data <- read_excel("Population data.xlsx")
full_data$treatment<-as.factor(full_data$treatment)

# create extra variables (i.e., convert counts to densities, calculate ratios)
# calculate the number of organisms per liter (from raw count data)
full_data$UA_L<-((full_data$'UA'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # uninfected adults
full_data$UJ_L<-((full_data$'UJ'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # uninfected juveniles
full_data$UM_L<-((full_data$'UM'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # uninfected males
full_data$OPA_L<-((full_data$'OPA'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii infected adults
full_data$OPJ_L<-((full_data$'OPJ'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii infected juveniles
full_data$OPM_L<-((full_data$'OPM'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii infected males
full_data$MBA_L<-((full_data$'MBA'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # M. bicuspidata infected adults
full_data$MBJ_L<-((full_data$'MBJ'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # M. bicuspidata infected juveniles
full_data$MBM_L<-((full_data$'MBM'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # M. bicuspidata infected males
full_data$OPMBA_L<-((full_data$'OPMBA'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii and M. bicuspidata coinfected adults
full_data$OPMBJ_L<-((full_data$'OPMBJ'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii and M. bicuspidata coinfected juveniles
full_data$OPMBM_L<-((full_data$'OPMBM'*full_data$'total volume ml')/full_data$'subsample volume ml')/5 # O. pajunii and M. bicuspidata coinfected males

full_data$density<-rowSums(full_data[,c('UA_L','UJ_L','UM_L','OPA_L','OPJ_L','OPM_L',
                                          'MBA_L','MBJ_L','MBM_L','OPMBA_L','OPMBJ_L','OPMBM_L' )],na.rm=T) # total density of animals per pitcher

full_data$uninfected<-rowSums(full_data[,c('UA_L','UJ_L','UM_L')],na.rm=T) # total density of uninfected
full_data$OPinfected<-rowSums(full_data[,c('OPA_L','OPJ_L','OPM_L')],na.rm=T) # total density of O. pajunii infected
full_data$MBinfected<-rowSums(full_data[,c('MBA_L','MBJ_L','MBM_L')],na.rm=T) # total density of M. bicuspidata infected
full_data$coinfected<-rowSums(full_data[,c('OPMBA_L','OPMBJ_L','OPMBM_L' )],na.rm=T) # total density of coinfected infected
full_data$adults<-rowSums(full_data[,c('UA_L','OPA_L','MBA_L','OPMBA_L' )],na.rm=T) # total density of adults
full_data$juveniles<-rowSums(full_data[,c('UJ_L','OPJ_L','MBJ_L','OPMBJ_L' )],na.rm=T) # total density of juveniles
full_data$males<-rowSums(full_data[,c('UM_L','OPM_L','MBM_L','OPMBM_L' )],na.rm=T) # total density of males
full_data$OPprev<-full_data$OPinfected/full_data$density # O. pajunii prevalence
full_data$MBprev<-full_data$MBinfected/full_data$density # M. bicuspidata prevalence
full_data$COprev<-full_data$coinfected/full_data$density # coinfection prevalence
full_data$adultOPprev<-full_data$OPA_L/full_data$adults # O. pajunii prevalence among adults
full_data$adultMBprev<-full_data$MBA_L/full_data$adults # M. bicuspidata prevalence among adults
full_data$adultCOprev<-full_data$OPMBA_L/full_data$adults # coinfection prevalence among adults
full_data$OPtotalprev<-rowSums(full_data[,c('OPprev','COprev')], na.rm = T) # all Daphnia infected with O. pajunii (single and coinfection), total prevalence
full_data$MBtotalprev<-rowSums(full_data[,c('MBprev','COprev')], na.rm = T) # all Daphnia infected with M. bicuspidata (single and coinfection), total prevalence
full_data$adultOPtotalprev<-rowSums(full_data[,c('adultOPprev','adultCOprev')], na.rm = T) # all Daphnia infected with O. pajunii (single and coinfection), prevalence among adults
full_data$adultMBtotalprev<-rowSums(full_data[,c('adultMBprev','adultCOprev')], na.rm = T) # all Daphnia infected with M. bicuspidata (single and coinfection), prevalence among adults


###############################################################################
########### POPULATION GROWTH RATE ############################################
###############################################################################
# Subsetting database to get first 30 days this time
data_first23<- full_data[full_data$day<30,]
data_first23$treatment<- as.factor(data_first23$treatment)


# initial plot of growth curves
ggplot(data=data_first23, aes(x=day, y=density, group=treatment, colour=treatment))+
  geom_smooth(span = 0.8, aes(ymin = ifelse(after_stat(ymin) <0,0,after_stat(ymin))))+
  #geom_point()+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','#FF00FF'))+
  #scale_y_continuous(limits = c(0,2000))+
  theme_bw()


# growth rate calculation
# create new dataframe with the subset of variables
data_GR<-data.frame("treatment"=data_first23$treatment,
                    "day"=data_first23$day,
                    "rep"=data_first23$rep,
                    "density"=data_first23$density)
data_GR$day<-as.factor(data_GR$day) #change day into factor 

#pivot table wide (to get density in each day in a separate column)
#that helps doing math on columns
data_GRwide <- pivot_wider(data_GR,
                           id_cols = c("treatment","rep"),
                           names_from = day, 
                           values_from = density)

names(data_GRwide)[names(data_GRwide) == "0"] <- "day0" #change column name
names(data_GRwide)[names(data_GRwide) == "23"] <- "day23" #change column name

# calculation of pop. growth rate r with equation r=log(density_t(n)-density_t(0))/time
data_GRwide$r<-log(data_GRwide$day23-data_GRwide$day0)/23 

# because the growth rates are measured before day 31 (dosing M. bicuspidata)
# groups C and MB are effectively the same, and so are groups OP and OPMB
# so for the growth rate computation these two pairs will be merged and 
# treated as one treatment

# adding variable parasite to discriminate between groups dosed and not dosed with O. pajunii
data_GRwide$parasite <- with(data_GRwide, ifelse(treatment == "C", "C",
                                          ifelse(treatment == "MB", "C", 
                                          ifelse(treatment == "OP", "OP",
                                          ifelse(treatment == "OPMB", "OP", NA)))))

# Growth rate comparison between groups
# attemnt on using a linear model
GRm1<-lm(r~parasite, data = data_GRwide)
anova(GRm1)
plot(GRm1) # trouble with the outliers

# linear model after removing the obvious outlier #39
data_GRwide2<-data_GRwide[-39,]
GRm2<-lm(r~parasite, data = data_GRwide2)
anova(GRm2)
plot(GRm2) # non-normal distribution of residuals

# a non-parametric test due to influential outliers and problems with the
# assumptions of normality of resid. distribution and heterogeneity of variances
# use full data (with the outlier #39)
KWr<-kruskal.test(r~parasite, data = data_GRwide)
KWr

##### FIGURE ####

#summarize r accross replicates (to obtain error estimates)
GR_SE<-summarySE(data = data_GRwide, measurevar = "r", groupvars = c("parasite"))

#plot
GRplot<-
  ggplot(data = GR_SE, aes(x = parasite, y = r))+
  geom_point(data=data_GRwide,shape=21, aes(x=parasite,y=r, color = treatment, fill=treatment),
             position = position_jitter(width = 0.15, height = 0))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=r-se, ymax=r+se), width=.1)+
  scale_color_manual(values = c('dodgerblue', 'red','darkorange','green3'))+
  scale_fill_manual(values = c('dodgerblue', 'red','darkorange','green3'))+
  scale_y_continuous(limits = c(0, 0.4), name = "Population growth rate r")+
  scale_x_discrete(name="Parasite")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(fill="Treatment", color = "Treatment")
GRplot


###############################################################################
########### POPULATION DENSITY ################################################
###############################################################################

# loess figure of Daphnia density in treatments
# vertical line indicates the day at which the MB and OPMB pitchers were dosed with M. bicuspidata spores
ggplot(data=full_data, aes(x=day, y=density, group=treatment, colour=treatment))+
  geom_smooth(span = 0.6, aes(ymin = ifelse(after_stat(ymin) <0,0,after_stat(ymin))))+
  geom_vline(xintercept = 31)+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','#FF00FF'))+
  theme_bw()


# create new dataframe that is a subset of the old database trimmed after day 23 
# start at day 30 corresponds to when populations reached peak densities
# and to the timepoint when M. bicuspidata was dosed to the system (day 31)

data_past23<- full_data[full_data$day>23,] 
data_past23$treatment<- as.factor(data_past23$treatment)

# GAM model on Daphnia density
# interaction between treatment and day and separate smoothers for each treatment
dens_gam <- gam(density ~ treatment *day+
             s(day, by = as.numeric(treatment == "C")) + 
             s(day, by = as.numeric(treatment == "MB")) +
             s(day, by = as.numeric(treatment == "OP")) + 
             s(day, by = as.numeric(treatment == "OPMB")), 
           data = data_past23)

summary(dens_gam)
anova(dens_gam)
# significant interaction between day and treatment; 
# and significant smoothers for C and MB

vis.gam(dens_gam, theta = 120, color = "heat") # visualisation of gam model
vis.gam(dens_gam, theta = 220, color = "heat") # visualisation of gam model, another angle

# model diagnostics
op <- par(mfrow = c(2, 2))
gam.check(dens_gam)
par(mfrow = c(1, 1))

res_dens_gam <- resid(dens_gam)
coplot(res_dens_gam ~ day | treatment, data = data_past23)

# Pairwise comparisons and contrasts between the models accros time
EM_dens_gam <- emmeans(dens_gam, ~ treatment | day, at = list(day = seq(min(data_past23$day), max(data_past23$day), by = 1)))
contrast(EM_dens_gam, method = "pairwise", adjust = "bonferroni") # pairwise contrast comparison with bonferroni correction
contrasts_dens_gam<-contrast(EM_dens_gam, method = "pairwise", adjust = "bonferroni")

# Exporting post-hoc comparisons into .csv files
pairwise_comparisons_dens_gam <- as.data.frame(summary(contrasts_dens_gam))
write.csv(pairwise_comparisons_dens_gam, file = "pairwise_comparisons_dens_gam.csv", row.names = FALSE)
write.csv(EM_dens_gam, file = "EM_dens_gam.csv", row.names = FALSE)

##### FIGURE ####

# create dataframe with a sequence of time from day 30 to 93 (last day of the experiment) for each treatment
fitted_23 <- expand.grid(day = seq(min(data_past23$day), max(data_past23$day), length.out = 1000),
                      treatment = levels(data_past23$treatment))
# create new column with densities predicted for each timepoint using the GAM model
predictions_dens <- predict(dens_gam, newdata = fitted_23, se.fit = TRUE)
fitted_23$density <- predictions_dens$fit
fitted_23$se <- predictions_dens$se.fit

# plot past 23 days
fig2b<-
ggplot(fitted_23, aes(x=day,y=density,group = treatment,color = treatment))+
  geom_smooth(data=fitted_23,aes(ymin=density-se, ymax=density+se),stat='identity',alpha = .3)+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'))+
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=64, ymin=5, ymax=10), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=61, ymin=13, ymax=18), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=71, xmax=83, ymin=320, ymax=325), fill="red", color="red", alpha=0.5) +
  scale_y_continuous(expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.8,0.2),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")
fig2b


#### FIGURE FOR THE ENTIRE SPAN OF THE EXPERIMENT ###
#this figure will have two parts:
#first part - gam modeled pop. increase curves (first 30 days - not statistically analyzed, no error estimates)
#second part - gam modeled pop. density after the peaks (post day 30) used in the analysis

# remove the extreme outlier from the model for the early stage growth curve modelling
full_data2<-full_data[!(full_data$treatment =="OPMB" & full_data$rep == "9"),]

# generate gam for the entire experiment
dens_gam_full <- gam(density ~ treatment *day+
                  s(day, by = as.numeric(treatment == "C")) + 
                  s(day, by = as.numeric(treatment == "MB")) +
                  s(day, by = as.numeric(treatment == "OP")) + 
                  s(day, by = as.numeric(treatment == "OPMB")), 
                data = full_data2)

# create a new dataframe for fitted density data based on the dens_gam_full model
fitted_early <- expand.grid(day = seq(min(full_data2$day), max(full_data2$day), length.out = 1000),
                         treatment = levels(full_data2$treatment))

# create new column with densities predicted for each timepoint using the GAM model
fitted_early$density <- predict(dens_gam_full, newdata = fitted_early)

#change negative densities into 0 (densities cannot be negative, but
# gam smoothing can generate negative fitted values when the actual values are close to 0 )
fitted_early$density<-with(fitted_early, ifelse(density<=0,0,
                                              ifelse(density>0,density,0)))

# trim the fitted_early data to days 0 to <30
fitted_early<-fitted_early[fitted_early$day<30,]

# plotting two sets of GAMs, one that was used for Daphnia density analysis (starting at day 30)
# and one for days 0-30 to visualize the population increase curves
DENSplot<-
ggplot(fitted_23, aes(x=day,y=density,group = treatment,color = treatment))+
  geom_smooth(data=fitted_23,aes(ymin=density-se, ymax=density+se),stat='identity',alpha = .3)+
  geom_smooth(data=fitted_early,stat='identity',alpha = .3,linetype="longdash")+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'))+
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=64, ymin=5, ymax=10), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=61, ymin=13, ymax=18), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=71, xmax=83, ymin=320, ymax=325), fill="red", color="red", alpha=0.5) +
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  scale_y_continuous(expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.1,0.7),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")
DENSplot

combined_plot <- GRplot + DENSplot +  plot_layout(widths = c(1, 3)) + 
  plot_annotation(tag_levels = 'A')
combined_plot <- combined_plot & theme(plot.tag = element_text(size = 28))
combined_plot
ggsave("Fig2.png",combined_plot, width = 12, height = 6) 


###############################################################################
################ M. bicuspidata prevalence ####################################
###############################################################################

# Create new dataframe with just data after day 31 when M. bicuspidata was doesed
# and just treatments dosed with M. bicuspidata
data_MBprev<-full_data[full_data$day>31 & full_data$treatment=="MB" | full_data$day>31 &full_data$treatment == "OPMB",]
data_MBprev$day<-as.numeric(data_MBprev$day)
data_MBprev$treatment<-as.factor(data_MBprev$treatment)

# beta regression model
PMbeta<-betareg(MBtotalprev~day*treatment, data = data_MBprev)
summary(PMbeta)
plot(PMbeta) # diagnostics
hist(resid(PMbeta)) # histogram of residuals

# post-hoc pairwise comparison at each day with emmeans
emmsPMbeta <- emmeans(PMbeta, ~ treatment | day, at = list(day = seq(min(data_MBprev$day), max(data_MBprev$day), by = 1)))
contrast(emmsPMbeta, method = "pairwise", adjust = "bonferroni")
# post-hoc visualization
emmip(emmsPMbeta, treatment ~ day, at = list(day = seq(min(data_MBprev$day), max(data_MBprev$day), by = 1)), CIs = TRUE)
summary(contrast(emmsPMbeta, method = "pairwise", adjust = "bonferroni"))

##### FIGURE ####
#prepare the new dataframe for predicted prevalence
data_MBprev$treatment<-as.factor(data_MBprev$treatment)
MBprev <- expand.grid(day = seq(min(data_MBprev$day), max(data_MBprev$day), length.out = 100),
                      treatment = c("MB","OPMB"))

# Get predictions from the model
MBprev <- MBprev %>%
  mutate(predicted = predict(PMbeta, MBprev, type = "response"))

MBplot<-
ggplot(data_MBprev, aes(x=day, y=MBprev, group=treatment, color=treatment, shape = treatment)) +
  geom_point(size=1.5, alpha=0.5) +
  scale_shape_manual(values = c(15,16))+
  geom_line(data=MBprev, aes(x=day, y=predicted), size=1) + 
  scale_colour_manual(values=c('red','green3')) +
  scale_x_continuous(limits=c(29.7, 94), expand=c(0, 0), name="Day") +
  scale_y_continuous(limits=c(NA, 1), expand=c(0, 0), name="M. bicuspidata prevalence") +
  geom_rect(data=data_MBprev, mapping=aes(xmin=56, xmax=93, ymin=0.965, ymax=0.98), fill="red", color="red", alpha=0.5) +
  theme_classic() +
  theme(
    legend.position=c(0.1, 0.82),
    legend.background=element_blank(),
    axis.title=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_text(size=14, color="black"),
    legend.text=element_text(size=14, color="black")
  ) +
  labs(color="Treatment", shape="Treatment") 
MBplot
################################################################################
######################### O. pajunii prevalence ################################
################################################################################

#new dataframe with just O. pajunii exposed treatments
data_OPprev<-full_data[ full_data$treatment=="OP" | full_data$treatment == "OPMB",]
data_OPprev$treatment<-as.factor(data_OPprev$treatment)
# Drop unused factor levels
data_OPprev$treatment <- droplevels(data_OPprev$treatment)
#remove days <30 to match the density analysis and the M. bicuspidata analysis
data_OPprev<-data_OPprev[data_OPprev$day>23,]

#gam model with single smoother for both groups
POPgam<-gam(OPtotalprev~s(day)+day*treatment,data = data_OPprev)
summary(POPgam)
anova(POPgam)

#model diagnostics
gam.check(POPgam)
op <- par(mfrow = c(2, 2))
gam.check(POPgam)
par(op)

# Post-hoc pairwise comparisons
emmsPOPgam <- emmeans(POPgam, ~ treatment | day, 
          at = list(day = seq(min(data_OPprev$day), max(data_OPprev$day), by = 1)))
contrast(emmsPOPgam, method = "pairwise", adjust = "bonferroni")
emmip(POPgam, treatment ~ day, at = list(day = seq(min(data_OPprev$day), max(data_OPprev$day), by = 1)), CIs = TRUE)

##### FIGURE ####
OPprev <- expand.grid(day = seq(min(data_OPprev$day), max(data_OPprev$day), length.out = 100),
                      treatment = levels(data_OPprev$treatment))

OPprev <- OPprev %>%
  mutate(predicted = predict(POPgam, OPprev, type = "response"))

OPplot<-
ggplot(data_OPprev, aes(x=day,y=OPprev,group = treatment, color = treatment, shape = treatment))+
  geom_point(size = 1.5, alpha = .5)+
  scale_shape_manual(values = c(17,16))+
  geom_line(data = OPprev, aes(x = day, y = predicted), linewidth = 1) +
  scale_colour_manual(values = c( 'darkorange','green3'))+
  scale_x_continuous(limits=c(29.7, 94),expand = c(0,0), name = "Day")+
  scale_y_continuous(limits=c(NA, 1),expand = c(0,0),name = "O. pajunii prevalence")+
  geom_rect(data=data_OPprev, mapping=aes(xmin=44, xmax=93, ymin=0.965, ymax=0.98), fill="darkorange", color="darkorange", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.1,0.82),
        legend.background = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment", shape="Treatment")

combined_plot2 <- MBplot + OPplot + plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = 'A')
combined_plot2 <- combined_plot2 & theme(plot.tag = element_text(size = 28))
combined_plot2
ggsave("Fig3.png",combined_plot2, width = 6, height = 8) 


###########################################################################
########### DEMOGRAPHICS ##################################################
###########################################################################

#new dataframe to use for demographics analysis, includes density of adults and juveniles
demo<-full_data[c("treatment","rep","day","adults","juveniles")]
demo$juvrat<-demo$juveniles/(demo$juveniles+demo$adults) # calculate juvenile ratio in populations
demo_past23<-demo[demo$day>23,] #remove days <30

# linear model of treatment*day interaction
juvratM1<-lm(juvrat~treatment*day,data=demo_past23)
summary(juvratM1)
anova(juvratM1)
par(mfrow=c(2,2))
plot(juvratM1)
par(mfrow=c(1,1))


emmjuvrat<-emmeans(juvratM1,~treatment)
contrast(emmjuvrat, method = "pairwise", adjust = "bonferroni")

# Get predictions from the model
JRM1<- expand.grid(day = seq(min(demo_past23$day), max(demo_past23$day), length.out = 100),
                   treatment = levels(demo_past23$treatment))

predictions_juvrat <- predict(juvratM1, newdata = JRM1, se.fit = TRUE)
JRM1$juvrat <- predictions_juvrat$fit
JRM1$se <- predictions_juvrat$se.fit

JRplot<-
ggplot(JRM1, aes(x = day, y = juvrat, group=treatment,fill = "grey")) +
  geom_ribbon(aes(ymin = juvrat - se, ymax = juvrat + se, fill = "grey"), alpha = 0.2) +
  geom_line(aes(y = juvrat, color = treatment), linewidth = 1) +
  ylab(label = "proportion of juveniles") +
  scale_color_manual(values = c('C' = 'dodgerblue', 'MB' = 'red', 'OP' = 'darkorange', 'OPMB' = 'green3')) +
  scale_fill_manual(values = c('C' = 'dodgerblue', 'MB' = 'red', 'OP' = 'darkorange', 'OPMB' = 'green3')) +
  scale_x_continuous(expand = c(0,0), name = "Day")+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
   theme_classic()+
  theme(legend.position = c(0.8,0.8),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")
JRplot
ggsave("Fig4.png",JRplot, width = 6, height = 4) 


###############################################################################
############################### ALGAE #########################################
###############################################################################

full_data$algae<- as.numeric(full_data$algae)
# Make new dataframe for algae analysis
algaedat<-full_data[c("treatment","day","rep","density","algae")]

# Make new variable "day" as factor
algaedat$fday<-as.factor(algaedat$day)

#get rid of day 0
algaedatNo0<-algaedat[algaedat$day>0,]

# GAM model with Tweedie distribution
gam_mod <- gam(algae ~  treatment +
                 s(day, by = as.numeric(treatment == "C")) + 
                 s(day, by = as.numeric(treatment == "MB")) +
                 s(day, by = as.numeric(treatment == "OP")) + 
                 s(day, by = as.numeric(treatment == "OPMB")), 
               family = tw(link = "log"), 
               data = algaedatNo0)
summary(gam_mod)
anova(gam_mod)
plot(gam_mod$residuals)

alg_gam <- emmeans(gam_mod, ~ treatment | day, at = list(day = seq(min(algaedatNo0$day), max(algaedatNo0$day), by = 1)))
contrast(alg_gam, method = "pairwise", adjust = "bonferroni") # pairwise contrast comparison with bonferroni correction
contrasts_alg_gam<-contrast(alg_gam, method = "pairwise", adjust = "bonferroni")

# Exporting post-hoc comparisons into .csv files
pairwise_comparisons_algae_gam <- as.data.frame(summary(contrasts_alg_gam))
write.csv(pairwise_comparisons_algae_gam, file = "pairwise_comparisons_algae_gam.csv", row.names = FALSE)

##### FIGURE ####

fitalg <- expand.grid(day = seq(min(algaedatNo0$day), max(algaedatNo0$day), length.out = 1000),
                      treatment = levels(algaedatNo0$treatment))
# create new column with densities predicted for each timepoint using the GAM model
predictions_alg <- predict(gam_mod, newdata = fitalg, se.fit = TRUE, type="response")
fitalg$alg <- predictions_alg$fit
fitalg$se <- predictions_alg$se.fit

algGAMplot<-
  ggplot(fitalg, aes(x=day,y=alg,group = treatment,color = treatment))+
  geom_smooth(data=fitalg,aes(ymin=alg-se, ymax=alg+se),stat='identity',alpha = .3)+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'))+
  geom_rect(data=fitalg, mapping=aes(xmin=25, xmax=64, ymin=16.8, ymax=17), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=77.5, xmax=78.5, ymin=16.8, ymax=17), fill="dodgerblue", color="green3", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=29, xmax=56, ymin=16.4, ymax=16.6), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=67, xmax=84, ymin=16.4, ymax=16.6), fill="dodgerblue", color="darkorange", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=70, xmax=84, ymin=16, ymax=16.2), fill="dodgerblue", color="red", alpha=0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(0,17), name = "Algae concentration [μg chl-a/l]")+
  scale_x_continuous(breaks=c(9,30,60,90),expand = c(0,0),limits = c(9,93), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.8,0.7),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")

algGAMplot
ggsave("FigS2.png",algGAMplot, width = 8, height = 6) 
