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


# Let's fit another model, that uses clutch as smoother fitted separately for each parasite treatment
MGvM_M2<-gam(NB~parasite+s(clutch, by=parasite),data = repMGvM)
anova(MGvM_M2) #stats report
summary(MGvM_M2)
concurvity(MGvM_M2, full = TRUE) # High concurvity values for MB are due to no data past clutch 4
                                 # They don't really indicate concurvity, can be ignored
op <- par(mfrow = c(2, 2))
gam.check(MGvM_M2)
par(mfrow = c(1, 1))

# Use emmeans to get the emmeans object for parasite groups across levels of clutch
emmsRep <- emmeans(MGvM_M2, ~ parasite | clutch, at = list(clutch = seq(min(repMGvM$clutch), max(repMGvM$clutch), by = 1)))

# Perform pairwise comparisons
contrast(emmsRep, method = "pairwise", adjust = "bonferroni")
contrastRep<-contrast(emmsRep, method = "pairwise", adjust = "bonferroni")
summary(contrastRep)
emmsRep
# Exporting post-hoc comparisons into .csv files
pairwise_comparisons_Rep <- as.data.frame(summary(contrastRep))
write.csv(pairwise_comparisons_Rep, file = "pairwise_comparisons_Rep.csv", row.names = FALSE)
write.csv(emmsRep, file = "emmsRep.csv", row.names = FALSE)

#### FIUGRE #######

# a few lines of code to design new dataframe with values predicted with the GAM model
#create new dataframe
fitREP_MGvM2 <- data.frame()

# Loop through each level of the parasite group
for (p in levels(repMGvM$parasite)) {
  # Subset the data for the current parasite level
  subset_data <- subset(repMGvM, parasite == p)
  
  # Create a sequence of clutch values specific to the current parasite level
  clutch_seq <- seq(min(subset_data$clutch), max(subset_data$clutch), length.out = 1000)
  
  # Create a data frame for the current parasite level
  temp_df <- expand.grid(clutch = clutch_seq, parasite = p)
  
  # Combine with the main data frame
  fitREP_MGvM2 <- rbind(fitREP_MGvM2, temp_df)
}

# Create new column with densities predicted for each timepoint using the GAM model
predNB <- predict(MGvM_M2, newdata = fitREP_MGvM2, se.fit = TRUE)
fitREP_MGvM2$NB <- predNB$fit
fitREP_MGvM2$se <- predNB$se.fit



# plot the fitted values and errorbars
REP1<-
ggplot(fitREP_MGvM2, aes(x=clutch,y=NB,group = parasite,color = parasite))+
  geom_smooth(data=fitREP_MGvM2,aes(ymin=NB-se, ymax=NB+se),stat='identity',alpha = .3)+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange'),breaks=c("C","MB","OP"),
                      labels=c("C","AM","OP"))+
  scale_y_continuous(limits = c(0, 20), expand = c(0,0), name = "# offspring in clutch")+
  scale_x_continuous(breaks = c(0,3,6,9,12),expand=c(0,0),name="Clutch (reproduction event)")+
  geom_rect(data=fitREP_MGvM2, mapping=aes(xmin=4.8, xmax=14, ymin=0.75, ymax=1.05), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitREP_MGvM2, mapping=aes(xmin=1.8, xmax=4, ymin=0.2, ymax=0.5), fill="red", color="red", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.15,0.88),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"))+
  labs( color = "Parasite")
REP1

# SI figure

fitREP_C<-fitREP_MGvM2[fitREP_MGvM2$parasite=="C",]
repMGvM_C<-repMGvM[repMGvM$parasite=="C",]

repC<-
  ggplot(fitREP_C, aes(x=clutch,y=NB))+
  geom_line(data=fitREP_C,stat='identity',color = "dodgerblue", size=3)+
  geom_line(data=repMGvM_C,(aes(x=clutch,y=NB, group = rep)),color = "dodgerblue",alpha=0.2)+
    scale_y_continuous(limits = c(0, 30), expand = c(0,0), name = "# offspring in clutch")+
    scale_x_continuous(limits = c(1, 14),breaks = c(0,3,6,9,12),expand=c(0,0),name="Clutch (reproduction event)")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
repC

fitREP_AM<-fitREP_MGvM2[fitREP_MGvM2$parasite=="MB",]
repMGvM_AM<-repMGvM[repMGvM$parasite=="MB",]

repAM<-
  ggplot(fitREP_AM, aes(x=clutch,y=NB))+
  geom_line(data=fitREP_AM,stat='identity',color = "red", size=3)+
  geom_line(data=repMGvM_AM,(aes(x=clutch,y=NB, group = rep)),color = "red",alpha=0.2)+
  scale_y_continuous(limits = c(0, 30), expand = c(0,0), name = "# offspring in clutch")+
  scale_x_continuous(limits = c(1, 14),breaks = c(0,3,6,9,12),expand=c(0,0),name="Clutch (reproduction event)")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
repAM


fitREP_OP<-fitREP_MGvM2[fitREP_MGvM2$parasite=="OP",]
repMGvM_OP<-repMGvM[repMGvM$parasite=="OP",]

repOP<-
  ggplot(fitREP_OP, aes(x=clutch,y=NB))+
  geom_line(data=fitREP_OP,stat='identity',color = "darkorange", size=3)+
  geom_line(data=repMGvM_OP,(aes(x=clutch,y=NB, group = rep)),color = "darkorange",alpha=0.2)+
  scale_y_continuous(limits = c(0, 30), expand = c(0,0), name = "# offspring in clutch")+
  scale_x_continuous(limits = c(1, 14),breaks = c(0,3,6,9,12),expand=c(0,0),name="Clutch (reproduction event)")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
repOP


design_indSI <- "
A
B
C
"
combined_plot_indSI <- repC + repAM + repOP + plot_layout(design = design_indSI)+
  plot_annotation(tag_levels = 'A')
combined_plot_indSI
ggsave("Fig1SI.tiff",combined_plot_indSI, width = 8, height = 8,dpi = 600)

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
hist(resid(LSM3))
plot(LSM3)

#contrasts
emmLSM <- emmeans(LSM3, ~ parasite)
contrast(emmLSM, method = "pairwise", adjust = "bonferroni")

####### FIGURE ############
LS_SE<-summarySE(data = mortMGvM, measurevar = "lifespan", groupvars = "parasite", na.rm=T)

LSplot<-
  ggplot(data = LS_SE, aes(x = parasite, y = lifespan))+
  geom_point(data=mortMGvM,shape=21, aes(x=parasite,y=lifespan, color = parasite, fill=parasite),
             position = position_jitter(width = 0.15, height = 0))+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin=lifespan-se, ymax=lifespan+se), width=.1)+
  scale_color_manual(values = c('dodgerblue', 'red','darkorange'),breaks=c("C","MB","OP"),
                     labels=c("C","AM","OP"))+
  scale_fill_manual(values = c('dodgerblue', 'red','darkorange'),breaks=c("C","MB","OP"),
                    labels=c("C","AM","OP"))+
  scale_y_continuous(limits= c(0,80),expand = c(0, 0), name = "Lifespan [days]")+
  scale_x_discrete(name="Parasite", ,breaks=c("C","MB","OP"),
                   labels=c("C","AM","OP"))+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(fill="Treatment", color = "Treatment")

REP1 <- REP1 + theme(plot.tag = element_text(size = 28))
LSplot <- LSplot + theme(plot.tag = element_text(size = 28))

combined_plot_LT <- REP1 + LSplot + plot_layout(widths = c(2, 1))+
  plot_annotation(tag_levels = 'A')
combined_plot_LT
ggsave("Fig1.tiff",combined_plot_LT, width = 8, height = 4,dpi=600) 

############ AGE AT REPRODUCTION ###############################

repMGvM_4C<-repMGvM[repMGvM$clutch<5,] # create a new database by subsetting first four clutches
                                       # of the individual-level experiment
repageSE<-summarySE(data=repMGvM_4C, measurevar = "age", groupvars = c("clutch","parasite"))

#mixed effects variance covariate
varId<-varIdent(form = ~1 | parasite)
M_lme_cv <- lme(age ~ parasite ,
             random = ~1 | clutch,         # use your grouping variable here
             data = repMGvM_4C,
             weights = varId)
summary(M_lme_cv)
anova(M_lme_cv)

plot(M_lme_cv)
qqPlot(M_lme_cv$resid)

emm_M_lme <- emmeans(M_lme_cv, ~ parasite)
contrast(emm_M_lme, method = "pairwise", adjust = "bonferroni")
emmip(M_lme_cv, ~ parasite)

### FIGURE ###
repageSE$clutch<-as.factor(repageSE$clutch)
AgePlot<-
  ggplot(repageSE,aes(x=clutch,y=age,color=parasite))+
  geom_point(position = position_dodge(.8))+
  geom_errorbar(aes(ymin=age-se,ymax=age+se),position = position_dodge(.8))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange'),breaks=c("C","MB","OP"),
                      labels=c("C","AM","OP"))+
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
ggsave("FigS1.tiff",AgePlot, width = 4, height = 4, dpi=600) 

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
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
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
dens_gam<-gam(density~treatment + s(day, by=treatment),data=data_past23)
anova(dens_gam)
summary(dens_gam)

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

concurvity(dens_gam)

# check for temporal autocorrelation (i.e., whether the observation of Daphnia density
# in time t+1 is correlated to the observation on time t)

data_past23$resid_dens <- residuals(dens_gam, type = "pearson")
data_past23$pitcherID <- interaction(data_past23$treatment,data_past23$rep,drop = TRUE,sep = "_")

series <- split(
  data_past23,
  interaction(
    data_past23$treatment,
    data_past23$pitcherID,
    drop = TRUE
  )
)

# Keep eligible groups and sort alphabetically
series <- series[vapply(series, nrow, integer(1)) >= 3]
series <- series[order(names(series))]

# Create an output folder with cross correlation figures
output_dir <- "residual_acf_jpg"
dir.create(output_dir, showWarnings = FALSE)

plots_per_page <- 12

if (length(series) > 0) {
  pages <- split(
    seq_along(series),
    ceiling(seq_along(series) / plots_per_page)
  )
  
  for (page in seq_along(pages)) {
    idx <- pages[[page]]
    
    jpeg(
      filename = file.path(
        output_dir,
        sprintf("residual_acf_page_%02d.jpg", page)
      ),
      width = 12,
      height = 9,
      units = "in",
      res = 300,
      quality = 95
    )
    
    par(
      mfrow = c(3, 4),
      mar = c(2.5, 2.5, 2, 0.5)
    )
    
    for (i in idx) {
      nm <- names(series)[i]
      d <- series[[i]]
      d <- d[order(d$day), ]
      
      acf(
        d$resid_dens,
        main = "",
        na.action = na.fail,
        ylim = c(-1, 1.25)
      )
      
      mtext(
        paste("Residual ACF:", nm),
        side = 3,
        line = -1.2,
        cex = 0.65,
        font = 2
      )
    }
    
    dev.off()
  }
}


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
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=64, ymin=5, ymax=10), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=61, ymin=13, ymax=18), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=71, xmax=83, ymin=320, ymax=325), fill="red", color="red", alpha=0.5) +
  #geom_line(data=data_past23,(aes(x=day,y=density,color=treatment, group = interaction(treatment, rep))),alpha=0.2,position = position_dodge(1))+
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
dens_gam_full <- gam(density ~ treatment + s(day, by=treatment),data = full_data2)

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
fitted_early2<-fitted_early
fitted_early2<-fitted_early2[fitted_early2$day<30,]

# plotting two sets of GAMs, one that was used for Daphnia density analysis (starting at day 30)
# and one for days 0-30 to visualize the population increase curves
DENSplot<-
ggplot(fitted_23, aes(x=day,y=density,group = treatment,color = treatment))+
  geom_smooth(data=fitted_23,aes(ymin=density-se, ymax=density+se),stat='identity',alpha = .3)+
  geom_smooth(data=fitted_early2,stat='identity',alpha = .3,linetype="longdash")+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=64, ymin=5, ymax=10), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=61, ymin=13, ymax=18), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=71, xmax=83, ymin=335, ymax=340), fill="red", color="red", alpha=0.5) +
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  scale_y_continuous(limits = c(0,350),expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.1,0.8),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 9, color = "black"))+
  labs(color="Treatment")
DENSplot

GRplot <- GRplot + theme(plot.tag = element_text(size = 28))
DENSplot <- DENSplot + theme(plot.tag = element_text(size = 28))
combined_plot <- GRplot + DENSplot +  plot_layout(widths = c(1, 3)) + 
  plot_annotation(tag_levels = 'A')
combined_plot
ggsave("Fig2.tiff",combined_plot, width = 8, height = 4,dpi = 600) 

#### SI figures #############
# just control
fitdensC<-fitted_early[fitted_early$treatment=="C",]
fullC<-full_data[full_data$treatment=="C",]

denC<-
ggplot(fitdensC, aes(x=day,y=density))+
  geom_line(data=fitdensC,stat='identity',color = "dodgerblue", size=3)+
  geom_line(data=fullC,(aes(x=day,y=density, group = rep)),color = "dodgerblue",alpha=0.2)+
  scale_y_continuous(limits = c(0,510),expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
denC

fitdensMB<-fitted_early[fitted_early$treatment=="MB",]
fullMB<-full_data[full_data$treatment=="MB",]

denMB<-
ggplot(fitdensMB, aes(x=day,y=density))+
  geom_line(data=fitdensMB,stat='identity',color = "red", size=3)+
  geom_line(data=fullMB,(aes(x=day,y=density, group = rep)),color = "red",alpha=0.2)+
  scale_y_continuous(limits = c(0,510),expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")


fitdensOP<-fitted_early[fitted_early$treatment=="OP",]
fullOP<-full_data[full_data$treatment=="OP",]

denOP<-
ggplot(fitdensOP, aes(x=day,y=density))+
  geom_line(data=fitdensOP,stat='identity',color = "darkorange", size=3)+
  geom_line(data=fullOP,(aes(x=day,y=density, group = rep)),color = "darkorange",alpha=0.2)+
  scale_y_continuous(limits = c(0,510),expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")

fitdensOPMB<-fitted_early[fitted_early$treatment=="OPMB",]
fullOPMB<-full_data[full_data$treatment=="OPMB",]

denOPMB<-
ggplot(fitdensOPMB, aes(x=day,y=density))+
  geom_line(data=fitdensOPMB,stat='identity',color = "green3", size=3)+
  geom_line(data=fullOPMB,(aes(x=day,y=density, group = rep)),color = "green3",alpha=0.2)+
  scale_y_continuous(limits = c(0,510),expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  geom_vline(xintercept = 30,size=3, color = "grey")+
  geom_vline(xintercept = 30,size=1, color = "black")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")

designSI <- "
A
B
C
D
"
combined_plotSI <- denC + denMB + denOP + denOPMB + plot_layout(design = designSI)+
  plot_annotation(tag_levels = 'A')
combined_plotSI
ggsave("Fig2SI.tiff",combined_plotSI, width = 8, height = 10,dpi = 600) 

########## transgenerational virulence analysis ###############################

#plotting early reproduction
earlyonlySE<-summarySE(data=data_first23, measurevar = "density",groupvars = c("treatment","day"))

ggplot(earlyonlySE, aes(x=day,y=density,group = treatment,color = treatment))+
  geom_point(position = position_dodge(1))+
  geom_errorbar(aes(ymin=density-se,ymax=density+se),width=0,position = position_dodge(1))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")

# Day 16 is the first day that we can test TGV for

dataTGV<-data_first23[data_first23$day==16,]

# merging MB and C together and OP and OPMB together because day 16 is before AM exposure
dataTGV$parasite <- with(dataTGV, ifelse(treatment == "C", "C",
                                  ifelse(treatment == "MB", "C", 
                                  ifelse(treatment == "OP", "OP",
                                  ifelse(treatment == "OPMB", "OP", NA)))))

ADearlyonlySE<-summarySE(data=dataTGV, measurevar = "adults",groupvars = c("parasite"))

ggplot(ADearlyonlySE, aes(x=parasite,y=adults,group = parasite,color = parasite))+
  geom_point(position = position_dodge(1))+
  geom_errorbar(aes(ymin=adults-se,ymax=adults+se),width=0,position = position_dodge(1))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")

JUearlyonlySE<-summarySE(data=dataTGV, measurevar = "juveniles",groupvars = c("parasite"))

ggplot(JUearlyonlySE, aes(x=parasite,y=juveniles,group = parasite,color = parasite))+
  geom_point(position = position_dodge(1))+
  geom_errorbar(aes(ymin=juveniles-se,ymax=juveniles+se),width=0,position = position_dodge(1))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")

# analyze day 16 to find TGV

TGVadultM1<-aov(adults~parasite,data = dataTGV)
anova(TGVadultM1)

op <- par(mfrow = c(2, 2))
plot(TGVadultM1)
par(mfrow = c(1, 1))

#calculating % reduction of adult density in OP vs C
(ADearlyonlySE$adults[ADearlyonlySE$parasite=="C"]-ADearlyonlySE$adults[ADearlyonlySE$parasite=="OP"])/ADearlyonlySE$adults[ADearlyonlySE$parasite=="C"]

# juveniles
TGVjuvenilesM1<-aov(juveniles~parasite,data = dataTGV)
anova(TGVjuvenilesM1)

op <- par(mfrow = c(2, 2))
plot(TGVjuvenilesM1)
par(mfrow = c(1, 1))

#calculating % reduction of juveniles density in OP vs C
(JUearlyonlySE$juveniles[JUearlyonlySE$parasite=="C"]-JUearlyonlySE$juveniles[JUearlyonlySE$parasite=="OP"])/JUearlyonlySE$juveniles[JUearlyonlySE$parasite=="C"]

###############################################################################
################ A. monospora prevalence ######################################
###############################################################################

# Create new dataframe with just data after day 31 when A. monospora was dosed
# and just treatments dosed with A. monospora
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
MBpred<-predict(PMbeta,newdata = MBprev,type = "link")

######## a bunch of code to generate standard errors for the fitted values from betareg

# Design matrix for the mean model
MBprevMat <- model.matrix(~ day * treatment, data = MBprev)

# Mean-model coefficients and their covariance matrix
b <- coef(PMbeta, model = "mean")
V <- vcov(PMbeta, model = "mean")

# Match columns to coefficient order
MBprevMat <- MBprevMat[, names(b), drop = FALSE]

# Predictions and standard errors on the logit scale
MBprev_fit_link <- drop(MBprevMat %*% b)
MBprev_se_link <- sqrt(rowSums((MBprevMat %*% V) * MBprevMat))

# Fitted means and approximate 95% confidence intervals
MBprev <- MBprev %>%
  mutate(
    predicted = plogis(MBprev_fit_link),
    se = MBprev_se_link * predicted * (1 - predicted),
    lower = predicted - se,
    upper = predicted + se
  )

####### Figure

MBplot<-
ggplot(MBprev, aes(x=day, y=predicted, group=treatment, color=treatment, shape = treatment)) +
  #geom_line(data=data_MBprev,(aes(x=day,y=MBprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5,fill = "grey", colour=NA) + 
  geom_line( aes(x=day, y=predicted), size=1) + 
  scale_colour_manual(values=c('red','green3'),breaks=c("MB","OPMB"),
                      labels=c("AM","OPAM")) +
  scale_x_continuous(limits=c(29.7, 94), expand=c(0, 0), name="Day") +
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0), name="A. monospora prevalence") +
  geom_rect(MBprev, mapping=aes(xmin=56, xmax=93, ymin=0.965, ymax=0.98), fill="red", color="red", alpha=0.5) +
  theme_classic()+ 
  theme(
    legend.position=c(0.15, 0.82),
    legend.background=element_blank(),
    axis.title=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank(),
    legend.text=element_text(size=14, color="black"),
    legend.key = element_rect(fill = scales::alpha("grey", 0.5), colour = NA)
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

# gam with separate smoothers
POPgam2<-gam(OPtotalprev~treatment + s(day, by=treatment),data = data_OPprev)
summary(POPgam2)
anova(POPgam2)

#model diagnostics
op <- par(mfrow = c(2, 2))
gam.check(POPgam2)
par(op)

# Post-hoc pairwise comparisons
emmsPOPgam <- emmeans(POPgam2, ~ treatment | day, 
          at = list(day = seq(min(data_OPprev$day), max(data_OPprev$day), by = 1)))
contrast(emmsPOPgam, method = "pairwise", adjust = "bonferroni")
emmip(POPgam2, treatment ~ day, at = list(day = seq(min(data_OPprev$day), max(data_OPprev$day), by = 1)), CIs = TRUE)

##### FIGURE ####
OPpred <- expand.grid(day = seq(min(data_OPprev$day), max(data_OPprev$day), length.out = 100),
                      treatment = levels(data_OPprev$treatment))

# create new column with densities predicted for each timepoint and SE using the GAM model
predictions_OP <- predict(POPgam2, newdata = OPpred, se.fit = TRUE)
OPpred$OPprev <- predictions_OP$fit
OPpred$se <- predictions_OP$se.fit

OPplot<-
ggplot(OPpred, aes(x=day,y=OPprev,group = treatment, color = treatment, shape = treatment))+
  #geom_line(data=data_OPprev,(aes(x=day,y=OPprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_smooth(data=OPpred,aes(ymin=OPprev-se, ymax=OPprev+se),stat='identity',alpha = .3, size=1)+
  scale_colour_manual(values = c( 'darkorange','green3'),breaks=c("OP","OPMB"),
                      labels=c("OP","OPAM"))+
  scale_x_continuous(limits=c(29.7, 94),expand = c(0,0), name = "Day")+
  scale_y_continuous(limits=c(NA, 1),expand = c(0,0),name = "O. pajunii prevalence")+
  geom_rect(data=data_OPprev, mapping=aes(xmin=37, xmax=46, ymin=0.965, ymax=0.98), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=data_OPprev, mapping=aes(xmin=63, xmax=93, ymin=0.965, ymax=0.98), fill="darkorange", color="darkorange", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.15,0.82),
        legend.background = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment", shape="Treatment")
OPplot
MBplot <- MBplot + theme(plot.tag = element_text(size = 28))
OPplot <- OPplot + theme(plot.tag = element_text(size = 28))
combined_plot2 <- MBplot + OPplot + plot_layout(heights = c(1, 1)) + 
  plot_annotation(tag_levels = 'A')
combined_plot2
ggsave("Fig3.tiff",combined_plot2, width = 5, height = 7, dpi=600) 


### SI figure ##############################################

MBprevMB<-MBprev[MBprev$treatment=="MB",]
data_MBprevMB<-data_MBprev[data_MBprev$treatment=="MB",]

MBplotMBSI<-
  ggplot(MBprevMB, aes(x=day, y=predicted, group=treatment, color=treatment, shape = treatment)) +
  geom_line(data=data_MBprevMB,(aes(x=day,y=MBprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5,fill = "grey", colour=NA) + 
  geom_line( aes(x=day, y=predicted), size=2) + 
  scale_colour_manual(values=c('red','green3'),breaks=c("MB","OPMB"),
                      labels=c("AM","OPAM")) +
  scale_x_continuous(limits=c(29.7, 94), expand=c(0, 0), name="Day") +
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0), name="A. monospora prevalence") +
  #geom_rect(MBprev, mapping=aes(xmin=56, xmax=93, ymin=0.965, ymax=0.98), fill="red", color="red", alpha=0.5) +
  theme_classic()+ 
  theme(
    legend.position=c(0.15, 0.82),
    legend.background=element_blank(),
    axis.title=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank(),
    legend.text=element_text(size=14, color="black"),
    legend.key = element_rect(fill = scales::alpha("grey", 0.5), colour = NA)
  ) +
  labs(color="Treatment", shape="Treatment") 
MBplotMBSI

MBprevOP<-MBprev[MBprev$treatment=="OPMB",]
data_MBprevOP<-data_MBprev[data_MBprev$treatment=="OPMB",]

MBplotOPSI<-
  ggplot(MBprevOP, aes(x=day, y=predicted, group=treatment, color=treatment, shape = treatment)) +
  geom_line(data=data_MBprevOP,(aes(x=day,y=MBprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5,fill = "grey", colour=NA) + 
  geom_line( aes(x=day, y=predicted), size=2) + 
  scale_colour_manual(values=c('red','green3'),breaks=c("MB","OPMB"),
                      labels=c("AM","OPAM")) +
  scale_x_continuous(limits=c(29.7, 94), expand=c(0, 0), name="Day") +
  scale_y_continuous(limits=c(0, 1), expand=c(0, 0), name="A. monospora prevalence") +
  #geom_rect(MBprev, mapping=aes(xmin=56, xmax=93, ymin=0.965, ymax=0.98), fill="red", color="red", alpha=0.5) +
  theme_classic()+ 
  theme(
    legend.position=c(0.15, 0.82),
    legend.background=element_blank(),
    axis.title=element_text(size=14, color="black"),
    axis.text=element_text(size=14, color="black"),
    legend.title=element_blank(),
    legend.text=element_text(size=14, color="black"),
    legend.key = element_rect(fill = scales::alpha("grey", 0.5), colour = NA)
  ) +
  labs(color="Treatment", shape="Treatment") 
MBplotOPSI

OPpredOP<-OPpred[OPpred$treatment=="OP",]
data_OPprevOP<-data_OPprev[data_OPprev$treatment=="OP",]

OPplotOPSI<-
  ggplot(OPpredOP, aes(x=day,y=OPprev,group = treatment, color = treatment, shape = treatment))+
  geom_line(data=data_OPprevOP,(aes(x=day,y=OPprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_smooth(data=OPpredOP,aes(ymin=OPprev-se, ymax=OPprev+se),stat='identity',alpha = .3, size=2)+
  scale_colour_manual(values = c( 'darkorange','green3'),breaks=c("OP","OPMB"),
                      labels=c("OP","OPAM"))+
  scale_x_continuous(limits=c(29.7, 94),expand = c(0,0), name = "Day")+
  scale_y_continuous(limits=c(NA, 1),expand = c(0,0),name = "O. pajunii prevalence")+
  #geom_rect(data=data_OPprev, mapping=aes(xmin=44, xmax=93, ymin=0.965, ymax=0.98), fill="darkorange", color="darkorange", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.15,0.82),
        legend.background = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment", shape="Treatment")
OPplotOPSI

OPpredOPMB<-OPpred[OPpred$treatment=="OPMB",]
data_OPprevOPMB<-data_OPprev[data_OPprev$treatment=="OPMB",]

OPplotOPMBSI<-
  ggplot(OPpredOPMB, aes(x=day,y=OPprev,group = treatment, color = treatment, shape = treatment))+
  geom_line(data=data_OPprevOPMB,(aes(x=day,y=OPprev, group = interaction(treatment,rep))),alpha=0.2)+
  geom_smooth(data=OPpredOPMB,aes(ymin=OPprev-se, ymax=OPprev+se),stat='identity',alpha = .3, size=2)+
  scale_colour_manual(values = c( 'darkorange','green3'),breaks=c("OP","OPMB"),
                      labels=c("OP","OPAM"))+
  scale_x_continuous(limits=c(29.7, 94),expand = c(0,0), name = "Day")+
  scale_y_continuous(limits=c(NA, 1),expand = c(0,0),name = "O. pajunii prevalence")+
  #geom_rect(data=data_OPprev, mapping=aes(xmin=44, xmax=93, ymin=0.965, ymax=0.98), fill="darkorange", color="darkorange", alpha=0.5) +
  theme_classic()+
  theme(legend.position = c(0.15,0.82),
        legend.background = element_blank(),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment", shape="Treatment")
OPplotOPMBSI


combined_plot2SI <- MBplotMBSI + MBplotOPSI +  OPplotOPSI+ OPplotOPMBSI + 
  plot_annotation(tag_levels = 'A')
combined_plot2SI
ggsave("Fig3SI.tiff",combined_plot2SI, width = 9, height = 7, dpi=600) 

###########################################################################
########### DEMOGRAPHICS ##################################################
###########################################################################

#new dataframe to use for demographics analysis, includes density of adults and juveniles
demo<-full_data[c("treatment","rep","day","adults","juveniles")]
demo$juvrat<-demo$juveniles/(demo$juveniles+demo$adults) # calculate juvenile ratio in populations
demo_past23<-demo[demo$day>23,] #remove days <30
demo_past23 <- demo_past23[!is.nan(demo_past23$juvrat), ]
demo_past23$repuniqID <- factor(paste(demo_past23$treatment, demo_past23$rep, sep = "_"))

#### GAM model for demography ####
juvgam2<-gam(juvrat~treatment+s(day, by=treatment),data = demo_past23)

anova(juvgam2)
summary(juvgam2)
concurvity(juvgam2, full = TRUE)
plot(juvgam2, pages = 1, residuals = TRUE)
vis.gam(juvgam2, theta = 220, color = "heat") 
# model diagnostics
op <- par(mfrow = c(2, 2))
gam.check(juvgam2)
par(mfrow = c(1, 1))

res_juv_gam <- resid(juvgam2)
coplot(res_juv_gam ~ day | treatment, data = demo_past23)

days_common <- seq(30, 93, by = 1)

EM_juv_average <- emmeans(juvgam2,~ treatment,at = list(day = days_common),weights = "equal")
pairs(EM_juv_average, adjust = "bonferroni")

############################# contrasts
#contrasts for the main effect
EM_juv_bytreat <- emmeans(juvgam2,~ treatment,at = list(day = sort(unique(demo_past23$day))),weights = "equal")
pairs(EM_juv_bytreat, adjust = "tukey")

#contrasts by day
EM_juvrat_gam <- emmeans(juvgam2, ~ treatment | day, at = list(day = seq(min(demo_past23$day), max(demo_past23$day), by = 1)))
contrast(EM_juvrat_gam, method = "pairwise", adjust = "bonferroni") # pairwise contrast comparison with bonferroni correction
contrasts_juv_gam<-contrast(EM_juvrat_gam, method = "pairwise", adjust = "bonferroni")

# Exporting post-hoc comparisons into .csv files
pairwise_comparisons_juv_gam <- as.data.frame(summary(contrasts_juv_gam))
write.csv(pairwise_comparisons_juv_gam, file = "pairwise_comparisons_juv_gam.csv", row.names = FALSE)
write.csv(EM_juvrat_gam, file = "EM_juv_gam.csv", row.names = FALSE)

# Get predictions from the model
JRM1<- expand.grid(day = seq(min(demo_past23$day), max(demo_past23$day), length.out = 100),
                   treatment = levels(demo_past23$treatment))

predictions_juvrat <- predict(juvgam2, newdata = JRM1, se.fit = TRUE)
JRM1$juvrat <- predictions_juvrat$fit
JRM1$se <- predictions_juvrat$se.fit

JRplot<-
ggplot(JRM1, aes(x = day, y = juvrat, group=treatment,fill = "grey")) +
  geom_ribbon(aes(ymin = juvrat - se, ymax = juvrat + se, fill = "grey"), alpha = 0.2) +
  geom_line(aes(y = juvrat, color = treatment), linewidth = 1) +
  ylab(label = "proportion of juveniles") +
  scale_color_manual(values = c('C' = 'dodgerblue', 'MB' = 'red', 'OP' = 'darkorange', 'OPMB' = 'green3'),breaks=c("C","MB","OP","OPMB"),
                     labels=c("C","AM","OP","OPAM")) +
  scale_fill_manual(values = c('C' = 'dodgerblue', 'MB' = 'red', 'OP' = 'darkorange', 'OPMB' = 'green3'),breaks=c("C","MB","OP","OPMB"),
                     labels=c("C","AM","OP","OPAM")) +
  geom_rect(data=JRM1, mapping=aes(xmin=30, xmax=31, ymin=0.01, ymax=0.03), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=JRM1, mapping=aes(xmin=59.9, xmax=60.1, ymin=0.01, ymax=0.03), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=JRM1, mapping=aes(xmin=90, xmax=91, ymin=0.01, ymax=0.03), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=JRM1, mapping=aes(xmin=31, xmax=32, ymin=0.034, ymax=0.054), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=JRM1, mapping=aes(xmin=86, xmax=93, ymin=0.034, ymax=0.054), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=JRM1, mapping=aes(xmin=79, xmax=81, ymin=0.98, ymax=1), fill="red", color="red", alpha=0.5) +
  scale_x_continuous(expand = c(0,0), name = "Day")+
  #geom_point(data=demo_past23, aes(x=day,y=juvrat,color=treatment),alpha=.4)+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
   theme_classic()+
  theme(legend.position = c(0.91,0.79),
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")
JRplot
ggsave("Fig4.tiff",JRplot, width = 6, height = 4, dpi=600) 

####### Figure SI #######################################

# just control
fitjuvC<-JRM1[JRM1$treatment=="C",]
demoC<-demo_past23[demo_past23$treatment=="C",]

juvC<-
  ggplot(fitjuvC, aes(x=day,y=juvrat))+
  geom_line(data=fitjuvC,stat='identity',color = "dodgerblue", size=3)+
  geom_line(data=demoC,(aes(x=day,y=juvrat, group = rep)),color = "dodgerblue",alpha=0.2)+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
juvC

fitjuvAM<-JRM1[JRM1$treatment=="MB",]
demoAM<-demo_past23[demo_past23$treatment=="MB",]

juvAM<-
  ggplot(fitjuvAM, aes(x=day,y=juvrat))+
  geom_line(data=fitjuvAM,stat='identity',color = "red", size=3)+
  geom_line(data=demoAM,(aes(x=day,y=juvrat, group = rep)),color = "red",alpha=0.2)+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
juvAM



fitjuvOP<-JRM1[JRM1$treatment=="OP",]
demoOP<-demo_past23[demo_past23$treatment=="OP",]

juvOP<-
  ggplot(fitjuvOP, aes(x=day,y=juvrat))+
  geom_line(data=fitjuvOP,stat='identity',color = "darkorange", size=3)+
  geom_line(data=demoOP,(aes(x=day,y=juvrat, group = rep)),color = "darkorange",alpha=0.2)+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
juvOP

fitjuvOPAM<-JRM1[JRM1$treatment=="OPMB",]
demoOPAM<-demo_past23[demo_past23$treatment=="OPMB",]

juvOPAM<-
  ggplot(fitjuvOPAM, aes(x=day,y=juvrat))+
  geom_line(data=fitjuvOPAM,stat='identity',color = "green3", size=3)+
  geom_line(data=demoOPAM,(aes(x=day,y=juvrat, group = rep)),color = "green3",alpha=0.2)+
  scale_y_continuous(limits=c(0, 1),expand = c(0,0),name = "Proportion of juveniles")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = c(0.88,0.22),
        legend.background = element_blank(),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 11, color = "black"),
        legend.text = element_text(size = 11, color = "black"))+
  labs(color="Treatment")
juvOPAM

designSIj <- "
A
B
C
D
"

combined_plot_juvSI <- juvC + juvAM + juvOP + juvOPAM +
  plot_layout(design = designSIj) +
  plot_annotation(tag_levels = "A")&
  theme(plot.tag = element_text(size = 14,face = "bold",family = "sans"))

combined_plot_juvSI

ggsave("Fig4SI.tiff",plot = combined_plot_juvSI,width = 8,height = 12,units = "in",dpi = 600)


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
gam_mod <- gam(algae ~  treatment + s(day, by=treatment),family = tw(link = "log"), data = algaedatNo0)
summary(gam_mod)
anova(gam_mod)
plot(gam_mod$residuals)

op <- par(mfrow = c(2, 2))
gam.check(gam_mod)
par(mfrow = c(1, 1))

concurvity(gam_mod,full=T)

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
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  geom_rect(data=fitalg, mapping=aes(xmin=25, xmax=64, ymin=16.8, ymax=17), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=88, xmax=93, ymin=16.8, ymax=17), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=77, xmax=79, ymin=16.8, ymax=17), fill="dodgerblue", color="green3", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=28, xmax=56, ymin=16.4, ymax=16.6), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=92, xmax=93, ymin=16.4, ymax=16.6), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=67, xmax=84, ymin=16.4, ymax=16.6), fill="dodgerblue", color="darkorange", alpha=0.5) +
  geom_rect(data=fitalg, mapping=aes(xmin=66, xmax=82, ymin=16, ymax=16.2), fill="dodgerblue", color="red", alpha=0.5) +
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
ggsave("FigS2.tiff",algGAMplot, width = 8, height = 6,dpi = 600) 
