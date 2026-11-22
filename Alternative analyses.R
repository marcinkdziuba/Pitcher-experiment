###############################################################################
################# ALTERNATIVE ANALYSES SCRIPT #################################
###############################################################################



###############################################################################
################ INDIVIDUAL LEVEL EXP #########################################
###############################################################################

############ REPRODUCTION #####################################################

# linear model
MGvM_LM<-lm(NB~parasite*clutch,data = repMGvM)
anova(MGvM_LM)
op <- par(mfrow = c(2, 2))
plot(MGvM_LM)
par(mfrow = c(1, 1))

# gam model with common smoother for all treatments
MGvM_M1<-gam(NB~parasite+s(clutch),data = repMGvM)
anova(MGvM_M1) #stats report
summary(MGvM_M1)

# model diagnostics
op <- par(mfrow = c(2, 2))
gam.check(MGvM_M1)
par(mfrow = c(1, 1))
vis.gam(MGvM_M1, theta = 150, color = "heat") # visualisation of gam model
concurvity(MGvM_M1, full = TRUE)

# The pattern in response vs Fitted is due to smaller variation in number of 
# offspring produced in first few clutches, when the clutches are small 

AIC(MGvM_LM,MGvM_M1,MGvM_M2) # lower AIC for the more complex model, we'll use MGvM_M2 then

#### SURVIVAL #################################################################

# linear model
MortM1<-lm(lifespan~parasite, data=mortMGvM)
anova(MortM1)
op <- par(mfrow = c(2, 2))
plot(MortM1)
par(mfrow = c(1, 1))

# Kruskal-Wallis test
KW_LSM<-kruskal.test(lifespan~parasite, data=mortMGvM)
KW_LSM
# Post-hoc
dunn.test(mortMGvM$lifespan,mortMGvM$parasite, method = "bonferroni")

############ AGE AT REPRODUCTION ###############################

#linear model
LMage<-lm(age~parasite*as.factor(clutch),data=repMGvM_4C)
anova(LMage)
op <- par(mfrow = c(2, 2))
plot(LMage)
par(mfrow = c(1, 1))

# mixed effects
M_lme <- lme(age~parasite, random = ~1|clutch,  data = repMGvM_4C)
anova(M_lme)
op <- par(mfrow = c(2, 2))
plot(M_lme)
par(mfrow = c(1, 1))

AIC(M_lme,M_lme_cv)


###############################################################################
################ POPULATION EXPERIMENT ########################################
###############################################################################

########### POPULATION GROWTH RATE ############################################

# Growth rate comparison between groups
# attemnt on using a linear model
GRm1<-lm(r~parasite, data = data_GRwide)
anova(GRm1)

op <- par(mfrow = c(2, 2))
plot(GRm1) 
par(mfrow = c(1, 1))
# linear model after removing the obvious outlier #39
data_GRwide2<-data_GRwide[-39,]
GRm2<-lm(r~parasite, data = data_GRwide2)
anova(GRm2)

op <- par(mfrow = c(2, 2))
plot(GRm2) 
par(mfrow = c(1, 1))# non-normal distribution of residuals

########### POPULATION DENSITY ################################################

# random eff
dens_gam_r<-gam(density~treatment + s(day, by=treatment)+ s(pitcherID,bs = "re"),data=data_past23)
anova(dens_gam_r)
concurvity(dens_gam_r)
op <- par(mfrow = c(2, 2))
gam.check(dens_gam_r)
par(mfrow = c(1, 1))
res_dens_gam_r <- resid(dens_gam_r)
coplot(res_dens_gam_r ~ day | treatment, data = data_past23)

AIC(dens_gam,dens_gam_r)

dens_gam_r2<-gam(density~treatment + s(day, by=treatment)+ s(pitcherID,bs = "re")+s(pitcherID, by=day, bs="re"),data=data_past23)
anova(dens_gam_r2)
concurvity(dens_gam_r2)
op <- par(mfrow = c(2, 2))
gam.check(dens_gam_r2)
par(mfrow = c(1, 1))
res_dens_gam_r2 <- resid(dens_gam_r2)
coplot(res_dens_gam_r2 ~ day | treatment, data = data_past23)

AIC(dens_gam,dens_gam_r,dens_gam_r2)

# plot mean and errorbar for each treatment x timepoint interaction based on raw data

densSErev<- summarySE(data=data_past23, measurevar = "density", groupvars=c("treatment","day"),na.rm=T)

ggplot(densSErev, aes(x=day,y=density,group = treatment,color = treatment))+
  geom_point(position=position_dodge(2))+
  geom_errorbar(aes(ymin=density-se,ymax=density+se),width=0,position=position_dodge(2))+
  scale_colour_manual(values = c('dodgerblue', 'red','darkorange','green3'),breaks=c("C","MB","OP","OPMB"),
                      labels=c("C","AM","OP","OPAM"))+
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=64, ymin=5, ymax=10), fill="green3", color="green3", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=30, xmax=61, ymin=13, ymax=18), fill="darkorange", color="darkorange", alpha=0.5) +
  geom_rect(data=fitted_23, mapping=aes(xmin=71, xmax=83, ymin=320, ymax=325), fill="red", color="red", alpha=0.5) +
  #geom_line(data=data_past23,(aes(x=day,y=density,color=treatment, group = interaction(treatment, rep))),alpha=0.2,position = position_dodge(1))+
  scale_y_continuous(expand = c(0,0), name = "Daphnia density/L")+
  scale_x_continuous(breaks=c(0,30,60,90),expand = c(0,0), name = "Day")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title = element_text(size = 14, color = "black"),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"))+
  labs(color="Treatment")

################ A. monospora prevalence ######################################

# Binomial data
data_MBprev$MBpositive<-rowSums(data_MBprev[, c("MBA", "OPMBA", "MBJ", "OPMBJ", "MBM", "OPMBM")],na.rm = TRUE)
data_MBprev$MBnegative <- rowSums(data_MBprev[, c("UA", "OPA", "UJ", "OPJ", "UM", "OPM")],na.rm = TRUE)

MBprevBinom <- gam(cbind(MBpositive, MBnegative) ~ treatment + s(day, by=treatment, k=8),
                   family = binomial(link = "logit"),data = data_MBprev)
anova(MBprevBinom)

op <- par(mfrow = c(2, 2))
gam.check(MBprevBinom)
par(op)

data_MBprev$pitcherID <- interaction(data_MBprev$treatment,data_MBprev$rep,drop = TRUE,sep = "_")

MBprevBinom2 <- gam(cbind(MBpositive, MBnegative) ~ treatment + s(day, by=treatment, k=8)+ s(pitcherID, bs = "re"),family = binomial(link = "logit"),
                    data = data_MBprev)

anova(MBprevBinom2)

op <- par(mfrow = c(2, 2))
gam.check(MBprevBinom2)
par(op)

MBprevBin_pred <- expand.grid(day = seq(min(data_MBprev$day),max(data_MBprev$day),length.out = 100),
                            treatment = c("MB", "OPMB"))

 predictions_MB <- predict(MBprevBinom,newdata = MBprevBin_pred,type = "link",se.fit = TRUE)

 MBprevBin_pred$predicted <- plogis(as.vector(predictions_MB$fit))
 MBprevBin_pred$lower <- plogis(as.vector(predictions_MB$fit - predictions_MB$se.fit))
 MBprevBin_pred$upper <- plogis(as.vector(predictions_MB$fit + predictions_MB$se.fit))

 ggplot(MBprevBin_pred, aes(x=day, y=predicted, group=treatment, color=treatment, shape = treatment)) +
   #geom_line(data=data_MBprev,(aes(x=day,y=MBprev, group = interaction(treatment,rep))),alpha=0.2)+
   geom_ribbon(aes(ymin = lower, ymax = upper),alpha = 0.5,fill = "grey", colour=NA) +
   geom_line( aes(x=day, y=predicted), size=1) +
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

 ######################### O. pajunii prevalence ################################ 
 
 #gam model with single smoother for both groups
 POPgam<-gam(OPtotalprev~treatment + s(day),data = data_OPprev)
 summary(POPgam)
 anova(POPgam)
 
 #model diagnostics
 op <- par(mfrow = c(2, 2))
 gam.check(POPgam)
 par(op)
 
 AIC(POPgam, POPgam2)
 
 ########### DEMOGRAPHICS ##################################################
 
 # linear model of treatment*day interaction
 juvratM1<-lm(juvrat~treatment*day,data=demo_past23)
 summary(juvratM1)
 anova(juvratM1)
 par(mfrow=c(2,2))
 plot(juvratM1)
 par(mfrow=c(1,1))
 
 emmjuvrat<-emmeans(juvratM1,~treatment)
 contrast(emmjuvrat, method = "pairwise", adjust = "bonferroni")
 
# gam model with common smoother 
 juvgam<-gam(juvrat~treatment+s(day), data = demo_past23)
 anova(juvgam)
 
 AIC(juvratM1,juvgam)
 anova(juvratM1,juvgam)
 concurvity(juvgam, full = TRUE)
 op <- par(mfrow = c(2, 2))
 gam.check(juvgam)
 par(mfrow = c(1, 1))
 
 AIC(juvratM1,juvgam,juvgam2)
 anova(juvratM1,juvgam,juvgam2)
 
 
 
 juvgam3<-gam(juvrat~treatment+s(day, by=treatment)+
                s(repuniqID, bs = "re"),data = demo_past23)
 
 anova(juvgam3)
 concurvity(juvgam3, full = TRUE)
 plot(juvgam3, pages = 1, residuals = TRUE)
 summary(juvgam3)
 vis.gam(juvgam3, theta = 220, color = "heat")
 # model diagnostics
 op <- par(mfrow = c(2, 2))
 gam.check(juvgam3)
 par(mfrow = c(1, 1))
 
 juvgam4 <-gam(juvrat ~ treatment +
                 s(day, by = treatment) +
                 s(repuniqID, bs = "re") +
                 s(repuniqID, by = day, bs = "re"),
               data = demo_past23)
 
 anova(juvgam4)
 concurvity(juvgam4, full = TRUE)
 plot(juvgam4, pages = 1, residuals = TRUE)
 summary(juvgam4)
 vis.gam(juvgam4, theta = 220, color = "heat")
 # model diagnostics
 op <- par(mfrow = c(2, 2))
 gam.check(juvgam4)
 par(mfrow = c(1, 1))
 
 AIC(juvgam2,juvgam3,juvgam4)
 