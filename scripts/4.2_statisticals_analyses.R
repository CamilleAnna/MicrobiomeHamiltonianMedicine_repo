# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#             Using evolutionary theory to predict microbes’ effect on host health
#                           Camille Simonet & Luke McNally
#         
#                            SCRIPT: Statistical analyses
#                               last edit: 21/08/2021
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#local_project_dir="path/to/repo"
#local_project_dir="/Users/s1687811/Documents/PhD/Research/MicrobiomeHamiltonianMedicine/MicrobiomeHamiltonianMedicine_repo"



# SETUP ----
# Sourcing packages ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(local_project_dir)
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(vegan)
library(lme4)
library(lmerTest)
library(irr)
library(effects)
library(jtools)
library(ROCR)
library(ggrepel)

# Load environment preped by script 4.1_data_prep.R 
load('./output/stats_runs/data_preps_23062021.RData')

# Data overview
table(dat.rescaled$caseControl)
nrow(dat.rescaled)


# MICROBIOME ANALYSIS ----

# Binomial regressions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We put disease as random effect allowing both intercepts and slope to vary per disease
# Let the model estimate ful;l covariance matrix
m.glmm.disease <- glmer(caseControl ~ 1+ NR + NRSPO + (1 + NR + NRSPO|disease_type),
                        family = binomial(link = "logit"),
                        data = dat.rescaled,
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=5e5)))

summary(m.glmm.disease)


# Not enough cohort per disease to distinguish the two effect, so partly confounded
# re-run analysis with cohort as random effect
# find similar result
m.glmm.cohort <- glmer(caseControl ~ 1+ NR + NRSPO + (1 + NR + NRSPO|cohortNew),
                       family = binomial(link = "logit"),
                       data = dat.rescaled,
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5)))

summary(m.glmm.cohort)


# Test random slopes (LRT) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Re-run model with the random slopes
# use LRT tests to verify if indeed allowing random slopes give a better fit

m.glmm.disease.fixedSlopes <- glmer(caseControl ~ 1+ NR + NRSPO + (1 |disease_type),
                                    family = binomial(link = "logit"),
                                    data = dat.rescaled,
                                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=5e5)))

summary(m.glmm.disease.fixedSlopes) # same results for NR, NRSPO comes out significant too



m.glmm.cohort.fixedSlopes <- glmer(caseControl ~ 1+ NR + NRSPO + (1 |cohortNew),
                                   family = binomial(link = "logit"),
                                   data = dat.rescaled,
                                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=5e5)))

summary(m.glmm.cohort.fixedSlopes) # same results


anova(m.glmm.disease.fixedSlopes, m.glmm.disease) # model with random slopes gives sig. better fit
anova(m.glmm.cohort.fixedSlopes, m.glmm.cohort)   # model with random slopes gives sig. better fit


# Check with MCMC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Verify with MCMCglmm if I get the same result for the choosen model (per disease, with random slopes)

prior.0<-list(R=list(V=1, fix=1),
              G=list(G1=list(V=diag(3), nu=3, alpha.mu=c(0,0,0), alpha.V=diag(3)*1000))) # variance in fixed effect across diseases

mcmc.disease<- MCMCglmm(caseControl ~ 1+ NR+NRSPO,
                        random = ~us(1+NR+NRSPO):disease_type,
                        data = dat.rescaled,
                        family = 'threshold',
                        prior = prior.0,
                        nitt = 150000, burnin = 50000, thin = 100)

summary(mcmc.disease) # MCMCglmm gives similar results



# MAKE SUPPLEMENTARY TABLES ----


source('./scripts/5.1_suppMat_make_modelSummaries.R')



# PLOT MICROBIOME MODEL ----

# Boxplots NR/NRSPO ----


mygrey<- colorRampPalette(c('white', 'lightgrey'))(10)[4]

theme.bp<-
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5, colour = 'black'),
        panel.grid = element_blank(),
        #legend.position = 'top', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4))


dat.rescaled$cohortNewPlot<- factor(dat.rescaled$cohortNew,
                                    levels = (c('YeZ2018_BD', 'LiJ2017_HYPER', 'Zhang2015_RA','JieZ2017_ACVD', 'Karlsson2012_ACVD',
                                                'QinN2014_LC', 'QinJ2012_T2D', 'Karlsson2013_T2D&IGT', 
                                                'FengQ2015_CRC', 'YuJ2015_CRC', 'ZellerG2014_CRC', 'ThomasAM2018_CRC', 'VogtmannE2016_CRC',
                                                'He2017_IBD', 'NielsenH2014_IBD', 'iHMP_IBD', 'Costea2017_MetS', 'LeChatelierE2013_OBE', 'Liu2017_OBE')))

dat.rescaled$cohortNewPlot2<- do.call('rbind', strsplit(dat.rescaled$cohortNew, '_'))[,1]


dat.rescaled$cohortNewPlot2<- factor(dat.rescaled$cohortNewPlot2,
                                     levels = (c('YeZ2018', 'LiJ2017', 'Zhang2015','JieZ2017', 'Karlsson2012',
                                                 'QinN2014', 'QinJ2012', 'Karlsson2013', 
                                                 'FengQ2015', 'YuJ2015', 'ZellerG2014', 'ThomasAM2018', 'VogtmannE2016',
                                                'He2017', 'NielsenH2014', 'iHMP', 'Costea2017', 'LeChatelierE2013', 'Liu2017')))

# Add samples sizes under cohort names on y-axis ticks
nsizes<- dat.rescaled %>%
  group_by(cohortNewPlot2, caseControl) %>%
  summarise(n = n()) %>%
  as.data.frame() %>%
  ungroup() %>%
  spread(key = caseControl, value = n) %>%
  as.data.frame() %>%
  rename(case = `0`,
         control = `1`)
  

dat.rescaled<- left_join(dat.rescaled, nsizes, by = 'cohortNewPlot2')

dat.rescaled<- dat.rescaled %>%
  mutate(cohortNewPlot3 = paste0(cohortNewPlot2, '\n(c:', case, ', h:', control, ')'))
  

dat.rescaled$cohortNewPlot3<- factor(dat.rescaled$cohortNewPlot3,
                                     levels = rev(c("Liu2017\n(c:104, h:101)",
                                                "LeChatelierE2013\n(c:130, h:79)",
                                                "Costea2017\n(c:60, h:26)",
                                                "iHMP\n(c:89, h:26)",
                                                "NielsenH2014\n(c:60, h:41)",
                                                "He2017\n(c:47, h:41)",
                                                "VogtmannE2016\n(c:24, h:19)",
                                                "ThomasAM2018\n(c:84, h:45)",
                                                "ZellerG2014\n(c:111, h:45)",
                                                "YuJ2015\n(c:74, h:53)",
                                                "FengQ2015\n(c:93, h:43)",
                                                "Karlsson2013\n(c:98, h:33)",
                                                "QinJ2012\n(c:182, h:182)",
                                                "QinN2014\n(c:118, h:111)",
                                                "Karlsson2012\n(c:12, h:10)",
                                                "JieZ2017\n(c:202, h:159)" ,
                                                "Zhang2015\n(c:88, h:89)",
                                                "LiJ2017\n(c:155, h:41)",
                                                "YeZ2018\n(c:19, h:40)")))


p.NR.boxes<- ggplot(dat.rescaled, aes(x = cohortNewPlot3, y = NR, fill = as.factor(caseControl)))+
  ylab('Gut average relatedness')+xlab('Cohort')+
  scale_y_continuous(limits = c(0.35, 1), expand = c(0, 0))+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  #scale_fill_manual(values = c('lightgrey', 'white'))+
  scale_fill_manual(values = c('firebrick', 'dodgerblue'))+
  coord_flip()+
  theme.bp+
  theme(legend.position = 'none')


ymin.rect.NR = 0.35
ymax.rect.NR = 1
text.pos = 0.38
p.NR<- p.NR.boxes +
 annotate("rect", xmin = 17.5, xmax = 19.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR,fill = mygrey)+
  annotate("rect", xmin = 16.5, xmax = 17.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR,fill = 'white')+
  annotate("rect", xmin = 14.5, xmax = 16.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 8.5, xmax = 13.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = 'white')+
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = 'white')+
  annotate("rect", xmin = 3.5, xmax = 5.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = 'white')+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = mygrey)+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = ymin.rect.NR, ymax = ymax.rect.NR, fill = 'white')+
  
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  annotate('text', x = 18.5, y = text.pos, label = 'Obesity', angle = 90, size = 2)+
  annotate('text', x = 17, y = text.pos, label = 'MetS', angle = 90, size = 2)+
  #annotate('text', x = 15.5, y = text.pos, label = "CD", angle = 90, size = 2)+
  #annotate('text', x = 14, y = text.pos, label = 'UC', angle = 90, size = 2)+
  annotate('text', x = 15, y = text.pos, label = 'IBD', angle = 90, size = 2)+
  annotate('text', x = 11, y = text.pos, label = 'Colorectal Cancer', angle = 90, size = 2)+
  #annotate('text', x = 8, y = text.pos, label = 'IGT', angle = 90, size = 2)+
  #annotate('text', x = 7, y = text.pos, label = 'T2D', angle = 90, size = 2)+
  annotate('text', x = 7.5, y = text.pos, label = 'T2D & IGT', angle = 90, size = 2)+
  annotate('text', x = 6, y = text.pos, label = 'LC', angle = 90, size = 2)+
  annotate('text', x = 4.5, y = text.pos, label = 'ACVD', angle = 90, size = 2)+
  annotate('text', x = 3, y = text.pos, label = 'RA', angle = 90, size = 2)+
  annotate('text', x = 2, y = text.pos, label = 'HypT', angle = 90, size = 2)+
  annotate('text', x = 1, y = text.pos, label = 'BD', angle = 90, size = 2)

p.NR


p.NRSPO.boxes<- ggplot(dat.rescaled, aes(x = cohortNewPlot3, y = NRSPO, fill = as.factor(caseControl)))+
  xlab(' \ ')+ylab('Relatedness weighed mean\nsporulation score')+
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0))+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)+
  #scale_fill_manual(values = c('lightgrey', 'white'))+
  scale_fill_manual(values = c('firebrick', 'dodgerblue'))+
  coord_flip()+
  theme.bp+
  theme(legend.position = 'none',
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())


ymin.rect.NRSPO = 0
ymax.rect.NRSPO = 0.5
p.NRSPO<- p.NRSPO.boxes +
annotate("rect", xmin = 17.5, xmax = 19.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 16.5, xmax = 17.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = 'white')+
  annotate("rect", xmin = 14.5, xmax = 16.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 8.5, xmax = 13.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = 'white')+
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = 'white')+
  annotate("rect", xmin = 3.5, xmax = 5.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = 'white')+
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = mygrey)+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = ymin.rect.NRSPO, ymax = ymax.rect.NRSPO, fill = 'white')+
  geom_boxplot(width = 0.5, outlier.size = 0.05, lwd = 0.2)


# Effect plots ----

# *** To run this part need to run the models in script 4.2 ***

fm1 <- m.glmm.disease # store whatever model I want to use in fm1 variable and use this in code below

newdat <- expand.grid( # Make a new dataframe for predictions
  NR=seq(0.38, 0.94, 0.01),            # Make range of NR within range of observed values
  NRSPO = mean(dat.rescaled$NRSPO), # condition on typical (i.e. mean) value of NRSPO
  caseControl = NA
)

newdat$caseControl <- predict(fm1,newdat,re.form=NA) # Make predictions

plot(newdat$caseControl~ newdat$NR, type = 'l', ylim = c(-4, 1)) # Ok, what I get corresponds to what allEffect() wrapper gives
plot(allEffects(fm1), type = 'link')


# Ok now doing this manually (following Ben Bolkers tutorial on it) to customise plot
mm <- model.matrix(terms(fm1),newdat)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm)) # Thus is to produce CI based on fixed-effects uncertainty ONLY
tvar1 <- pvar1+VarCorr(fm1)$disease_type[1]+VarCorr(fm1)$disease_type[5]+VarCorr(fm1)$disease_type[9]  # This adds the variance of random effects. Adapted from Ben Bolker's code but here for more complex model. Because we let the model estimate variance in intercept, NR slope and NRSPO slope with disease_type, we must add those three sources of variance.
cmult <- 1.96 ## could use 1.96
newdat <- data.frame(
  newdat
  , plo = newdat$caseControl-cmult*sqrt(pvar1)
  , phi = newdat$caseControl+cmult*sqrt(pvar1)
  , tlo = newdat$caseControl-cmult*sqrt(tvar1)
  , thi = newdat$caseControl+cmult*sqrt(tvar1)
)
#plot confidence
g0 <- ggplot(newdat, aes(x=NR, y=plogis(caseControl)))+geom_point()
g0 + geom_pointrange(aes(ymin = plogis(plo), ymax = plogis(phi)))+
  labs(title="CI based on fixed-effects uncertainty ONLY")
#plot prediction
g0 + geom_pointrange(aes(ymin = plogis(tlo), ymax = plogis(thi)))+
  labs(title="CI based on FE uncertainty + RE variance")

plot(allEffects(fm1), type = 'response') # ok, what I get manually looks similar to what allEffects() wrapper gives



funFormat<- function(x, digits){ifelse(x<0.01,
                              formatC(x, digit = digits, format = 'e'),
                              formatC(x, digit = digits, format = 'f'))}

out.fix.format<- summary(m.glmm.disease)[[10]] %>%
  as.data.frame() %>%
  add_rownames('Effect') %>%
  as.data.frame() %>%
  mutate_if(is.numeric, funs(funFormat(., 2)))

b1 = out.fix.format[2,2]
pvalb1 = out.fix.format[2,5]


p.predict.NR<- ggplot(newdat, aes(x=NR, y=plogis(caseControl)))+
  #geom_segment(data = dat.rescaled, aes(x = NR, xend = NR, y = 0, yend = 0.02))+
  geom_segment(data = dat.rescaled[dat.rescaled$caseControl == '0',], aes(x = NR, xend = NR, y = 0, yend = 0.02), alpha = .2, col = 'firebrick')+
  geom_segment(data = dat.rescaled[dat.rescaled$caseControl == '1',], aes(x = NR, xend = NR, y = 0.02, yend = 0.04), alpha = .2, col = 'dodgerblue')+
  geom_ribbon(aes(NR, ymin = plogis(plo), ymax = plogis(phi)), fill = 'lightgrey')+
  geom_line()+
  ylab('Probability of healthy')+xlab('Mean gut relatedness')+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.grid = element_blank(),
        #legend.position = 'top', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))+
  annotate('text', x = 0.4, y = 0.75, label = paste0('slope = ', b1, ', pval = ', pvalb1, '***'), size = 1.8, hjust = 0)+
  ylim(c(0,0.75))
  
p.predict.NR



# Let's do the same for NRSPO now

newdat.2 <- expand.grid( # Make a new dataframe for predictions
  NRSPO=seq(0.08, 0.5, 0.01),         # Make range of NRSPO within range of observed values
  NR = mean(dat.rescaled$NR), # condition on typical (i.e. mean) value of NR
  caseControl = NA
)

newdat.2$caseControl <- predict(fm1,newdat.2,re.form=NA) # Make predictions

plot(newdat.2$caseControl~ newdat.2$NRSPO, type = 'l', ylim = c(-2, 0.5)) # Ok, what I get corresponds to what allEffect() wrapper gives
plot(allEffects(fm1), type = 'link')


mm.2 <- model.matrix(terms(fm1),newdat.2)
## or newdat$distance <- mm %*% fixef(fm1)
pvar1.2 <- diag(mm.2 %*% tcrossprod(vcov(fm1),mm.2)) # Thus is to produce CI based on fixed-effects uncertainty ONLY
tvar1.2 <- pvar1.2+VarCorr(fm1)$disease_type[1]+VarCorr(fm1)$disease_type[5]+VarCorr(fm1)$disease_type[9]  # This adds the variance of random effects. Adapted from Ben Bolker's code but here for more complex model. Because we let the model estimate variance in intercept, NR slope and NRSPO slope with disease_type, we must add those three sources of variance.
cmult <- 1.96 ## could use 1.96
newdat.2 <- data.frame(
  newdat.2
  , plo = newdat.2$caseControl-cmult*sqrt(pvar1.2)
  , phi = newdat.2$caseControl+cmult*sqrt(pvar1.2)
  , tlo = newdat.2$caseControl-cmult*sqrt(tvar1.2)
  , thi = newdat.2$caseControl+cmult*sqrt(tvar1.2)
)
#plot confidence
g0 <- ggplot(newdat.2, aes(x=NRSPO, y=plogis(caseControl)))+geom_point()
g0 + geom_pointrange(aes(ymin = plogis(plo), ymax = plogis(phi)))+
  labs(title="CI based on fixed-effects uncertainty ONLY")
#plot prediction
g0 + geom_pointrange(aes(ymin = plogis(tlo), ymax = plogis(thi)))+
  labs(title="CI based on FE uncertainty + RE variance")

plot(allEffects(fm1), type = 'response') # ok, what I get manually looks similar to what allEffects() wrapper gives



out.fix.format<- summary(m.glmm.disease)[[10]] %>%
  as.data.frame() %>%
  add_rownames('Effect') %>%
  as.data.frame()


b1 = round(out.fix.format[3,2],2)
pvalb1 = round(out.fix.format[3,5],2)

p.predict.NRSPO<- ggplot(newdat.2, aes(x=NRSPO, y=plogis(caseControl)))+
  #geom_segment(data = dat.rescaled, aes(x = NRSPO, xend = NRSPO, y = 0.1, yend = 0.116), alpha = .2)+
  geom_segment(data = dat.rescaled[dat.rescaled$caseControl == '0',], aes(x = NRSPO, xend = NRSPO, y = 0, yend = 0.02), alpha = .2, col = 'firebrick')+
  geom_segment(data = dat.rescaled[dat.rescaled$caseControl == '1',], aes(x = NRSPO, xend = NRSPO, y = 0.02, yend = 0.04), alpha = .2, col = 'dodgerblue')+
  geom_ribbon(aes(NRSPO, ymin = plogis(plo), ymax = plogis(phi)), fill = 'lightgrey')+
  geom_line()+
  ylab('Probability of healthy')+xlab('Relatedness weighed mean\nsporulation score')+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.grid = element_blank(),
        #legend.position = 'top', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))+
  annotate('text', x = 0.1, y = 0.62, label = paste0('slope = ', b1, ', pval = ', pvalb1), size = 1.8, hjust = 0)

p.predict.NRSPO



# Assemble ----

library(patchwork)
#pdf('./output/figures/Figure2_microbiome.pdf', width = 16/2.55, height = 13/2.55)
#(p.NR | p.NRSPO | (p.all.NR|p.all.NRSPO)/p.predict.NR / p.predict.NRSPO) +
#  plot_annotation(tag_levels = 'A') &
#  theme(plot.tag = element_text(size = 7, face = 'bold'))
#dev.off()


pdf('./output/figures/Figure2_microbiome.pdf', width = 14/2.55, height = 11.5/2.55)
(p.NR | p.NRSPO | (p.predict.NR / p.predict.NRSPO/plot_spacer())) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 7, face = 'bold'))
dev.off()



library(cowplot)
leg<- get_legend(
p.NR+
  theme(legend.position = 'right',
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "cm"))+
        scale_fill_manual(values = c('firebrick', 'dodgerblue'),
                    labels = c('case (sick)', 'control (healthy)'))
)
  
pdf('./output/figures/legend_microbiomePlot.pdf', width = 3/2.55, height = 1.8/2.55)
plot(leg)
dev.off()




# AUC (cvAUC package analysis) ----

library(pROC)
library(caret)
library(metRology)

# USING cvAUC PACKAGE, WITH FUNCTION ci.cvAUC() TO COMPUTE A CONFIDENCE INTERVAL ON THE AUC

# Package has a useful documentation page to understand exactly what they are doing:
# https://www.rdocumentation.org/packages/cvAUC/versions/1.1.0/topics/ci.cvAUC

# What they do is:
# 1) define the folds (as in: select the row ids to be used in each fold). Ideal is to stratify the folds by outcome, to keep balance of 1/0 cases in each fold
# 2) fit the model on all folds but one
# 3) using that model fit, predict values in the left-out fold, on response scale
# do for each fold in turn
# then in the ci.cvAUC() function:
#   predictions is the vector of all predictions, on response scale (probability)
#   label is the true observed value (on 1/0 scale)
#   fold is the list of row ids for each fold


# Let's use their demo code wrappers. Works for any dataset as long as
# input dataset has a first 'Y' variables of 1/0 outcome
# other variables are those to use in the model
# so if want only 'NR', subset dataset to get only the Y and NR columns
# and of course everytime, subset to keep all or specific disease

# in the code wrapper below, utility functons to get the folds and train/test are defined within the code wrapper
# then the cross validated AUC is passed through ci.cvAUC() to get the mean AUC and confidence intervals

cvAUC_demoWrapper <- function(data, V=10){
  
  .cvFolds <- function(Y, V){  #Create CV folds (stratify by outcome)
    Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
    Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
    folds <- vector("list", length=V)
    for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}		
    return(folds)
  }
  .doFit <- function(v, folds, data){  #Train/test glm for each fold
    fit <- glm(Y~., data=data[-folds[[v]],], family=binomial)
    pred <- predict(fit, newdata=data[folds[[v]],], type="response")
    return(pred)
  }
  folds <- .cvFolds(Y=data$Y, V=V)  #Create folds
  predictions <- unlist(sapply(seq(V), .doFit, folds=folds, data=data))  #CV train/predict
  predictions[unlist(folds)] <- predictions  #Re-order pred values
  # Get CV AUC and confidence interval
  out <- ci.cvAUC(predictions=predictions, labels=data$Y, folds=folds, confidence=0.95)
  return(out)
}

# data here must be dat.rescaled, which has the caseControl variable (0 = sick, 1 = healthy), and the disease_type variable which specific which disease the "cases" were
# in principle could modify the wrapper to select specific cohorts instead
# focal.predictor must be a vector of the predictors you want to use for the glm model fit
# focal.disease must be a vector of disease you want to subset the dataset for (it will select the cases and controls only from the studies focusing on those specific disease)
# V is number of folds to use for the k-fold cross validation, if set V = nrow(data), then that's a leave-one-out analysis

# Note: I'm setting set.seed(42) in the code wrapper, just before computing the folds
# that means that the output will be reproducible because that will reselect the same folds for a given dataset
# so when re-running the cvAUC (or running it with a different predictor), the folds are the same, so the AUC can be compared
# as in, the difference in AUC is NOT due to sampling, it's due to the different predictor
run.cvAUC_demoWrapper<- function(data, focal.predictor, focal.disease , V=10){
  
  d.focal<- data[data$disease_type %in% focal.disease, c('caseControl', focal.predictor)] %>%
    rename(Y = caseControl)
  
  set.seed(42) # by setting the seed, the folds selection will be the same. MUST set it everytime running the AUCcv for a different predictor, so that the AUC computed are comparable, i.e. computed over the same data split
  out<- cvAUC_demoWrapper(d.focal, V = V)
  
  out.df<- data.frame(predictor = paste(focal.predictor, collapse = '_'),
                      disease = paste(focal.disease, collapse = '_'),
                      cvAUC = out$cvAUC,
                      se = out$se,
                      ci.low = out$ci[1],
                      ci.up = out$ci[2])
  
  return(out.df)
  
}

# example
run.cvAUC_demoWrapper(dat.rescaled, c('NR', 'NRSPO'), 'MetS', V = 10)


# Ok, now, just have to run it for each of my four predictors (NR/NRSPO, or Shannon, invSimpson, Richness)
# And for each disease as well as all diseases together
# that's 44 runs
unique(dat.rescaled$disease_type) # 10 diseases
(10*4)+4


ls.diseases<- as.list(unique(dat.rescaled$disease_type))
ls.diseases[[11]]<- unique(dat.rescaled$disease_type)

ls.predictors<- as.list(c('NR', 'Shannon', 'invSimpson', 'Richness'))
#ls.predictors[[5]]<- c('NR', 'NRSPO')


ls.cvAUC.out<- list()

counter = 0
for(i in 1:length(ls.diseases)){
  for(j in 1:length(ls.predictors)){
    
    counter = counter + 1
    
    run.out<- run.cvAUC_demoWrapper(dat.rescaled, ls.predictors[[j]], ls.diseases[[i]], V = 10)
    
    print(counter)
    ls.cvAUC.out[[counter]] = run.out
    
  }
}


ls.cvAUC.out<- do.call('rbind', ls.cvAUC.out)
ls.cvAUC.out$disease[ls.cvAUC.out$disease == 'CRC_HYPER_IBD_ACVD_LC_OBE_T2D_IGC_RA_BD_MetS']<- 'All diseases'
#ls.cvAUC.out<- ls.cvAUC.out[ls.cvAUC.out$predictor != 'NR_NRSPO',]

ls.cvAUC.out$disease<- factor(ls.cvAUC.out$disease, 
                              levels = rev(c('All diseases', 'OBE', 'MetS', 'IBD', 'CRC', 'T2D_IGC', 'LC', 'ACVD', 'RA', 'HYPER', 'BD')))


# Add samples sizes under disease names on y-axis ticks
nsizes.diseases<- dat.rescaled %>%
  group_by(disease_type, caseControl) %>%
  summarise(n = n()) %>%
  as.data.frame() %>%
  ungroup() %>%
  spread(key = caseControl, value = n) %>%
  as.data.frame() %>%
  rename(case = `0`,
         control = `1`) %>%
  rename(disease = disease_type)


nsizes.diseases<- rbind(nsizes.diseases, 
                        data.frame(disease = 'All diseases',
                                   case = table(dat.rescaled$caseControl)[[1]],
                                   control = table(dat.rescaled$caseControl)[[2]]))


ls.cvAUC.out2<- left_join(ls.cvAUC.out, nsizes.diseases, by = 'disease')


ls.cvAUC.out2<- ls.cvAUC.out2 %>%
  mutate(diseasePlot2 = paste0(disease, '\n(c:', case, ', h:', control, ')'))

ls.cvAUC.out2$diseasePlot2[ls.cvAUC.out2$diseasePlot2 == 'T2D_IGC\n(c:280, h:215)']<- "T2D & IGT\n(c:280, h:215)"

unique(ls.cvAUC.out2$diseasePlot2)


ls.cvAUC.out2$diseasePlot2<- factor(ls.cvAUC.out2$diseasePlot2, 
                                    levels = rev(c("All diseases\n(c:1750, h:1184)",
                                                   "OBE\n(c:234, h:180)",
                                                   "MetS\n(c:60, h:26)",
                                                   "IBD\n(c:196, h:108)",
                                                   "CRC\n(c:386, h:205)",
                                                   "T2D & IGT\n(c:280, h:215)",
                                                   "LC\n(c:118, h:111)",
                                                   "ACVD\n(c:214, h:169)",          
                                                   "RA\n(c:88, h:89)",
                                                   "HYPER\n(c:155, h:41)",              
                                                   "BD\n(c:19, h:40)")))

# by definition, a 10-fold will split it in 10, meaning 90% of data goes to training and 10% goes to testing set
# So by knowing the sample size of the cohorts, the size of train/test sets is known


# redefining predictors as factor to have them in right orders on plot
ls.cvAUC.out2$predictor[ls.cvAUC.out2$predictor == 'NR']<- 'Mean Relatedness'
ls.cvAUC.out2$predictor<- factor(ls.cvAUC.out2$predictor,
                                 levels = rev(c('Mean Relatedness', 'Shannon', 'invSimpson', 'Richness')))


# Setting pallette for background coloring area
pal<- c(piratepal('basel', trans = 0.4), piratepal('pony', trans = .5)[9], piratepal('info2', trans = .5)[14])
# Redfining that 'pal' palette with colors that do not rely on transparency settings (for correct pdf conversion)
cols<- NULL
for(i in 1:length(piratepal('basel'))){
  cols<- c(cols, colorRampPalette(c("white", piratepal('basel')[i]))(10)[3])
}
cols<- c(cols, colorRampPalette(c("white", piratepal('pony')[9]))(10)[3])
cols<- c(cols, colorRampPalette(c("white", piratepal('info2')[14]))(10)[3])
pal<- cols


brewer.pal(9, 'BrBG')[9] # seablue
brewer.pal(2, 'BrBG')[2] # marron
brewer.pal(10, 'Paired')[10] # purple
brewer.pal(5, 'Paired')[5] # salmon



# plot it
p.auc<- ggplot(ls.cvAUC.out2, aes(x = cvAUC, y = diseasePlot2, col = predictor, xmin = ci.low, xmax = ci.up))+
  geom_point(position = position_dodge(.8))+
  geom_errorbarh(position = position_dodge(.8), height = 0.2)+
  #scale_color_manual(values = rev(c('darkred','firebrick', 'dodgerblue', 'purple', 'goldenrod')))+
  
  scale_color_manual(values = rev(c(brewer.pal(9, 'BrBG')[9], # seablue
                                    brewer.pal(7, 'Accent')[7], # marron
                                    brewer.pal(10, 'Paired')[10], # purple
                                    brewer.pal(5, 'Paired')[5])))+ # salmon
                                    
  
  #annotate("rect", ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = +Inf,fill = pal[12])+
  #annotate("rect", ymin = 1.5, ymax = 2.5, xmin = -Inf, xmax = +Inf,fill = pal[7])+
  #annotate("rect", ymin = 2.5, ymax = 3.5, xmin = -Inf, xmax = +Inf,fill = pal[10])+
  #annotate("rect", ymin = 3.5, ymax = 4.5, xmin = -Inf, xmax = +Inf,fill = pal[9])+
  #annotate("rect", ymin = 4.5, ymax = 5.5, xmin = -Inf, xmax = +Inf,fill = pal[8])+
  #annotate("rect", ymin = 5.5, ymax = 6.5, xmin = -Inf, xmax = +Inf,fill = pal[11])+
  #annotate("rect", ymin = 6.5, ymax = 7.5, xmin = -Inf, xmax = +Inf,fill = pal[5])+
  #annotate("rect", ymin = 7.5, ymax = 8.5, xmin = -Inf, xmax = +Inf,fill = pal[3])+
  #annotate("rect", ymin = 8.5, ymax = 9.5, xmin = -Inf, xmax = +Inf,fill = pal[2])+
  #annotate("rect", ymin = 9.5, ymax = 10.5, xmin = -Inf, xmax = +Inf,fill = pal[1])+
  #annotate("rect", ymin = 10.5, ymax = 11.5, xmin = -Inf, xmax = +Inf,fill = 'white', colour = 'black')+
  
  # for grey/white shaded background
  annotate("rect", ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = +Inf,fill = 'white')+
  annotate("rect", ymin = 1.5, ymax = 2.5, xmin = -Inf, xmax = +Inf,fill = mygrey)+
  annotate("rect", ymin = 2.5, ymax = 3.5, xmin = -Inf, xmax = +Inf,fill = 'white')+
  annotate("rect", ymin = 3.5, ymax = 4.5, xmin = -Inf, xmax = +Inf,fill = mygrey)+
  annotate("rect", ymin = 4.5, ymax = 5.5, xmin = -Inf, xmax = +Inf,fill = 'white')+
  annotate("rect", ymin = 5.5, ymax = 6.5, xmin = -Inf, xmax = +Inf,fill = mygrey)+
  annotate("rect", ymin = 6.5, ymax = 7.5, xmin = -Inf, xmax = +Inf,fill = 'white')+
  annotate("rect", ymin = 7.5, ymax = 8.5, xmin = -Inf, xmax = +Inf,fill = mygrey)+
  annotate("rect", ymin = 8.5, ymax = 9.5, xmin = -Inf, xmax = +Inf,fill = 'white')+
  annotate("rect", ymin = 9.5, ymax = 10.5, xmin = -Inf, xmax = +Inf,fill = mygrey)+
  annotate("rect", ymin = 10.5, ymax = 11.5, xmin = -Inf, xmax = +Inf,fill = 'white', colour = 'black')+
  
  
  
  geom_point(position = position_dodge(.8), size = 1)+
  geom_errorbarh(position = position_dodge(.8), height = 0.3, size = 0.4)+
  ylab('Disease')+xlab('AUC')+
  theme_bw()+
  geom_vline(xintercept = c(0.5, 0.7), linetype = 'dashed')+
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.4))

p.auc


pdf('./output/figures/Figure3_AUC_CrossValidation.pdf', width = 10/2.55, height = (15/2.55))
p.auc
dev.off()


















# THEORY PLOT ----

r = seq(0, 1, 0.01)
v = seq(0, 1, 0.01)
b = -0.5
c = 1
v_effect = 1


df<- expand.grid(r, v, b, c) %>%
  rename(r = Var1, v = Var2, b = Var3, c = Var4) %>%
  mutate(x = r*((b+(v_effect*v))/(2*c)))


# Basic theory plot
p.map.theory<- ggplot(df, aes(x = v, y = r))+
  geom_tile(aes(fill = x, col = x), size = 0.2)+
  scale_fill_gradient2('Predicted effect\non host',low = 'firebrick', mid = 'white', high='dodgerblue', midpoint = 0)+
  scale_colour_gradient2('Predicted effect\non host',low = 'firebrick', mid = 'white', high='dodgerblue', midpoint = 0)+
  ylab('Relatedness')+xlab('Probability of vertical transmission')+
  scale_y_continuous(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0))+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(size = 0.2, color = 'black'),
        panel.grid = element_line(size = 0.1),
        legend.position = 'right', #c(0.9, 0.1),
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(size = 4),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))

p.map.theory




# Plot of observed data

r.df<- read.table('./data/relatedness_sporulation.txt', header = TRUE, sep = '\t')

p.map.obs<- ggplot(r.df, aes(x = (1-sporulation_score), y = mean_relatedness))+
  geom_point(size = 0.3)+
  #geom_tile(aes(fill = x))+
  #scale_fill_gradient2('Predicted effect\non host',low = 'firebrick', mid = 'white', high='springgreen4', midpoint = 0)+
  ylab('Relatedness estimates')+xlab('(1 - sporulation score)')+
  #scale_y_continuous(expand = c(0, 0))+
  #scale_x_continuous(expand = c(0, 0))+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(size = 0.2, color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'right', #c(0.9, 0.1),
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(size = 4),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))


library(patchwork)

# adding different colors to point for better plotting
keep<- c("Ruminococcus_bromii_62047",
         "Clostridium_sp_61482",
         "Eubacterium_eligens_61678",
         "Guyana_massiliensis_60772",
         "Faecalibacterium_prausnitzii_61481",
         "Roseburia_intestinalis_56239",
         "Coprococcus_comes_61587",
         "Dorea_longicatena_61473",
         "Lachnospiraceae_bacterium_56833",
         "Erysipelotrichaceae_bacterium_55770",
         "Escherichia_coli_58110",
         "Akkermansia_muciniphila_55290",
         "Bacteroides_intestinalis_61596",
         "Alistipes_onderdonkii_55464",
         "Bacteroides_finegoldii_57739",
         "Parabacteroides_goldsteinii_56831",
         "Parabacteroides_merdae_56972",
         "Bacteroides_fragilis_54507",
         "Odoribacter_laneus_62216",
         "Acidaminococcus_intestini_54097",
         "Clostridium_leptum_61499")


r.df$cat<- 'a'
r.df[which(r.df$species_id%in% keep),'cat']<- 'b'

r.df$names<- ''
r.df<- r.df %>% mutate(names = ifelse(species_id %in% keep, species_id, ''))
r.df$names<- gsub('_', ' ', substr(r.df$names,1,nchar(r.df$name)-6))


p.map.obs<- ggplot(r.df, aes(x = (1-sporulation_score), y = mean_relatedness, col = cat))+
  geom_point(size = 0.3)+
  #geom_tile(aes(fill = x))+
  #scale_fill_gradient2('Predicted effect\non host',low = 'firebrick', mid = 'white', high='springgreen4', midpoint = 0)+
  #geom_text(aes(label = names), size = 1)+
  #geom_text_repel(aes(label = names), box.padding = unit(0.1, 'line'), segment.colour = 'black', size = 0.9, segment.size = 0.5)+
  ylab('Relatedness estimates')+xlab('(1 - sporulation score)')+
  scale_color_manual(values = c('darkgrey', 'black'))+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.ticks = element_line(size = 0.2, color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none', #c(0.9, 0.1),
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(size = 4),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))



pdf('./output/figures/Figure1_heatmap.pdf', width = 12/2.55, height = 5/2.55)
(p.map.theory | p.map.obs) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 7, face = 'bold'))
dev.off()




# Version with all names on plot
p.heatmap.obs.names<- ggplot(r.df, aes(x = (1-sporulation_score), y = mean_relatedness, label = species_id))+
  #geom_point(size = 0.1)+
  #geom_text_repel(aes(label = species_id), box.padding = unit(0.1, 'line'), fontface = 'bold', segment.colour = 'black', size = 1.5, segment.size = 0.2)+
  geom_text_repel(aes(label = species_id), box.padding = unit(0.1, 'line'), segment.colour = 'black', size = 1.3, segment.size = 0.2)+
  
  #geom_tile(aes(fill = x))+
  #scale_fill_gradient2('Predicted effect\non host',low = 'firebrick', mid = 'white', high='springgreen4', midpoint = 0)+
  ylab('Relatedness')+xlab('(1 - sporulation score)')+
  #scale_y_continuous(expand = c(0, 0))+
  #scale_x_continuous(expand = c(0, 0))+
  #scale_color_manual(values = rep('black', nrow(r.df)))+
  theme_bw()+
  theme(axis.title = element_text(face = 'bold', size = 9),
        axis.text = element_text(size = 9, color = 'black'),
        axis.ticks = element_line(size = 0.6, color = 'black'),
        panel.grid = element_blank(),
        legend.position = 'right', #c(0.9, 0.1),
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_text(size = 4),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.6))




#library(ggplot2)
#library(plotly)
#library(gapminder)
#install.packages('plotly') # if you haven't installed the package
#library(plotly)
#py <- plotly(username="r_user_guide", key="mw5isa4yqp")
#ggplotly(p.heatmap.obs.names)
#py$ggplotly(p.heatmap.obs.names)
#l <- plotly::ggplotly(p.heatmap.obs.names)
#htmlwidgets::saveWidget(l, file = "./output/figures/index.html")


pdf('./output/figures/Supp_heatmap_names.pdf', width = 15/2.55, height = 20/2.55)
p.heatmap.obs.names
dev.off()




# Map species on theory plot & other versions

r.df<- read.table('./data/relatedness_sporulation.txt', header = TRUE, sep = '\t')
r.df$revSporulation<- 1- r.df$sporulation_score
r.df$spoRescaled = 0.01 + ((0.99 - 0.01) / (max(r.df$revSporulation) - min(r.df$revSporulation))) * (r.df$revSporulation - min(r.df$revSporulation))

p.map.theory+
  geom_point(data = r.df, aes(x = revSporulation, y = mean_relatedness))#+
#geom_text_repel(data = r.df, aes(x = spoRescaled, y = mean_relatedness, label = species_id), box.padding = unit(0.45, 'line'), fontface = 'bold', segment.colour = 'black')

ggplot(r.df, aes(x = revSporulation, y = mean_relatedness))+
  geom_point()
#geom_text_repel(data = r.df, aes(x = spoRescaled, y = mean_relatedness, label = species_id), box.padding = unit(0.45, 'line'), fontface = 'bold', segment.colour = 'black')


leg <- get_legend(p.hm.leg) # library(ggpubr)
p.legend<- as_ggplot(leg)
p.legend


p.fan.theory<- ggplot(df, aes(x = r, y = x, col = v))+
  geom_point()+
  scale_color_gradient2(low = 'firebrick', mid = 'white', high='springgreen4', midpoint = 0.5)+
  geom_point(data = r.df, aes(x = mean_relatedness, y = pred_x), col = 'black')+
  geom_text_repel(data = r.df[r.df$species_id %in% keep,], aes(x = mean_relatedness, y = pred_x, label = species_id), box.padding = unit(0.45, 'line'), fontface = 'bold', segment.colour = 'black', col = 'black')

p.fan.theory



# MCMC general MMMC model ----

# Not included in MS. Another possible approach to analyse the data, using Multimembership model.

load('./output/stats_runs/data_preps_23062021.RData')


mf = 100


prior.1<-list(R=list(V=1, fix=1),
              G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100),
                     G2=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m_1<-MCMCglmm(caseControl~NR,
              random=~idv(prob)+idv(probr),
              data=dat.rescaled,
              family="threshold",
              prior=prior.1,
              nitt = 13000*mf,
              thin = 10*mf,burnin=3000*mf)

summary(m_2)



prior.2<-list(R=list(V=1, fix=1),
              G=list(G1=list(V=1, nu=1, alpha.mu=0, alpha.V=100)))

m_2<-MCMCglmm(caseControl~NR,
              random=~idv(prob),
              data=dat.rescaled,
              family="threshold",
              prior=prior.2,
              nitt = 13000*mf,
              thin = 10*mf,burnin=3000*mf)


save.image('./output/stats_runs/MMMC_HM2_run_July2021.RData')



# ARCHIVED CODE ----

# [deprecated] AllData boxplots ----

p.all.NR<- ggplot(dat.rescaled, aes(x = as.factor(caseControl), y = NR, fill = as.factor(caseControl)))+
  geom_boxplot(width = 0.5, outlier.size = 0.1, lwd = 0.3)+
  ylab('Mean gut relatedness')+xlab('')+
  scale_fill_manual(values = c('lightgrey', 'white'))+
  scale_x_discrete(labels = c('Sick', 'Healthy'))+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1),        panel.grid = element_line(size = 0.1),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))


p.all.NRSPO<- ggplot(dat.rescaled, aes(x = as.factor(caseControl), y = NRSPO, fill = as.factor(caseControl)))+
  geom_boxplot(width = 0.5, outlier.size = 0.1, lwd = 0.3)+
  ylab('Relatedness weighed mean\nsporulation score')+xlab('')+
  scale_x_discrete(labels = c('Sick', 'Healthy'))+
  scale_fill_manual(values = c('lightgrey', 'white'))+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1),
        panel.grid = element_line(size = 0.1),
        legend.position = 'none', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))


# [deprecated] ROC/AUC analysis ----
library(pROC)
library(caret)
library(metRology)

# Using caret package to run cross-validation analysis

# Example to see how it works
# Code wrapper to run the cross validation
rocCV<- function(data, nb.k, reps, focal.disease, focal.predictor){
  
  # Defining the cross validation algorithm parameters
  
  fit.control <- trainControl(method = "repeatedcv", # repeated k-fold validation
                              number = nb.k, repeats = reps, # use: 10 folds, repeat: 10 times
                              summaryFunction = twoClassSummary, classProbs = TRUE) # Use AUC summary stat
  
  
  #fit.control <- trainControl(method = "LOOCV", # for leave-one-out analysis
  #                            summaryFunction = twoClassSummary, classProbs = TRUE)
  
  d2<- data %>%
    filter(disease_type %in% focal.disease) %>%
    select(caseControl, NR, Shannon, invSimpson, Richness) %>%
    rename(focal.predictor = focal.predictor)
  
  d2$caseControl<- ifelse(d2$caseControl == 0, 'sick', 'healthy')
  
  set.seed(42)  
  fit<- train(caseControl ~ focal.predictor, data = d2, method = "glm", family = "binomial", trControl = fit.control)
  #fit$results # ROC corresponds to mean(fit$resample[,'ROC]), and ROCSD = sd(fit$resample[,'ROC])
  
  return(fit)
  
}

# run it
out<- rocCV(dat.rescaled, 10, 10, 'MetS', 'NR')

out$results # ROC = 0.53, ROCSD = 0.2328


mean(out$resample[,'ROC']) # this is an emprical sampling distribution of the AUC
sd(out$resample[,'ROC']) # and this is the standard deviation of the sampling distribution of our mean estimate of the AUC, i.e.the standard error
# then taking the 0.025 and 97.5 percentiles of this empirical distribution gives as the confidence interval

# So as far as I understand, the k-fold cross validation is a way to obtain a mean estimate of the AUC
# and the vector of AUC obtained from each fold/repeat sample of the data is essentially an empirical sampling distribution
# So the standard deviation of that is the standard error of the mean AUC
# and the percentiles of that is the confidence interval.

# Here the CI would be ~ 1.96*0.2328 = 0.456, so the mean AUC would be [0.074-0.986]
quantile(out$resample[,'ROC'], c(0.025, 0.975))
mean(out$resample[,'ROC'])-(1.96*sd(out$resample[,'ROC']))
mean(out$resample[,'ROC'])+(1.96*sd(out$resample[,'ROC']))
# Ok so don't find exactly the same, probably because the exact normal approximation here is not the best
# Like for certain case the t-distribution (= a slightly flattened) normal is better
#
qt.scaled(c(0.025, 0.975), 
          mean = mean(out$resample[,'ROC']), 
          sd = sd(out$resample[,'ROC']), 
          df = length(out$resample[,'ROC']))
# Yep, what I get there is similar to what I get with mean +/- 1.96 SD, but a tad different due to qt.scaled being a t distribution
# The empirical quantiles are a bit shifted towards lower estimates



# Quick sensitivity analysis test for k and reps values, see with MetS which is a small dataset
df.run.rocCV<- data.frame(nb.k = rep(seq(3,10,1), each = 4),
                          reps = rep(c(10, 50, 100, 1000), 8),
                          focal.disease = 'MetS', 
                          ROC = NA, 
                          ROCSD = NA)

for(i in 1:nrow(df.run.rocCV)){
  
  print(i)
  
  t<- rocCV(dat.rescaled, nb.k = df.run.rocCV$nb.k[i], reps = df.run.rocCV$reps[i], focal.disease = df.run.rocCV$focal.disease[i])$results[,c('ROC', 'ROCSD')]
  df.run.rocCV$ROC[i] = t[,'ROC']
  df.run.rocCV$ROCSD[i] = t[,'ROCSD']
  
}


# Not much effect, variance goes larger as k increases, which is expected
ggplot(df.run.rocCV, aes(x = nb.k, y = ROC, ymin = ROC-ROCSD, ymax = ROC+ROCSD, col = as.factor(reps)))+
  geom_point(position = position_dodge(.8))+
  geom_errorbar(position = position_dodge(.8), width = 0.2)+
  scale_color_manual(values = c('red', 'dodgerblue', 'green', 'goldenrod'))



# actual run
# keep k = 10, reps = 1000
nb.k = 10
reps = 1000


df.run.rocCV.fin<- data.frame(nb.k = nb.k,
                              reps = reps,
                              focal.disease = unique(dat.rescaled$disease_type), 
                              ROC = NA, 
                              ROCSD = NA)

# Make dataframe smaller (faster run?)
d.clean = dat.rescaled %>%
  select(caseControl, NR, Shannon, invSimpson, Richness, disease_type) %>%
  mutate(caseControl = ifelse(caseControl == 0, 'sick', 'healthy'))

head(d.clean)

fit.control <- trainControl(method = "repeatedcv", number = nb.k, repeats = reps,
                            summaryFunction = twoClassSummary, classProbs = TRUE)

# re-define code wrappers to work with actual run
rocCV<- function(data, focal.predictor, focal.disease, trainControl){
  
  data<- rename(data, focal.predictor = focal.predictor)
  
  set.seed(42) # re-setting the seed everytime ensures that the dataset is split the same way when running the ROCcv for the different predictors  
  fit<- train(caseControl ~ focal.predictor, data = data[data$disease_type %in% focal.disease,], method = "glm", family = "binomial", trControl = trainControl)
  
  
  d.out<- data.frame(disease = focal.disease,
                     predictor = focal.predictor,
                     ROC = fit$results[,'ROC'],
                     ROCSD = fit$results[,'ROCSD'],
                     mean.train.n = round(mean(lengths(fit$control$index)),0), # mean training set size
                     mean.test.n = round(mean(lengths(fit$control$indexOut)),0)) # mean test set size
  
  return(list(d.out, fit))
}


rocCV.severalDiseases<- function(data, focal.predictor, focal.disease, disease.tag = 'several', trainControl){
  
  data<- rename(data, focal.predictor = focal.predictor)
  
  set.seed(42) # re-setting the seed everytime ensures that the dataset is split the same way when running the ROCcv for the different predictors  
  fit<- train(caseControl ~ focal.predictor, data = data[data$disease_type %in% focal.disease,], method = "glm", family = "binomial", trControl = trainControl)
  
  
  d.out<- data.frame(disease = disease.tag,
                     predictor = focal.predictor,
                     ROC = fit$results[,'ROC'],
                     ROCSD = fit$results[,'ROCSD'],
                     mean.train.n = round(mean(lengths(fit$control$index)),0), # mean training set size
                     mean.test.n = round(mean(lengths(fit$control$indexOut)),0)) # mean test set size
  
  return(list(d.out, fit))
  
}


preds = c('NR', 'Shannon', 'invSimpson', 'Richness') # four predictors to run the cross-validation over
dis = unique(d.clean$disease_type) # 19 diseases to run the cross validation over

# to store results
dfin<- data.frame(disease = character(),
                  predictor = character(),
                  ROC = numeric(),
                  ROCSD = numeric(),
                  mean.train.n = numeric(),
                  mean.test.n = numeric())

dfin.list<- vector('list')

# In the loop, redefine the set.seed() everytime ensures the same fold split is done for each test over each predictor
# can be verified when looking within $result output, can see the seeds of the out/in rows are the same

counter = 0
for(p in 1:length(preds)){
  
  foc.predictor = preds[p]
  
  for(d in 1:length(dis)){
    
    counter = counter + 1
    print(counter)
    
    foc.disease = dis[d]
    
    roc.cv<- rocCV(d.clean, focal.predictor = foc.predictor, focal.disease = foc.disease, trainControl = fit.control)
    
    dfin<- rbind(dfin, roc.cv[[1]])
    
    dfin.list[[counter]]<- roc.cv[[2]]$resample
    
  }
  
}


# Finally, add the cross validation when pooling all diseases
roc.cv.all_NR<- rocCV.severalDiseases(d.clean, focal.predictor = 'NR', focal.disease = dis, disease.tag = 'All diseases', trainControl = fit.control)
roc.cv.all_Shannon<- rocCV.severalDiseases(d.clean, focal.predictor = 'Shannon', focal.disease = dis, disease.tag = 'All diseases', trainControl = fit.control)
roc.cv.all_invSimpson<- rocCV.severalDiseases(d.clean, focal.predictor = 'invSimpson', focal.disease = dis, disease.tag = 'All diseases', trainControl = fit.control)
roc.cv.all_Richness<- rocCV.severalDiseases(d.clean, focal.predictor = 'Richness', focal.disease = dis, disease.tag = 'All diseases', trainControl = fit.control)



dfin2<- rbind(dfin,
              roc.cv.all_NR[[1]],
              roc.cv.all_Shannon[[1]],
              roc.cv.all_invSimpson[[1]],
              roc.cv.all_Richness[[1]])


dfin.list2<- dfin.list
dfin.list2[[41]]<-  roc.cv.all_NR[[2]]$resample
dfin.list2[[42]]<-  roc.cv.all_Shannon[[2]]$resample
dfin.list2[[43]]<-  roc.cv.all_invSimpson[[2]]$resample
dfin.list2[[44]]<-  roc.cv.all_Richness[[2]]$resample


hist(dfin.list2[[11]][,'ROC'])


rm(dfin.list)
rm(dat.rescaled) # need only d.clean for the AUC analysis
rm(prob)
rm(prob.rescaled)
rm(roc.cv.all_invSimpson)
rm(roc.cv.all_NR)
rm(roc.cv.all_Richness)
rm(roc.cv.all_Shannon)
save.image('./output/stats_runs/AUC_crossValidation_runClean.RData')
load('./output/stats_runs/AUC_crossValidation_runClean.RData')

# looking at an example of resampling distribution
quantile(dfin.list2[[41]][,'ROC'], c(0.025, 0.975))
hist(dfin.list2[[41]][,'ROC'])

qt.scaled(c(0.025), mean = roc.cv.all_NR[[1]]$ROC, sd = roc.cv.all_NR[[1]]$ROCSD, df = 999)
qt.scaled(c(0.975), mean = roc.cv.all_NR[[1]]$ROC, sd = roc.cv.all_NR[[1]]$ROCSD, df = 999)

# Ok, those are super close, and the resampling distribution looks very well normal on this example
# so think should be fine, use the "empricial 95%CI", that is, the percentile distribution of the resampling distribution
# this way what this is is transparent

# Plot it
dfin2$disease<- factor(dfin2$disease, 
                       levels = rev(c('All diseases', 'OBE', 'MetS', 'IBD', 'CRC', 'T2D_IGC', 'LC', 'ACVD', 'RA', 'HYPER', 'BD')))

dfin2$predictor[dfin2$predictor == 'NR']<- 'Mean Relatedness'


dfin2$predictor<- factor(dfin2$predictor,
                         levels = rev(c('Mean Relatedness', 'Shannon', 'invSimpson', 'Richness')))

dfin3<- dfin2
head(dfin3)
dfin3$cilower.t<- qt.scaled(c(0.025), mean = dfin3$ROC, sd = dfin3$ROCSD, df = 999)
dfin3$ciupper.t<- qt.scaled(c(0.975), mean = dfin3$ROC, sd = dfin3$ROCSD, df = 999)


getci.emp<- function(x){quantile(x[,'ROC'], c(0.025, 0.975))}
cis.emp<- as.data.frame(do.call('rbind', lapply(dfin.list2, getci.emp))) %>%
  rename(cilower.emp = `2.5%`,
         ciupper.emp = `97.5%`)

dfin3<- cbind(dfin3, cis.emp)


# Using percentile intervals of the resampling distribution rather than the qt.scaled quantiles
# because not 100% sure about the statistical property of the repeated k-fold resampling
# the percentile interval will be 100% clear and transparent about the method
# The two correlate so not worried about it anyway

plot(cilower.emp ~ cilower.t, data = dfin3)


# Add samples sizes under disease names on y-axis ticks
nsizes.diseases<- dat.rescaled %>%
  group_by(disease_type, caseControl) %>%
  summarise(n = n()) %>%
  as.data.frame() %>%
  ungroup() %>%
  spread(key = caseControl, value = n) %>%
  as.data.frame() %>%
  rename(case = `0`,
         control = `1`) %>%
  rename(disease = disease_type)


nsizes.diseases<- rbind(nsizes.diseases, 
                        data.frame(disease = 'All diseases',
                                   case = table(dat.rescaled$caseControl)[[1]],
                                   control = table(dat.rescaled$caseControl)[[2]]))


dfin3<- left_join(dfin3, nsizes.diseases, by = 'disease')

dfin3<- dfin3 %>%
  mutate(diseasePlot2 = paste0(disease, '\n(sick:', case, ', healthy:', control, ')'))


dfin3$diseasePlot2[dfin3$diseasePlot2 == 'T2D_IGC\n(sick:280, healthy:215)']<- "T2D & IGT\n(sick:280, healthy:215)"


dfin3$diseasePlot2<- factor(dfin3$diseasePlot2, 
                            levels = rev(c("All diseases\n(sick:1750, healthy:1184)",
                                           "OBE\n(sick:234, healthy:180)",
                                           "MetS\n(sick:60, healthy:26)",
                                           "IBD\n(sick:196, healthy:108)",
                                           "CRC\n(sick:386, healthy:205)",
                                           "T2D & IGT\n(sick:280, healthy:215)",
                                           "LC\n(sick:118, healthy:111)",
                                           "ACVD\n(sick:214, healthy:169)",          
                                           "RA\n(sick:88, healthy:89)",
                                           "HYPER\n(sick:155, healthy:41)",              
                                           "BD\n(sick:19, healthy:40)")))

# by definition, a 10-fold will split it in 10, meaning 90% of data goes to training and 10% goes to testing set
# So by knowing the sample size of the cohorts, the size of train/test sets is known

library(yarrr)
pal<- c(piratepal('basel', trans = 0.4), piratepal('pony', trans = .5)[9], piratepal('info2', trans = .5)[14])

# Plot
p.auc<- ggplot(dfin3, aes(x = ROC, y = diseasePlot2, col = predictor, xmin = cilower.emp, xmax = ciupper.emp))+
  geom_point(position = position_dodge(.8))+
  geom_errorbarh(position = position_dodge(.8), height = 0.2)+
  scale_color_manual(values = rev(c('firebrick', 'dodgerblue', 'purple', 'goldenrod')))+
  annotate("rect", ymin = 0.5, ymax = 1.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[12])+
  annotate("rect", ymin = 1.5, ymax = 2.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[7])+
  annotate("rect", ymin = 2.5, ymax = 3.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[10])+
  annotate("rect", ymin = 3.5, ymax = 4.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[9])+
  annotate("rect", ymin = 4.5, ymax = 5.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[8])+
  annotate("rect", ymin = 5.5, ymax = 6.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[11])+
  annotate("rect", ymin = 6.5, ymax = 7.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[5])+
  annotate("rect", ymin = 7.5, ymax = 8.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[3])+
  annotate("rect", ymin = 8.5, ymax = 9.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[2])+
  annotate("rect", ymin = 9.5, ymax = 10.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = pal[1])+
  annotate("rect", ymin = 10.5, ymax = 11.5, xmin = -Inf, xmax = +Inf, alpha = .2,fill = 'white', colour = 'black')+
  geom_point(position = position_dodge(.8))+
  geom_errorbarh(position = position_dodge(.8), height = 0.3)+
  ylab('Disease')+
  theme_bw()+
  geom_vline(xintercept = c(0.5, 0.7), linetype = 'dashed')+
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 6),
        panel.grid = element_line(size = 0.1),
        legend.position = 'right',
        legend.justification = 'top',
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.35))
p.auc

pdf('./output/figures/Figure3_AUC_CrossValidation.pdf', width = 10/2.55, height = (15/2.55))
p.auc
dev.off()

