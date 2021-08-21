# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#             Using evolutionary theory to predict microbes’ effect on host health
#                           Camille Simonet & Luke McNally
#         
#                            SCRIPT: Scripts sourced in 4.2_statistical_analyses
#                 Makes supplementary tables andoutput a Latex compatible format of them in ./output/figures
#                               last edit: 21/08/2021
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #




library(kableExtra)

funFormat<- function(x, digits){ifelse(abs(x)<0.01,
                                       formatC(x, digit = digits, format = 'e'),
                                       formatC(x, digit = digits, format = 'f'))}




# TABLE S2 ----
# Random binomial regression summary, with disease type fitted as random effects on intercept and slopes

mod<- m.glmm.disease

# formatting random effects
out.rand.format<- summary(mod)[[13]] %>%
  as.data.frame() %>%
  filter(is.na(var2)) %>%
  select(c('grp','var1', 'vcov', 'sdcor')) %>%
  rename(Group = grp,
         Effect = var1,
         Variance = vcov,
         Std.Dev = sdcor) %>%
 mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
 #mutate_if(is.numeric, funs(funFormat(., 2)))



out.rand.format$Group<- c(paste0(out.rand.format$Group[1], ' (', summary(mod)[[9]], ')'), rep('', nrow(out.rand.format)-1))
out.rand.format$Group<- gsub('_', ' ', x = out.rand.format$Group)

# fixed effects
out.fix.format<- summary(mod)[[10]] %>%
  as.data.frame() %>%
  add_rownames('Effect') %>%
  as.data.frame() %>%
  #mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'e')))
  mutate_if(is.numeric, funs(funFormat(., 2)))


# AIC & Co

modChecks<- summary(mod)[[14]] %>%
  as.data.frame() %>%
  add_rownames('Measure') %>%
  spread(Measure, 2) %>%
  mutate_if(is.numeric, funs(formatC(., digit = 2, format = 'f'))) %>%
  as.data.frame() 



# Format for kable
# formatting for it to work with kable: putting colnames as a row and adding foo columns to align with fixed effects

out.rand.format2<- rbind(colnames(out.rand.format), out.rand.format)
out.rand.format2$foo1 = ''
out.rand.format2$foo2 = ''
out.rand.format2$foo3 = ''
out.rand.format2<- out.rand.format2 %>%
  select(foo1, Group, Effect, Variance, Std.Dev, foo2, foo3)


colnames(out.fix.format)[5]<- 'p-value'
colnames(out.fix.format)[4]<- 'z-value'
out.fix.format2<- rbind(colnames(out.fix.format), out.fix.format)
out.fix.format2$foo = ''
out.fix.format2$foo0 = ''
out.fix.format2<- out.fix.format2 %>%
  select(foo, foo0, Effect, Estimate, `Std. Error`, `z-value`, `p-value`)


modChecks<- rbind(colnames(modChecks), modChecks)
modChecks$f1 = ''
modChecks$f2 = ''
modChecks<- modChecks %>% select(f1, f2, AIC, BIC, deviance, df.resid, logLik)


colnames(out.fix.format2)<- LETTERS[1:7]
colnames(out.rand.format2)<- LETTERS[1:7]
colnames(modChecks)<- LETTERS[1:7]

out.fin<- rbind(out.fix.format2, out.rand.format2, modChecks)



tab.s2<- kable(out.fin, "latex", booktabs = T, escape = FALSE,
               col.names = NULL,
               #longtable = T,
               caption = "Random binomial regression summary, with disease type fitted as random effects on intercept and slopes (model syntax: caseControl $\\sim$ 1 + NR + NRSPO + (1 + NR + NRSPO $\\mid$ disease type), family = binomial (link = logit)).") %>%
  kable_styling() %>%
  pack_rows("Fixed effects", 1, 4) %>%
  pack_rows("Random effects", 5, 8) %>%
  pack_rows("Model quality", 9, 9)


fileConn<-file("./output/figures/TABLE_HM2_S2.tex")
writeLines(tab.s2, fileConn)
close(fileConn)






# TABLE S3 ----
# LRT test to test for random slope in model with diseases as random effects


mod.full<- m.glmm.disease
mod.reduced<- m.glmm.disease.fixedSlopes

anova.tab<- as.data.frame(anova(mod.reduced, mod.full)) %>%
  #mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f'))) %>%
  mutate_if(is.numeric, funs(funFormat(., 2))) %>%
  mutate(Model = c('Reduced model (fixed slopes)',
                   'Full model (random slopes)')) %>%
  rename(`N.o. parameters` = npar,
         `p-value` = `Pr(>Chisq)`) %>%
  select(Model, `N.o. parameters`, AIC, BIC, logLik, deviance, Chisq, Df, `p-value`)

anova.tab$Chisq<- gsub('NA', '', anova.tab$Chisq)
anova.tab$Df<- gsub('NA', '', anova.tab$Df)
anova.tab$`p-value`<- gsub('NA', '', anova.tab$`p-value`)

row.names(anova.tab)<- NULL


tab.s3<- kable(anova.tab, "latex", booktabs = T, escape = FALSE,
               #col.names = NULL,
               #longtable = T,
               caption = "Likelihood Ratio Test comparing model with and without (disease) random slopes (reduced model syntax: caseControl $\\sim$ 1 + NR + NRSPO + (1 $\\mid$ disease type), family = binomial (link = logit)).") %>%
  kable_styling() %>%
  pack_rows("LRT fixed vs. random slopes (disease type)", 1, 2)


fileConn<-file("./output/figures/TABLE_HM2_S3.tex")
writeLines(tab.s3, fileConn)
close(fileConn)








# TABLE S4 ----

# Random binomial regression summary, with COHORT type fitted as random effects on intercept and slopes


mod<- m.glmm.cohort

# formatting random effects
out.rand.format<- summary(mod)[[13]] %>%
  as.data.frame() %>%
  filter(is.na(var2)) %>%
  select(c('grp','var1', 'vcov', 'sdcor')) %>%
  rename(Group = grp,
         Effect = var1,
         Variance = vcov,
         Std.Dev = sdcor) %>%
  #mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))
  mutate_if(is.numeric, funs(funFormat(., 2)))
  


out.rand.format$Group<- c(paste0(out.rand.format$Group[1], ' (', summary(mod)[[9]], ')'), rep('', nrow(out.rand.format)-1))
out.rand.format$Group<- gsub('_', ' ', x = out.rand.format$Group)



# fixed effects
out.fix.format<- summary(mod)[[10]] %>%
  as.data.frame() %>%
  add_rownames('Effect') %>%
  as.data.frame() %>%
  mutate_if(is.numeric, funs(funFormat(., 2)))
  #mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f')))


# AIC & Co

modChecks<- summary(mod)[[14]] %>%
  as.data.frame() %>%
  add_rownames('Measure') %>%
  spread(Measure, 2) %>%
  #mutate_if(is.numeric, funs(formatC(., digit = 2, format = 'f'))) %>%
  mutate_if(is.numeric, funs(funFormat(., 2))) %>%
  as.data.frame() 



# Format for kable
# formatting for it to work with kable: putting colnames as a row and adding foo columns to align with fixed effects

out.rand.format2<- rbind(colnames(out.rand.format), out.rand.format)
out.rand.format2$foo1 = ''
out.rand.format2$foo2 = ''
out.rand.format2$foo3 = ''
out.rand.format2<- out.rand.format2 %>%
  select(foo1, Group, Effect, Variance, Std.Dev, foo2, foo3)


colnames(out.fix.format)[5]<- 'p-value'
colnames(out.fix.format)[4]<- 'z-value'
out.fix.format2<- rbind(colnames(out.fix.format), out.fix.format)
out.fix.format2$foo = ''
out.fix.format2$foo0 = ''
out.fix.format2<- out.fix.format2 %>%
  select(foo, foo0, Effect, Estimate, `Std. Error`, `z-value`, `p-value`)


modChecks<- rbind(colnames(modChecks), modChecks)
modChecks$f1 = ''
modChecks$f2 = ''
modChecks<- modChecks %>% select(f1, f2, AIC, BIC, deviance, df.resid, logLik)


colnames(out.fix.format2)<- LETTERS[1:7]
colnames(out.rand.format2)<- LETTERS[1:7]
colnames(modChecks)<- LETTERS[1:7]

out.fin<- rbind(out.fix.format2, out.rand.format2, modChecks)



tab.s4<- kable(out.fin, "latex", booktabs = T, escape = FALSE,
               col.names = NULL,
               #longtable = T,
               caption = "Random binomial regression summary, with cohort  fitted as random effects on intercept and slopes (model syntax: caseControl $\\sim$ 1 + NR + NRSPO + (1 + NR + NRSPO $\\mid$ cohort), family = binomial (link = logit)).") %>%
  kable_styling() %>%
  pack_rows("Fixed effects", 1, 4) %>%
  pack_rows("Random effects", 5, 8) %>%
  pack_rows("Model quality", 9, 9)


fileConn<-file("./output/figures/TABLE_HM2_S4.tex")
writeLines(tab.s4, fileConn)
close(fileConn)


# TABLE S5 ----

# LRT test to test for random slope in model with cohorts as random effects

mod.full<- m.glmm.cohort
mod.reduced<- m.glmm.cohort.fixedSlopes

anova.tab<- as.data.frame(anova(mod.reduced, mod.full)) %>%
  #mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f'))) %>%
  mutate_if(is.numeric, funs(funFormat(., 2))) %>%
  mutate(Model = c('Reduced model (fixed slopes)',
                   'Full model (random slopes)')) %>%
  rename(`N.o. parameters` = npar,
         `p-value` = `Pr(>Chisq)`) %>%
  select(Model, `N.o. parameters`, AIC, BIC, logLik, deviance, Chisq, Df, `p-value`)

anova.tab$Chisq<- gsub('NA', '', anova.tab$Chisq)
anova.tab$Df<- gsub('NA', '', anova.tab$Df)
anova.tab$`p-value`<- gsub('NA', '', anova.tab$`p-value`)

row.names(anova.tab)<- NULL


tab.s5<- kable(anova.tab, "latex", booktabs = T, escape = FALSE,
               #col.names = NULL,
               #longtable = T,
               caption = "Likelihood Ratio Test comparing model with and without (cohorts) random slopes (reduced model syntax: caseControl $\\sim$ 1 + NR + NRSPO + (1 $\\mid$ cohort), family = binomial (link = logit)).") %>%
  kable_styling() %>%
  pack_rows("LRT fixed vs. random slopes (cohort)", 1, 2)


fileConn<-file("./output/figures/TABLE_HM2_S5.tex")
writeLines(tab.s5, fileConn)
close(fileConn)









