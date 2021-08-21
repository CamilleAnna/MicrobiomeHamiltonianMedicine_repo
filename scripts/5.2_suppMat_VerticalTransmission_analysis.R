# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
#             Using evolutionary theory to predict microbes’ effect on host health
#                           Camille Simonet & Luke McNally
#         
#                            SCRIPT: supplementary materia analysis:
#       Correlation of vertical transmission and sporulation scores, check with phylogenetic mixed model
#                                 last edit: 21/08/2021
#
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


local_project_dir="/Users/s1687811/Documents/PhD/Research/MicrobiomeHamiltonianMedicine/MicrobiomeHamiltonianMedicine_repo"
setwd(local_project_dir)
source('./scripts/sourced_packages.R')
library(ape)
library(MCMCglmm)
library(readxl)



midas.info<- read.table('./data/species_info.txt', header=TRUE)
head(midas.info)
midas.info[grep('55206',midas.info$species_id), ]


spos<- read.table('./data/SPORULATION_SCORES_all_midas_species.txt', header=TRUE, sep = '\t') %>%
  rename(species_id = species)


# As in Browne's 2016 PAPER: FOCUS ON SPECIES DETECTED IN 3 OR MORE MOTHER-INFANT PAIRS

vt<- read_excel('./data/Midas_transmission_data.xls') %>%
  filter(tax_level == 'species',
         time_point == '12M') %>%
  mutate(species_id = gsub(')', '', fixed = TRUE, gsub('._', '_', fixed = TRUE, gsub('(id:', '', gsub(' ', '_', tax_group), fixed = TRUE)))) %>%
  filter(n1 > 3) %>%
  select(species_id, prop1) %>%
  rename(vt = prop1)


vt$species_id[grep('55206', vt$species_id)]<- midas.info[grep('55206',midas.info$species_id), 1]
vt<- left_join(vt, spos, by ='species_id')
vt$genus<- do.call('rbind', strsplit(vt$species_id, '_'))[,1]


tree<- read.tree('./data/midas_tree_renamed.newick')
tree<- drop.tip(tree, tree$tip.label[!tree$tip.label %in% vt$species_id])
phylogeny<- chronopl(tree, lambda = 0)
phylogeny<- makeNodeLabel(phylogeny)
Ainv<-inverseA(phylogeny, scale=FALSE)$Ainv


prior.1<- list(R=list(R1 = list(V=diag(1), nu = 0.002)),  # a single residual variance
               G=list(G1 = list(V=diag(1), nu = 0.002)))  # a single random effect

div = 1
nitt = 10500000/div
burnin = 500000/div
thin = ceiling(50000/div)
nitt
(nitt-burnin)/thin

vt<- as.data.frame(vt)


# Log transform response and run as gaussian

vt$logvtp1<- log(vt$vt+1)

hist(vt$logvtp1)
hist(vt$vt)


m.logvtp1<- MCMCglmm(logvtp1  ~ 1 + sporulation_score,
                random = ~species_id,
                ginverse = list(species_id=Ainv),
                data = vt,
                prior = prior.1,
                family=c("gaussian"),
                start=list(QUASI=FALSE),
                #pl = TRUE, pr = TRUE, nodes = 'ALL',
                DIC = TRUE,
                nitt=nitt, thin=thin, burnin=burnin)


plot(m.logvtp1)
save.image('./output/stats_runs/VerticalTransmission_phylogenetic_model.RData')



# Output summary and figure for supplementary material ----

get.effects.M<- function(model.output, model = ''){
  extract.d<- rbind(
    rbind(summary(model.output)$Gcovariances, summary(model.output)$Rcovariances) %>% as.data.frame() %>% mutate(effect = rownames(.)) %>% mutate(pMCMC = NA),
    summary(model.output)$solutions %>% as.data.frame() %>% mutate(effect = rownames(.)))
  extract.d<- extract.d %>%
    select(effect, post.mean, `l-95% CI`, `u-95% CI`, eff.samp, pMCMC) %>%
    mutate(model_name = model)
  
  extract.d[2:6]<- round(extract.d[2:6], 5)
  
  return(extract.d)
}

tab.VT.phylo<- get.effects.M(m.logvtp1) %>%
  select(-model_name) %>%
  mutate_if(is.numeric, funs(formatC(., digit = 3, format = 'f'))) %>%
  rename(Effect = effect,
         `Posterior mean` = post.mean,
         `CI95% lower` = `l-95% CI`,
         `CI95% upper` = `l-95% CI`,
         `Effective sampling` = eff.samp)

row.names(tab.VT.phylo)<- NULL
tab.VT.phylo$pMCMC<- gsub('NA', '', tab.VT.phylo$pMCMC)

tab.supp.text<- kable(tab.VT.phylo, "latex", booktabs = T,
                      caption = 'Model summary') %>%
  footnote(c("CI95%: 95% credible interval of the posterior distribution",
             "pMCMC: taken as twice the posterior probability that the estimate is negative"),
           fixed_small_size = TRUE, general_title = "") %>%
  kable_styling(latex_options = c("HOLD_position"))
  

fileConn<-file("./output/figures/TABLE_SuppText.tex")
writeLines(tab.supp.text, fileConn)
close(fileConn)


mycolors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", 'green', 'black')


pdf('./output/figures/FigSuppText.pdf', width = 7/2.55, height = 5/2.55)
ggplot(vt, aes(x = sporulation_score, y =  vt, col = genus))+
  ylim(0,1)+
  geom_point(size = 1)+
  #scale_color_brewer(palette = c('Dark2', 'firebrick', 'dodgerblue'))+
  scale_color_manual(values = mycolors)+
  theme_bw()+
  theme(axis.title = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.grid = element_line(size = 0.1),
        #legend.position = 'top', #c(0.9, 0.1),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 3),
        legend.title = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.3))

dev.off()

