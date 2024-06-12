##  ###################################################  ##
##  Beta-diversity measures - btwn site and depths       ##
##                                                       ##
##  Author: Ming Sheng - 15 Nov, 2022                    ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  broom v 0.8.0                                        ##
##  patchwork v 1.1.1                                    ##
##  microbiome v 1.18.0                                  ##
##  purrr v 0.3.4                                        ##
##  corncob v 0.2.0                                      ##
##  indicspecies v 1.7.12                                ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(purrr); packageVersion("purrr")
library(corncob); packageVersion("corncob")
library(indicspecies); packageVersion("indicspecies")
library(pairwiseAdonis); packageVersion("pairwiseAdonis")
library(ggfortify)
library(metagMisc)


# read in phyloseq, filter, format, and transform to robust clr ####
ps <- readRDS('./Output/noncontam_ps_object.RDS')

# filter to archaea and to set threshold
# threshold: prevalence of > 1% AND reads > 5
filt_list <- (microbiome::prevalence(ps) > 0.01) & (microbiome::abundances(ps) %>% rowSums() > 5)

ps_filtered <- readRDS('./Output/noncontam_ps_object.RDS') %>% 
  subset_samples(Depth != 'NA') %>% # depth == NA are mock/blanks
  subset_taxa(filt_list) %>%
  subset_taxa(Class != 'NA')

# subset to archaea group
ps_arch <- ps_filtered %>% subset_taxa(Kingdom == 'Archaea') %>% subset_taxa(Class != 'NA') %>% subset_samples(Site != 'NA')

# transform to robust clr
arch_ra <- ps_arch@otu_table %>% decostand(method = 'rclr')

# metadata file
full_meta <- ps_arch@sam_data %>% data.frame()
mdf_meta <- full_meta[4:11] %>% mutate_all(., function(x) {as.numeric(x)})

# set depth groupings
full_meta$Depth <- as.numeric(full_meta$Depth)
full_meta <- full_meta[order(full_meta$Depth),]
full_meta %>% count(Depth)
full_meta$grouping1 <- c(rep('10-30', 88), rep('40-60', 66), rep('70-100', 39))

#### PCA ####
otu_res <- prcomp(arch_ra, center = T)
en <- envfit(otu_res, mdf_meta, permutations = 999, na.rm = T)

summary(otu_res)
plot(otu_res)

ordiplot(otu_res)
plot(en)

variance = otu_res$sdev^2 / sum(otu_res$sdev^2)
qplot(c(1:5), variance[1:5]) + 
  geom_col() +
  xlab("Principal Component") +
  ylab("Variance Explained") + 
  ggtitle("Scree Plot") +
  ylim(0,0.5) + 
  theme_bw()

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
data.scores$Depth = mdf_meta$Depth
data.scores$Site = mdf_meta$Site

data.scores = data.scores[order(rownames(data.scores)),]
full_meta = full_meta[order(rownames(full_meta)),]
rownames(full_meta) == rownames(data.scores)

data.scores$grouping1 = full_meta$grouping1 # visualize depth groupings

PCAloadings <- data.frame(Variables = rownames(otu_res$rotation), otu_res$rotation)
PCAloadings <- PCAloadings %>% dplyr::select(PC1, PC2)
PCAloadings$strength <- sqrt((PCAloadings$PC1)^2 + (PCAloadings$PC2)^2)
PCAloadings <- PCAloadings[order(PCAloadings$strength, decreasing = TRUE),]

top5_loadings <- PCAloadings[1:5,]

gg = ggplot(data = data.scores, aes(x = PC1, y = PC2)) + 
  geom_point(data = data.scores, aes(colour = grouping1)) +
  geom_polygon(stat = 'ellipse', aes(fill = grouping1), alpha = 0.3) +
  #geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), #to include env_fit arrows
  #data = en_coord_cont, size = 1, colour = "darkred") +
  #geom_text(data = en_coord_cont, aes(x = PC1 + 0.05, y = PC2 + 0.10), colour = "darkred", size = 5, 
  #fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Depth") +
  #geom_segment(data = top5_loadings, aes(x=0,y=0, xend = (PC1*100), yend = (PC2*100)), #to include env_fit arrows
  #color = 'black', size = 1) +
  #geom_text(data = top5_loadings, aes(x = PC1 * 100, y = PC2 * 100), size = 5, 
  #fontface = "bold", label = row.names(top5_loadings)) +
  xlab('PC1 (explained var. 11.06%)') +
  ylab('PC2 (explained var. 5.62%)')

gg

# PERMANOVA to investigate difference in composition across sediment depths ####
set.seed(5678)
otu_tab <- arch_ra[order(rownames(arch_ra)),]
full_meta <- full_meta[order(rownames(full_meta)),]
rownames(otu_tab) == rownames(full_meta) # ensure that formatting of metadata table and otu table is the same

# full model, with euclidean distances which is equivalent to robust aitchison distances
# depth modelled as a categorical variable so that stratification can be observed
depth_model <- adonis2(otu_tab ~ as.factor(full_meta$Depth) + full_meta$Site, strata = full_meta$Core, method = 'euclidean')

# pairwise with fdr corrected
depth_model_pairwise <- pairwise.adonis(otu_tab, factors = full_meta$Depth, sim.method  = 'euclidean', p.adjust.m = 'fdr')
depth_model_pairwise[order(depth_model_pairwise$p.adjusted),]

write.csv(depth_model_pairwise, './Output/archaea_pairwise_permanova.csv')

# betadisper to assess homogeneity of group dispersion
bd_arch <- betadisper(vegdist(otu_tab,'euclidean'), group = full_meta$Depth, type = 'centroid')
anova(bd_arch)
plot(bd_arch)
boxplot(bd_arch)
permutest(bd_arch)

# which pairs are significantly different
arch_pair <- TukeyHSD(bd_arch, p.adjust = 'fdr')$group
arch_pair %>% as.data.frame() %>% filter("p adj" < 0.05) # no depth pairs with significantly different group dispersion

write.csv(anova(bd_arch), './Output/betadisper_results_arch.csv')
write.csv(arch_pair, './Output/betadisper_pairwise_arch.csv')


grouping1_pairwise <- pairwise.adonis(otu_tab, factors = full_meta$grouping1, sim.method = 'euclidean')

# corncob for differential analysis ####
# devtools::install_github("bryandmartin/corncob")
library(corncob)
ps_arch@sam_data <- sample_data(full_meta)
ps_class <- ps_arch %>% tax_glom('Class')
ps_class_sub <- ps_class %>% subset_samples(grouping1 != '10-30') # change to see different pairs

corncob_arch_2v3 <- differentialTest(formula = ~ grouping1, # change name to create two comparison groups
                                    phi.formula = ~ grouping1,
                                    formula_null = ~ 1,
                                    phi.formula_null = ~ grouping1,
                                    test = "Wald", boot = FALSE,
                                    data = ps_class_sub,
                                    fdr_cutoff = 0.05)


# what are the taxa that are different between surface vs subsurface, and subsurface vs deep=
overlap_taxa <- corncob_arch_1v2$significant_taxa[corncob_arch_1v2$significant_taxa %in% corncob_arch_2v3$significant_taxa]
ps_class@tax_table[overlap_taxa,'Class']

plot(corncob_arch_1v2, level = 'Class')
plot(corncob_arch_2v3, level = 'Class')

saveRDS(corncob_arch_1v2, './Output/corncob_arch_1v2.RDS')
saveRDS(corncob_arch_2v3, './Output/corncob_arch_2v3.RDS')

corncob$p_fdr %>% sort()
corncob$restrictions_DA
corncob$discriminant_taxa_DA

# arch db-RDA ####
mdf_meta = ps_arch@sam_data %>% data.frame()
rownames(arch_ra) == rownames(mdf_meta)

mdf_meta_complete <- mdf_meta[mdf_meta$Mean_PS !='NA',]
mdf_com_complete <- arch_ra[mdf_meta$Mean_PS !='NA',]

mdf_meta_complete <- mutate_all(mdf_meta_complete, function(x) as.numeric(as.character(x)))

mdf_meta_scale <- mdf_meta_complete[,5:12] %>% scale() %>% data.frame() %>% select(-c('X.H', 'Median_PS'))

arch_rda <- rda(mdf_com_complete ~ mdf_meta_scale$pH + mdf_meta_scale$P + mdf_meta_scale$X.N +  
                  mdf_meta_scale$X.S + mdf_meta_scale$X.C + mdf_meta_scale$Mean_PS,
               distance = 'euclidean')

overall<-anova(arch_rda) #overall model significance
axis<-anova(arch_rda, by="axis") #test axes for significance
terms<-anova(arch_rda, by="terms") #test for significant environmental variables

arch_rda_summary <- summary(arch_rda)
arch_rda_summary$concont #proportion of the variance explained by these explanatory variables (e.g. 55.33% (fitted) of 9.159%)
arch_rda_summary$cont # actual proportion of variance explained (total variation)

arch_rda_summary$biplot #coefficient associated with each explanatory variable
arch_rda_summary$sites #used to plot dbRDA

rda_df <- arch_rda_summary$sites %>% as.data.frame()
rda_df <- rda_df[order(rownames(rda_df)),]
mdf_meta_complete <- mdf_meta_complete[order(mdf_meta_complete$Depth),]
mdf_meta_complete$Depth <- as.numeric(mdf_meta_complete$Depth)
mdf_meta_complete %>% count(Depth)
mdf_meta_complete$grouping1 <- c(rep('10-30', 86), rep('40-60', 65), rep('70-100',39))

mdf_meta_complete <- mdf_meta_complete[order(rownames(mdf_meta_complete)),]
rownames(mdf_meta_complete) == rownames(rda_df)


rda_df$grouping1 <- mdf_meta_complete$grouping1

sig.en <- terms %>% as.data.frame() %>% filter(`Pr(>F)` < 0.05) %>% row.names()
en_arrows <- arch_rda_summary$biplot %>% as.data.frame() %>% .[sig.en,]
en_arrows$names <- rownames(en_arrows)

ggplot() +
  geom_point(data = rda_df, aes(x = RDA1, y = RDA2, colour = grouping1)) +
  geom_segment(data = en_arrows, aes(x = 0, y = 0, xend = RDA1*2.5, yend = RDA2*2.5), colour = '#898989') +
  geom_text(data = en_arrows, aes(x = RDA1*2.5, y= RDA2*2.5), label = row.names(en_arrows)) +
  theme_classic()

plot(arch_rda)
