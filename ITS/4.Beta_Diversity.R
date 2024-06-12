##  ###################################################  ##
##  Beta-diversity measures - btwn site and depth        ##
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
library(pairwiseAdonis); packageVersion("pairwiseAdonis")
library(ggfortify)
library(metagMisc)


# read in phyloseq, filter, format, and transform to robust clr ####
# threshold: prevalence of > 1% AND reads > 5
filt_list <- (microbiome::prevalence(ps) > 0.01) & (microbiome::abundances(ps) %>% rowSums() > 5)

ps_filtered <- readRDS('./Output/noncontam_ps_object.RDS') %>%
  subset_samples(SampleID != 'UL7_30') %>%
  subset_samples(Depth != 'NA') %>%
  subset_taxa(filt_list) %>%
  subset_taxa(Class != 'NA')

mdf_com <- ps_filtered %>% otu_table() %>% decostand(method = 'rclr')

# metadata file
full_meta = ps_filtered@sam_data %>% data.frame()
mdf_meta <- full_meta[4:11] %>% mutate_all(., function(x) {as.numeric(x)})

# set depth groupings
full_meta$Depth = as.numeric(full_meta$Depth)
full_meta <- full_meta[order(full_meta$Depth),]
full_meta %>% count(Depth)
full_meta$grouping1 <- c(rep('10-30',87),
                         rep('40-60',66),
                         rep('70-100',39))


#### PCA ####
otu_res <- prcomp(mdf_com, center = T)
en <- envfit(otu_res, mdf_meta, permutations = 999, na.rm = T)

ordiplot(otu_res)
plot(en)

summary(otu_res)
plot(otu_res)

data.scores = as.data.frame(scores(otu_res, display =  'sites'))
data.scores$Depth = full_meta$Depth
data.scores$Site = full_meta$Site

data.scores <- data.scores[order(rownames(data.scores)),]
full_meta <- full_meta[order(rownames(full_meta)),]
rownames(data.scores) == rownames(full_meta)

data.scores$grouping1 = full_meta$grouping1 # visualize depth groupings

PCAloadings <- data.frame(Variables = rownames(otu_res$rotation), otu_res$rotation)
PCAloadings <- PCAloadings %>% dplyr::select(PC1, PC2)
PCAloadings$strength <- sqrt((PCAloadings$PC1)^2 + (PCAloadings$PC2)^2)
PCAloadings <- PCAloadings[order(PCAloadings$strength, decreasing = TRUE),]

top5_loadings <- PCAloadings[1:5,]

gg = ggplot(data = data.scores, aes(x = PC1, y = PC2)) + 
  geom_point(data = data.scores, aes(colour = grouping1), size = 3, alpha = 1) +
  geom_polygon(stat = 'ellipse', aes(fill = grouping1), alpha = 0.3) +
  #geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               #data = en_coord_cont, size = 1, colour = "darkred") +
  #geom_text(data = en_coord_cont, aes(x = PC1, y = PC2), colour = "darkred", size = 4, 
            #fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Depth") +
  xlab('PC1 (explained var. 9.65%)') +
  ylab('PC2 (explained var. 3.79%)')

gg

ggsave('./Output/Figs/pca_plot.pdf', height = 12, width = 6)

#PERMANOVA to investigate difference in composition across sediment depths ####
set.seed(7890)
full_meta <- full_meta[order(rownames(full_meta)),]
rownames(mdf_com) == rownames(full_meta) # ensure that formatting of metadata table and otu table is the same

# full model, with euclidean distances which is equivalent to robust aitchison distances
# depth modelled as a categorical variable so that stratification can be observed
depth_model <- adonis2(mdf_com ~ as.factor(full_meta$Depth) + full_meta$Site, strata = full_meta$Core, method = 'euclidean')

# pairwise with fdr corrected
depth_model_pairwise <- pairwise.adonis(mdf_com, factors = full_meta$Depth, sim.method  = 'euclidean', p.adjust.m = 'fdr')
depth_model_pairwise[order(depth_model_pairwise$p.adjusted),]

write.csv(depth_model_pairwise, './Output/permanova_pairwise.csv')

# betadisper to assess homogeneity of group dispersion
bd <- betadisper(vegdist(mdf_com, 'euclidean'), group = full_meta$Depth)
anova(bd)
plot(bd)
boxplot(bd)
permutest(bd)

# which pairs are significantly different
bd_pair <- TukeyHSD(bd, p.adjust = 'fdr')$group
bd_pair %>% as.data.frame() %>% filter("p adj" < 0.05)

write.csv(anova(bd), './Output/betadisper_results.csv')
write.csv(bd_pair, './Output/betadisper_pairwise.csv')

# corncob for differential analysis ####
# devtools::install_github("bryandmartin/corncob")
library(corncob)
ps_genus@sam_data <- sample_data(full_meta)
ps_class <- ps_genus %>% tax_glom('Class')
ps_class_sub <- ps_class %>% subset_samples(grouping1 != '10-30') # change to see different pairs

corncob_genus_2v3 <- differentialTest(formula = ~ grouping1, # change name to create two comparison groups
                                    phi.formula = ~ grouping1,
                                    formula_null = ~ 1,
                                    phi.formula_null = ~ grouping1,
                                    test = "Wald", boot = FALSE,
                                    data = ps_class_sub,
                                    fdr_cutoff = 0.05)

# what are the taxa that are different between surface vs subsurface, and subsurface vs deep
overlap_taxa <- corncob_genus_1v2$significant_taxa[corncob_genus_1v2$significant_taxa %in% corncob_genus_2v3$significant_taxa]
ps_class@tax_table[overlap_taxa,'Class']

plot(corncob_genus_1v2, level = 'Class')
plot(corncob_genus_2v3, level = 'Class')

saveRDS(corncob_genus_1v2, './Output/corncob_genus_1v2.RDS')
saveRDS(corncob_genus_2v3, './Output/corncob_genus_2v3.RDS')

corncob$p_fdr %>% sort()
corncob$restrictions_DA
corncob$discriminant_taxa_DA

# db-RDA ####
mdf_meta_complete <- mdf_meta[complete.cases(mdf_meta),]
mdf_com_complete <- mdf_com[complete.cases(mdf_meta),]

rownames(mdf_meta_complete) == rownames(mdf_com_complete)

mdf_meta_scale <- mdf_meta_complete[,2:8] %>% scale() %>% data.frame() %>% select(-c('X.H'))

fung_rda <- rda(mdf_com_complete ~ mdf_meta_scale$pH + mdf_meta_scale$P + mdf_meta_scale$X.N + mdf_meta_scale$X.S + mdf_meta_scale$Mean_PS + mdf_meta_scale$X.C,
      distance = 'euclidean')

overall<-anova(fung_rda) #overall model significance
axis<-anova(fung_rda, by="axis") #test axes for significance
terms<-anova(fung_rda, by="terms") #test for significant environmental variables

fung_rda_summary <- summary(fung_rda)
fung_rda_summary$concont #proportion of the variance explained by these explanatory variables (e.g. 55.33% (fitted) of 9.159%)
fung_rda_summary$cont # actual proportion of variance explained (total variation)

fung_rda_summary$biplot #coefficient associated with each explanatory variable
fung_rda_summary$constraints #used to plot dbRDA

rda_df <- fung_rda_summary$constraints %>% as.data.frame()
rda_df <- rda_df[order(rownames(rda_df)),]
mdf_meta_complete <- mdf_meta_complete[order(mdf_meta_complete$Depth),]
mdf_meta_complete %>% count(Depth)
mdf_meta_complete$grouping1 <- c(rep('10-30',85),
                                 rep('40-60',65),
                                 rep('70-100',39))

mdf_meta_complete <- mdf_meta_complete[order(rownames(mdf_meta_complete)),]
rownames(mdf_meta_complete) == rownames(rda_df)

rda_df$grouping1 <- mdf_meta_complete$grouping1

sig.en <- terms %>% as.data.frame() %>% filter(`Pr(>F)` <0.05) %>% row.names()
en_arrows <- fung_rda_summary$biplot %>% as.data.frame() %>% .[sig.en,]
en_arrows$names <- rownames(en_arrows)

ggplot() +
  geom_point(data = rda_df, aes(x = RDA1, y = RDA2, colour = grouping1)) +
  geom_segment(data = en_arrows, aes(x = 0, y = 0, xend = RDA1*3, yend = RDA2*3), colour = '#898989') +
  geom_text(data = en_arrows, aes(x = RDA1*3, y= RDA2*3), label = row.names(en_arrows)) +
  theme_classic()

plot(fung_rda)
