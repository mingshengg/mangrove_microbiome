##  ###################################################  ##
##  Beta-diversity measures - btwn site and depths       ##
##                                                       ##
##  Author: Ming Sheng - 16 Nov, 2022                    ##
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
library(metagMisc)

# read in phyloseq, filter, format, and transform to robust clr ####
ps <- readRDS('./Output/noncontam_ps_object.RDS')
filt_list <- (microbiome::prevalence(ps) > 0.01) & (microbiome::abundances(ps) %>% rowSums() > 5)

# filter to bacteria and to set threshold
# threshold: prevalence of > 1% AND reads > 5
ps_filtered <- ps %>% 
  subset_samples(Depth != 'NA') %>%
  subset_taxa(filt_list) %>%
  subset_taxa(Class != 'NA')

# subset to bacteria group
ps_bac <- ps_filtered %>% subset_taxa(Kingdom == 'Bacteria') %>% subset_taxa(Class != 'NA') %>% subset_samples(Depth != 'NA')

# transform to robust clr
bac_ra <- ps_bac@otu_table %>% decostand(method = 'rclr')

# metadata file
full_meta <- ps_bac@sam_data %>% data.frame()
mdf_meta <- full_meta[4:11] %>% mutate_all(., function(x) {as.numeric(x)})

# set depth groupings
full_meta$Depth <- as.numeric(full_meta$Depth)
full_meta <- full_meta[order(full_meta$Depth),]
full_meta %>% count(Depth)
full_meta$grouping1 <- c(rep('10-30', 88), rep('40-60', 66), rep('70-100',39))

ps_bac@sam_data <- ps_bac@sam_data[order(ps_bac@sam_data$Depth),]
ps_bac@sam_data$grouping1 <- c(rep('10-30', 88), rep('40-60', 66), rep('70-100',39))

#### PCA ####
otu_res <- prcomp(bac_ra, center = T)

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

gg = ggplot(data = data.scores, aes(x = PC1, y = PC2)) + 
  geom_point(data = data.scores, aes(colour = grouping1)) +
  geom_polygon(stat = 'ellipse', aes(fill = grouping1), alpha = 0.3) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Depth") +
  xlab('PC1 (explained var. 5.06%)') +
  ylab('PC2 (explained var. 4.44%)')

gg

### PERMANOVA ###
set.seed(18005)
otu_tab <- bac_ra[order(rownames(bac_ra)),]
full_meta <- full_meta[order(rownames(full_meta)),]
rownames(otu_tab) == rownames(full_meta)

# full model, with euclidean distances which is equivalent to robust aitchison distances
# depth modelled as a categorical variable so that stratification can be observed
depth_model <- adonis2(otu_tab ~ as.factor(full_meta$Depth) + full_meta$Site, strata = full_meta$Core, method = 'euclidean')

# pairwise with fdr corrected
depth_model_pairwise <- pairwise.adonis(otu_tab, factors = full_meta$Depth, sim.method  = 'euclidean', p.adjust.m = 'fdr')
depth_model_pairwise[order(depth_model_pairwise$p.adjusted),]

write.csv(depth_model_pairwise, './Output/bacteria_pairwise_permanova.csv')

# betadisper to visualize homogeneity across groups 
bd_bac <- betadisper(vegdist(otu_tab, 'euclidean'), group = full_meta$grouping1)
anova(bd_bac)
boxplot(bd_bac)
bac_pair <- TukeyHSD(bd_bac, p.adjust = 'fdr')$group

write.csv(anova(bd_bac), './Output/betadisper_results_bac.csv')
write.csv(bac_pair, './Output/betadisper_pairwise_bac.csv')


# corncob for differential analysis ####
# devtools::install_github("bryandmartin/corncob")
library(corncob)
ps_bac@sam_data <- sample_data(full_meta)
ps_class <- ps_bac %>% tax_glom('Class')
ps_class_sub <- ps_class %>% subset_samples(grouping1 != '40-60') #  change to see different pairs

corncob_bac_1v2 <- differentialTest(formula = ~ grouping1, # change name to create two comparison groups
                                    phi.formula = ~ grouping1,
                                    formula_null = ~ 1,
                                    phi.formula_null = ~ grouping1,
                                    test = "Wald", boot = FALSE,
                                    data = ps_class_sub,
                                    fdr_cutoff = 0.05)

ps_class_sub <- ps_class %>% subset_samples(grouping1 != '10-30')

corncob_bac_2v3 <- differentialTest(formula = ~ grouping1,
                            phi.formula = ~ grouping1,
                            formula_null = ~ 1,
                            phi.formula_null = ~ grouping1,
                            test = "Wald", boot = FALSE,
                            data = ps_class_sub,
                            fdr_cutoff = 0.05)


# what are the taxa that are different between surface vs subsurface, and subsurface vs deep
overlap_taxa <- corncob_bac_1v2$significant_taxa[corncob_bac_1v2$significant_taxa %in% corncob_bac_2v3$significant_taxa]
ps_class@tax_table[overlap_taxa,'Class']

plot(corncob_bac_1v2, level = 'Class')
plot(corncob_bac_2v3, level = 'Class')

saveRDS(corncob_bac_1v2, './Output/corncob_bac_1v2.RDS')
saveRDS(corncob_bac_2v3, './Output/corncob_bac_2v3.RDS')

corncob$p_fdr %>% sort()
corncob$restrictions_DA
corncob$discriminant_taxa_DA

# bac db-RDA ####
mdf_meta = ps_bac@sam_data %>% data.frame()

bac_ra <- bac_ra[order(rownames(bac_ra)),]
mdf_meta <- mdf_meta[order(rownames(mdf_meta)),]

rownames(bac_ra) == rownames(mdf_meta)


mdf_meta_complete <- mdf_meta[mdf_meta$Mean_PS !='NA',]
mdf_com_complete <- bac_ra[mdf_meta$Mean_PS !='NA',]

mdf_meta_complete <- mutate_all(mdf_meta_complete[,4:12], function(x) as.numeric(as.character(x)))

mdf_meta_scale <- mdf_meta_complete[,2:9] %>% scale() %>% data.frame() %>% select(-c('X.H','Median_PS'))

corrplot::corrplot.mixed(cor(mdf_meta_scale,  method = 'spearman'))

bac_rda <- rda(mdf_com_complete ~ mdf_meta_scale$pH + mdf_meta_scale$P + mdf_meta_scale$X.N + 
                 mdf_meta_scale$X.S + mdf_meta_scale$X.C + mdf_meta_scale$Mean_PS,
                distance = 'euclidean')

overall <- anova(bac_rda) #overall model significance
axis <- anova(bac_rda, by="axis") #test axes for significance
terms <- anova(bac_rda, by="terms") #test for significant environmental variables

bac_rda_summary <- summary(bac_rda)
bac_rda_summary$concont #proportion of the variance explained by these explanatory variables (e.g. 55.33% (fitted) of 9.159%)
bac_rda_summary$cont # actual proportion of variance explained (total variation)

bac_rda_summary$biplot #coefficient associated with each explanatory variable
bac_rda_summary$sites #used to plot dbRDA

rda_df <- bac_rda_summary$sites %>% as.data.frame() # dataframe for plotting db-rda
rda_df <- rda_df[order(rownames(rda_df)),] 
mdf_meta_complete <- mdf_meta_complete[order(mdf_meta_complete$Depth),]
mdf_meta_complete$Depth <- as.numeric(mdf_meta_complete$Depth)
mdf_meta_complete %>% count(Depth)
mdf_meta_complete$grouping1 <- c(rep('10-30', 86), rep('40-60', 65), rep('70-100',39))

mdf_meta_complete <- mdf_meta_complete[order(rownames(mdf_meta_complete)),]
rownames(mdf_meta_complete) == rownames(rda_df) # ensure rownames are the same before plotting


rda_df$grouping1 <- mdf_meta_complete$grouping1

# adding only significant variables
sig.en <- terms %>% as.data.frame() %>% filter(`Pr(>F)` < 0.05) %>% row.names() 
en_arrows <- bac_rda_summary$biplot %>% as.data.frame() %>% .[sig.en,]
en_arrows$names <- rownames(en_arrows)

ggplot() +
  geom_point(data = rda_df, aes(x = RDA1, y = RDA2, colour = grouping1)) +
  geom_segment(data = en_arrows, aes(x = 0, y = 0, xend = RDA1*2.5, yend = RDA2*2.5), colour = '#898989') +
  geom_text(data = en_arrows, aes(x = RDA1*2.5, y= RDA2*2.5), label = row.names(en_arrows)) +
  theme_classic()

plot(bac_rda)
