##  ###################################################  ##
##  Alpha-diversity measures - btwn site and structure   ##
##                                                       ##
##  Author: Ming Sheng - June 17, 2022                   ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  broom v 0.8.0                                        ##
##  patchwork v 1.1.1                                    ##
##  microbiome v 1.18.0                                  ##
##  emmeans v 1.7.4.1                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages, data, and customizations ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(patchwork); packageVersion("patchwork")
library(microbiome); packageVersion("microbiome")
library(broom); packageVersion("broom")
library(lme4); packageVersion("lme4")
library(lmerTest); packageVersion("lmerTest")
library(emmeans)
library(metagMisc)
source("./plot_bar2.R")

se <- function(x) sd(x)/sqrt(length(x))

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#FAEBD7")

# Load ps objects
filt_list <- (microbiome::prevalence(ps) > 0.01) & (microbiome::abundances(ps) %>% rowSums() > 5)

ps_filtered <- readRDS('./Output/noncontam_ps_object.RDS') %>% 
  subset_samples(Depth != 'NA') %>%
  subset_taxa(filt_list) %>%
  subset_taxa(Class != 'NA') %>%
  subset_samples(SampleID != 'UL7_30')

# Model alpha diversity ####
meta <- microbiome::meta(ps_genus)                
meta$Shannon <- vegan::diversity(otu_table(ps_genus),index = "shannon")
meta$Richness <- vegan::specnumber(otu_table(ps_genus))
meta$Evenness <- microbiome::evenness(otu_table(ps_genus), index = "simpson") %>% .$simpson


# add to ps object
ps_genus@sam_data$Richness <- meta$Richness
ps_genus@sam_data$Shannon <- meta$Shannon
ps_genus@sam_data$Evenness <- meta$Evenness


# lme models
library(nlme)
shannon_mod <- lme(data = meta,
                   Shannon ~ as.numeric(Depth) + Site, random = ~1|Core)
richness_mod <- lme(data = meta,
                   Richness ~ as.numeric(Depth) + Site, random = ~1|Core)
evenness_mod <- lme(data = meta,
                    Evenness ~ as.numeric(Depth) + Site, random = ~1|Core)

# send to file
sink("./Output/Stats/Shannon_Diversity_Model.txt")
anova(shannon_mod)
plot(shannon_mod)
sink(NULL)

meta %>% group_by(Depth) %>%
  summarise(shannon_mean = mean(Shannon), shannon_se = se(Shannon))

sink("./Output/Stats/Richness_Model.txt")
anova(richness_mod)
plot(richness_mod)
sink(NULL)

meta %>% group_by(Depth) %>%
  summarise(richness_mean = mean(Richness), shannon_se = se(Richness))


sink("./Output/Stats/Evenness_Model.txt")
anova(evenness_mod)
plot(evenness_mod)
sink(NULL)

meta %>% group_by(Depth) %>%
  summarise(evenness_mean = mean(Evenness), Evenness_se = se(Evenness))

emmeans(shannon_mod, ~ Site) %>% pairs()
emmeans(richness_mod, ~ Site) %>% pairs()
emmeans(evenness_mod, ~ Site) %>% pairs()

meta %>% group_by(Site) %>%
  summarise(Shannon_mean = mean(Shannon), Evenness_se = se(Shannon))
meta %>% group_by(Site) %>%
  summarise(Richness_mean = mean(Richness), Evenness_se = se(Richness))
meta %>% group_by(Site) %>%
  summarise(Evenness_mean = mean(Evenness), Evenness_se = se(Evenness))
# violin plots for shannon, richness, and evenness by depth 
meta$Depth <- factor(meta$Depth, levels = c('10','20','30','40','50','60','70','80','90','100'))

cf <- coef(lm(Shannon ~ as.numeric(Depth), data = meta))
shannon_plot <- meta %>% 
  ggplot(aes(x = Depth, y = Shannon, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(c(0,8)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none") +
  geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8)

shannon_plot

ggsave("./Output/Figs/shannon_by_depth.pdf",dpi=400,width = 12,height = 6)

cf <- coef(lm(Richness ~ as.numeric(Depth), data = meta))
richness_plot <- meta %>% 
  ggplot(aes(x = Depth, y = Richness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width  = 0.3, alpha = 0.3) +
  theme_classic() +
  ylim(c(0, 1500)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none") +
  geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8)

richness_plot

ggsave("./Output/Figs/richness_by_depth.pdf",dpi=400,width = 12,height = 6)

richness_mean <- meta %>% group_by(Depth) %>% summarise(mean_richness = mean(Richness)) %>%
  ggplot(aes(x = Depth, y = mean_richness, fill = Depth)) + 
  geom_point(size = 3) +
  geom_point(data = meta, aes(x = Depth, y = Richness), colour = 'grey') +
  theme_classic() +
  ylim(c(0, 700)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none")

richness_mean

cf <- coef(lm(Evenness ~ as.numeric(Depth), data = meta))
evenness_plot <- meta %>%
  ggplot(aes(x = Depth, y = Evenness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) +
  ylim(c(0,1.0)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size= 20),
        legend.position = "none") +
  geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8)

evenness_plot

ggsave("./Output/Figs/evenness_by_depth.pdf",dpi=400,width = 12,height = 6)


evenness_mean <- meta %>% group_by(Depth) %>% summarise(mean_evenness = mean(Evenness)) %>%
  ggplot(aes(x = Depth, y = mean_evenness, fill = Depth)) + 
  geom_point(size = 3) +
  geom_point(data = meta, aes(x = Depth, y = Evenness), colour = 'grey') +
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none")

evenness_mean

## Take depth as continuous variable
shannon_main <- aov(data = meta,
                   formula = Shannon ~ as.numeric(Depth))
richness_main <- aov(data = meta,
                    formula = Richness ~ as.numeric(Depth))

shannon_predict <- data.frame(Shannon = predict(shannon_main, list(Depth = seq(10,100,10))), Depth = seq(10,100,10))

shannon_plot <- ggplot(data = meta, aes(x = as.numeric(Depth), y = Shannon, colour = Site)) + 
  geom_point() +
  stat_smooth(method = 'lm', linewidth = 1.3, alpha = 0.3, se = TRUE, colour = 'black') +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25)) +
  scale_x_continuous(name = 'Depth (cm)', breaks = seq(10,100,10)) +
  labs(y = 'Shannon diversity index')

shannon_plot

ggsave("./Output/Figs/shannon_by_depth_continuous.pdf",dpi=400,width = 12,height = 6)

richness_predict <- data.frame(Richness = (predict(richness_main, list(Depth = seq(10,100,10)))), Depth = seq(10,100,10))

richness_plot <- ggplot(data = meta, aes(x = as.numeric(Depth), y = Richness, colour = Site)) + 
  geom_point() +
  stat_smooth(method = 'lm', linewidth = 1.3, alpha = 0.3, se = TRUE, colour = 'black') +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25)) +
  scale_x_continuous(name = 'Depth (cm)', breaks = seq(10,100,10)) +
  labs(y = 'Richness')

richness_plot

ggsave("./Output/Figs/richness_by_depth_continuous.pdf",dpi=400,width = 12,height = 6)

# violin plots for shannon and richness by site
shannon_plot <- meta %>% 
  ggplot(aes(x = Site, y = Shannon, fill = Site)) + 
  geom_violin(scale = 'width') +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Shannon diversity index')

shannon_plot

ggsave("./Output/Figs/shannon_by_Site.pdf",dpi=400,width = 12,height = 6)

richness_plot <- meta %>% 
  ggplot(aes(x = Site, y = Richness, fill = Site)) + 
  geom_violin(scale = 'width') +
  geom_jitter(width  = 0.3, alpha = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none")

richness_plot

ggsave("./Output/Figs/richness_by_Site.pdf",dpi=400,width = 12,height = 6)

# lme models
shannon_mod <- aov(data = meta,
                   formula = Shannon ~ Site)
richness_mod <- aov(data = meta,
                    formula = Richness ~ Site)

emmeans(shannon_mod, pairwise ~ Site)
emmeans(richness_mod, pairwise ~ Site)

# send to file
sink("./Output/Stats/Shannon_Diversity_Model.txt")
summary(shannon_mod)
par(mfrow = c(2,2))
plot(shannon_mod)
sink(NULL)

sink("./Output/Stats/Richness_Model.txt")
summary(richness_mod)
par(mfrow = c(2,2))
plot(richness_mod)
sink(NULL)


### violin plots per site
violin_shannon_cj <- ggplot(aes(x = Depth, y = Shannon, fill = Depth), data = meta %>% filter(Site == 'Chek Jawa')) + 
  geom_violin(scale = 'width', trim = F) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Shannon diversity index')


violin_shannon_jd <- ggplot(aes(x = Depth, y = Shannon, fill = Depth), data = meta %>% filter(Site == 'Jalan Durian')) + 
  geom_violin(scale = 'width', trim = F) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Shannon diversity index')

violin_shannon_ul <- ggplot(aes(x = Depth, y = Shannon, fill = Depth), data = meta %>% filter(Site == 'Ubin Living Lab')) + 
  geom_violin(scale = 'width', trim = F) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Shannon diversity index')

### richness
violin_rich_cj <- ggplot(aes(x = Depth, y = Richness, fill = Depth), data = meta %>% filter(Site == 'Chek Jawa')) + 
  geom_violin(scale = 'width', trim = F) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(0,150) +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Richness diversity index')


violin_rich_jd <- ggplot(aes(x = Depth, y = Richness, fill = Depth), data = meta %>% filter(Site == 'Jalan Durian')) + 
  geom_violin(scale = 'width', trim = F) +
  ylim(0,250) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Richness diversity index')

violin_rich_ul <- ggplot(aes(x = Depth, y = Richness, fill = Depth), data = meta %>% filter(Site == 'Ubin Living Lab')) + 
  geom_violin(scale = 'width', trim = F) +
  ylim(0, 210) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  stat_summary(fun = mean, geom = 'point', shape = 18, size = 4) +
  theme(axis.text.x = element_text(face="bold",size=20),
        axis.text.y = element_text(face="bold",size=20),
        axis.title = element_text(face="bold",size=25),
        legend.position = "none") +
  labs(y = 'Richness diversity index')

