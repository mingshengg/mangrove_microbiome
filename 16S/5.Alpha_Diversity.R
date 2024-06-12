##  ###################################################  ##
##  Alpha-diversity measures - btwn site and depth       ##
##                                                       ##
##  Author: Ming Sheng Ng                                ##
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

se <- function (x) sd(x)/sqrt(length(x))

# custom palette
pal <- c("#d98416","#25802d","#664c13","#858585", "#FAEBD7")

# Load ps objects
ps <- readRDS("./Output/noncontam_ps_object_bac.RDS")

filt_list <- (microbiome::prevalence(ps) > 0.01) & (microbiome::abundances(ps) %>% rowSums() > 5)

ps_genus <- ps %>% 
  subset_samples(Depth != 'NA') %>%
  subset_taxa(filt_list) %>%
  subset_taxa(Class != 'NA')

# Model alpha diversity ####
ps_bact <- ps_genus %>% subset_taxa(Kingdom == 'Bacteria')

meta_bact <- ps_bact %>% microbiome::meta()
meta_bact$Shannon <- vegan::diversity(otu_table(ps_bact),index = "shannon")
meta_bact$Richness <- vegan::specnumber(otu_table(ps_bact))
meta_bact$Evenness <- microbiome::evenness(otu_table(ps_bact), index = "simpson") %>% .$simpson

ps_arch <- ps_genus %>% subset_taxa(Kingdom == 'Archaea')
meta_arch <- ps_arch %>% microbiome::meta()
meta_arch$Shannon <- vegan::diversity(otu_table(ps_arch),index = "shannon")
meta_arch$Richness <- vegan::specnumber(otu_table(ps_arch))
meta_arch$Evenness <- microbiome::evenness(otu_table(ps_arch), index = "simpson") %>% .$simpson


# add to ps object
ps_bact@sam_data$Richness <- meta_bact$Richness
ps_bact@sam_data$Shannon <- meta_bact$Shannon
ps_bact@sam_data$Evenness <- meta_bact$Evenness

ps_arch@sam_data$Richness <- meta_arch$Richness
ps_arch@sam_data$Shannon <- meta_arch$Shannon
ps_arch@sam_data$Evenness <- meta_arch$Evenness

# exploratory plots
plot(meta_bact$Depth, meta_bact$Shannon)
plot(meta_bact$Depth, meta_bact$Richness)

plot(meta_arch$Depth, meta_arch$Shannon) #looks linear
plot(meta_arch$Depth, meta_arch$Richness) #does not look linear

## bacteria
# violin plots for shannon and richness by depth 
meta_bact$Depth <- factor(meta_bact$Depth, levels = c('10','20','30','40','50','60','70','80','90','100'))

shannon_plot_bact <- meta_bact %>%
  ggplot(aes(x = Depth, y = Shannon, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(c(0,8)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size= 20),
        legend.position = "none")

shannon_plot_bact

ggsave("./Output/Figs/shannon_bact_by_depth.pdf",dpi=400,width = 12,height = 6)

richness_plot_bact <- meta_bact %>% 
  ggplot(aes(x = Depth, y = Richness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width  = 0.3, alpha = 0.3) +
  ylim(c(0, 1500)) +
  theme_classic() +
  theme(axis.text.x = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.title = element_text(face="bold",size=20),
        legend.position = "none")

richness_plot_bact

ggsave("./Output/Figs/richness_bact_by_depth.pdf",dpi=400,width = 12,height = 6)

evenness_plot_bact <- meta_bact %>%
  ggplot(aes(x = Depth, y = Evenness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size= 20),
        legend.position = "none")

evenness_plot_bact

ggsave("./Output/Figs/evenness_bact_by_depth.pdf",dpi=400,width = 12,height = 6)


## archaea
# violin plots for shannon and richness by depth 
cf <- coef(lm(Shannon ~ as.numeric(Depth), data = meta_arch))

meta_arch$Depth <- factor(meta_arch$Depth, levels = c('10','20','30','40','50','60','70','80','90','100'))
shannon_plot_arch <- meta_arch %>%
  ggplot(aes(x = Depth, y = Shannon, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(c(0,8)) +
  geom_smooth(aes(y=Shannon, x=as.numeric(as.character(Depth))),method = 'lm') +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none") +
  #geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8) + 
  facet_grid(~Site)

shannon_plot_arch

ggsave("./Output/Figs/shannon_arch_by_depth.pdf",dpi=400,width = 12,height = 6)

cf <- coef(lm(Richness ~ as.numeric(Depth), data = meta_arch))
richness_plot_arch <- meta_arch %>% 
  ggplot(aes(x = Depth, y = Richness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width  = 0.3, alpha = 0.3) +
  ylim(c(0,1500)) +
  theme_classic() +
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size=20),
        legend.position = "none") +
  geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8) +

richness_plot_arch
ggsave("./Output/Figs/richness_arch_by_depth.pdf",dpi=400,width = 12,height = 6)

cf <- coef(lm(Evenness ~ as.numeric(Depth), data = meta_arch))
evenness_plot_arch <- meta_arch %>%
  ggplot(aes(x = Depth, y = Evenness, fill = Depth)) + 
  geom_violin(scale = 'width', trim = FALSE) +
  geom_jitter(width = 0.3, alpha = 0.3) + 
  theme_classic() +
  ylim(c(0,1)) + 
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title = element_text(size= 20),
        legend.position = "none") +
  geom_segment(aes(x = 1, xend = 10, y = cf[1] + cf[2], yend = cf[1] + 9*cf[2]), lwd = .8)

evenness_plot_arch

ggsave("./Output/Figs/evenness_arch_by_depth.pdf",dpi=400,width = 12,height = 6)

# Bacteria
# Shannon models
library(nlme)
shannon_mod <- lme(data = meta_bact,
                   Shannon ~ as.numeric(Depth) + Site, random = ~ 1|Core) #site not significant

# send to file
sink("./Output/Shannon_Diversity_Model_bact.txt")
anova(shannon_mod)
summary(shannon_mod)
plot(shannon_mod)
sink(NULL)

mean(meta_bact$Shannon)
se(meta_bact$Shannon)

#Richness models
richness_mod <- lme(data = meta_bact,
                    log(Richness) ~ as.numeric(Depth) + Site, random = ~ 1|Core)

sink("./Output/Richness_Model_bact.txt")
anova(richness_mod)
summary(richness_mod)
plot(richness_mod)
sink(NULL)

mean(meta_bact$Richness)
se(meta_bact$Richness)
#Evenness models
evenness_mod <- lme(data = meta_bact,
                    (Evenness) ~ as.numeric(Depth) + Site, weights = ~as.numeric(Depth)^2, random = ~ 1|Core)

qqnorm(evenness_mod, ~ranef(., level=1))

qqnorm(resid(evenness_mod))
qqline(resid(evenness_mod))


sink("./Output/Evenness_Model_bact.txt")
anova(evenness_mod)
summary(evenness_mod)
plot(evenness_mod)
sink(NULL)

mean(meta_bact$Evenness)
se(meta_bact$Evenness)

# Archaea
# Shannon models
shannon_mod <- lme(data = meta_arch,
                   Shannon ~ as.numeric(Depth) + Site, random = ~ 1|Core)

# send to file
sink("./Output/Shannon_Diversity_Model_arch.txt")
anova(shannon_mod)
summary(shannon_mod)
plot(shannon_mod)
sink(NULL)

meta_arch %>% group_by(Depth) %>%
  summarise(Shannon_depth = mean(Shannon), shannon_se = se(Shannon))

#Richness models
richness_mod <- lme(data = meta_arch,
                    Richness ~ as.numeric(Depth) + Site , random = ~ 1|Core)

sink("./Output/Richness_Model_arch.txt")
anova(richness_mod)
summary(richness_mod)
plot(richness_mod)
sink(NULL)

meta_arch %>% group_by(Depth) %>%
  summarise(Richness_depth = mean(Richness), Richness_se = se(Richness))

#Evenness models
evenness_mod <- lme(data = meta_arch,
                    (Evenness) ~ as.numeric(Depth) + Site, random = ~ 1|Core)

qqnorm(evenness_mod, ~ranef(., level=1))
qqnorm(resid(evenness_mod))

sink("./Output/evenness_Model_arch.txt")
anova(evenness_mod)
summary(evenness_mod)
plot(evenness_mod)
sink(NULL)

meta_arch %>% group_by(Depth) %>%
  summarise(Evenness_depth = mean(Evenness), Evenness_se = se(Evenness))
