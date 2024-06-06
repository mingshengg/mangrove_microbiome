##  ###################################################  ##
##  Investigate and tidy up full phyloseq object         ##
##                                                       ##
##  Author: Ming Sheng - 14 Nov, 2022                    ##
##                                                       ##
##  Software versions:                                   ##
##  R v 4.2.0                                            ##
##  tidyverse v 1.3.1                                    ##
##  phyloseq v 1.40.0                                    ##
##  vegan v 2.6.2                                        ##
##  VennDiagram v 1.7.3                                  ##
##  patchwork v 1.1.1                                    ##
##                                                       ##
##  ###################################################  ##

# Load packages ####
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(vegan); packageVersion("vegan")
library(VennDiagram); packageVersion("VennDiagram")
library(patchwork); packageVersion("patchwork")
library(ggplot2)

`%notin%` <- Negate(`%in%`)

# custom palette
pal <- c("#d98416","#25802d","#664c13")
# names(pal) <- c("f","l","p","s")

# Load ps object with tree ####
full_ps <- readRDS("./Output/noncontam_ps_object.RDS")

# agglomerate taxa at genus level ####
full_ps_genus <- tax_glom(full_ps,"Genus")
full_ps_genus <- full_ps_genus %>% subset_taxa(taxa_sums(full_ps_genus)>0)
full_ps_genus <- full_ps_genus %>% subset_samples(sample_sums(full_ps_genus)>0)
full_ps_genus <- full_ps_genus %>% subset_samples(Site != "NA")

# output genus-level ps_object for convenience
saveRDS(full_ps_genus, "./Output/full_ps_object_genus-glom.RDS")
full_ps_genus <- readRDS("./Output/full_ps_object_genus-glom.RDS")

# VennDiagram of overlap for Species ####

A <- full_ps_genus %>% subset_samples(Site == "Chek Jawa") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Chek Jawa"))>0) %>% tax_table() %>% row.names()
B <- full_ps_genus %>% subset_samples(Site == "Ubin Living Lab") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Ubin Living Lab"))>0) %>% tax_table() %>% row.names()
C <- full_ps_genus %>% subset_samples(Site == "Jalan Durian") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Site == "Jalan Durian"))>0) %>% tax_table() %>% row.names()


full <- unique(c(A,B,C))
n12 <- sum(full %in% unique(A) & full %in% unique(B))
n23 <- sum(full %in% unique(B) & full %in% unique(C))
n13 <- sum(full %in% unique(A) & full %in% unique(C))
n123 <- sum(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))

venn.plot1 <- draw.triple.venn(length((A)),length((B)),length((C)),
                         n12,n23,n13,n123,
                         category = c('Chek Jawa','Ubin Living Lab','Jalan Durian'),
                         fill = pal)

dev.off()
png("./Output/Figs/VennDiagram_Shared_Genus-Level_Taxa_by_Site.png")
grid.draw(venn.plot1)
dev.off()

# Investigate shared species across
full_ps <- readRDS('./Output/clean_noncontam_ps_object.RDS')

A <- full_ps %>% subset_samples(Site == "Chek Jawa") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Chek Jawa"))>0) %>% tax_table() %>% row.names()
B <- full_ps %>% subset_samples(Site == "Ubin Living Lab") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Ubin Living Lab"))>0) %>% tax_table() %>% row.names()
C <- full_ps %>% subset_samples(Site == "Jalan Durian") %>% 
  subset_taxa(taxa_sums(full_ps %>% subset_samples(Site == "Jalan Durian"))>0) %>% tax_table() %>% row.names()

full <- unique(c(A,B,C))

shared_species <- full %>% subset(full %in% unique(A) & full %in% unique(B) & full %in% unique(C))

ps_comp <- full_ps %>% microbiome::transform('compositional')

species_shared_otu_table <- ps_comp %>% subset_samples(PCR_Negative == FALSE) %>% 
  .@otu_table %>% as.data.frame() %>% select(shared_species)

species_shared_tax_table <- ps_comp@tax_table %>% as.data.frame() %>% filter(row.names(ps_comp@tax_table) %in% shared_species)

species_max_otu <- species_shared_otu_table %>% apply(2, max) %>% sort(decreasing = T) 

species_top_5_otu <- head(species_max_otu, 5)

top_5_species <- species_shared_tax_table %>% filter(row.names(species_shared_tax_table) %in% names(species_top_5_otu))

species_mean_otu <- species_shared_otu_table %>% apply(2, mean) %>% sort(decreasing = T) %>% head(5)

mean_5_species <- species_shared_tax_table %>% filter(row.names(species_shared_tax_table) %in% names(species_mean_otu))

# Number of unique ASVs per depth

n1 <- full_ps_genus %>% subset_samples(Depth == "10") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "10"))>0) %>% tax_table() %>% row.names()

n2 <- full_ps_genus %>% subset_samples(Depth == "20") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "20"))>0) %>% tax_table() %>% row.names()

n3 <- full_ps_genus %>% subset_samples(Depth == "30") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "30"))>0) %>% tax_table() %>% row.names()

n4 <- full_ps_genus %>% subset_samples(Depth == "40") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "40"))>0) %>% tax_table() %>% row.names()

n5 <- full_ps_genus %>% subset_samples(Depth == "50") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "50"))>0) %>% tax_table() %>% row.names()

n6 <- full_ps_genus %>% subset_samples(Depth == "60") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "60"))>0) %>% tax_table() %>% row.names()

n7 <- full_ps_genus %>% subset_samples(Depth == "70") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "70"))>0) %>% tax_table() %>% row.names()

n8 <- full_ps_genus %>% subset_samples(Depth == "80") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "80"))>0) %>% tax_table() %>% row.names()

n9 <- full_ps_genus %>% subset_samples(Depth == "90") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "90"))>0) %>% tax_table() %>% row.names()

n10 <- full_ps_genus %>% subset_samples(Depth == "100") %>% 
  subset_taxa(taxa_sums(full_ps_genus %>% subset_samples(Depth == "100"))>0) %>% tax_table() %>% row.names()


n1_unique <- n1 %>% as.data.frame() %>% filter(. %notin% unique(c(n2,n3,n4,n5,n6,n7,n8,n9,n10)))
n2_unique <- n2 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n3,n4,n5,n6,n7,n8,n9,n10)))
n3_unique <- n3 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n4,n5,n6,n7,n8,n9,n10)))
n4_unique <- n4 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n5,n6,n7,n8,n9,n10)))
n5_unique <- n5 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n6,n7,n8,n9,n10)))
n6_unique <- n6 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n5,n7,n8,n9,n10)))
n7_unique <- n7 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n5,n6,n8,n9,n10)))
n8_unique <- n8 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n5,n6,n7,n9,n10)))
n9_unique <- n9 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n5,n6,n7,n8,n10)))
n10_unique <- n10 %>% as.data.frame() %>% filter(. %notin% unique(c(n1,n2,n3,n4,n5,n6,n7,n8,n9)))

temp <- as.data.frame(c(nrow(n1_unique),nrow(n2_unique),nrow(n3_unique),nrow(n4_unique),nrow(n5_unique),
                        nrow(n6_unique),nrow(n7_unique),nrow(n8_unique),nrow(n9_unique),nrow(n10_unique)),
                      row.names = c('10','20','30','40','50','60','70','80','90','100'))
colnames(temp) <- 'Number'
temp$Depth <- rownames(temp)
temp$Depth <- as.integer(temp$Depth)

unique_genera_depth <- ggplot(temp, aes(y = Number, fill = Depth, x = Depth)) + geom_bar(position = 'stack', stat = 'identity') + theme_bw()

ggsave('./Output/Figs/unique_genera_depth.pdf')
