library(phyloseq)
library(dplyr)
library(SpiecEasi)
library(vegan)
library(network)
library(igraph)
library(grid)
library(gridExtra)
library(reshape2)
library(plyr)

## Following 'Inference and Analysis of SPIEC-EASI Microbiome Networks' Birt & Dennis 2021
# Analysis of network data with permanova and visualization with Gephi
`%!in%` = Negate(`%in%`)
# Load data
ps.fung <- readRDS('./Output/noncontam_ps_object_fung.RDS') %>% subset_taxa(Class != 'NA') %>% tax_glom('Class')
ps.bac <- readRDS('./Output/noncontam_ps_object_bac.RDS') %>% subset_taxa(Class != 'NA') %>% tax_glom('Class')


# Ensure that the dataframes are in the same order (remove bacteria samples that were removed from fungi)
# Change tax_glom taxonomic resolution to desired
#ps.fung <- ps.fung %>% tax_glom("Genus")
#ps.bac <- ps.bac %>% tax_glom("Genus")

ps.fung.sub <- ps.fung %>%
  subset_samples(Depth == '10' | Depth == '20' | Depth == '30') %>%
  #subset_samples(Depth == '40' | Depth == '50' | Depth == '60') %>%
  #subset_samples(Depth == '70' | Depth == '80' | Depth == '90' | Depth == '100') %>%
  subset_samples(SampleID != 'UL7_30')

ps.fung.sub
ps.fung.sub@sam_data

ps.bac.sub <- ps.bac %>% # UL7_30 not in fungi data
  subset_samples(Depth == '10' | Depth == '20' | Depth == '30') %>%
  #subset_samples(Depth == '40' | Depth == '50' | Depth == '60') %>%
  #subset_samples(Depth == '70' | Depth == '80' | Depth == '90' | Depth == '100') %>%
  subset_samples(SampleID != 'UL7_30')

ps.fung.sub@sam_data$SampleID == ps.bac.sub@sam_data$SampleID #ensure that everything is correct

# Filtering ASVs to reduce ambiguous relationships and improve sparsity
thresh <- function(x){microbiome::prevalence(x) > 0.05 & colSums(microbiome::abundances(x)) > 5}

ps.bac.sub <- filter_taxa(ps.bac.sub, thresh, TRUE)

ps.fung.sub <- filter_taxa(ps.fung.sub, thresh, TRUE)

ps.bac.sub.f <- microbiomeutilities::format_to_besthit(ps.bac.sub)
ps.fung.sub.f <- microbiomeutilities::format_to_besthit(ps.fung.sub)

otu.bac.c <- t(otu_table(ps.bac.sub))
otu.fung.c <- t(otu_table(ps.fung.sub))

# Run SpiecEasi

# method = 'mb': meinshausen-buhlmann's neighborhood selection
# nlambda = 40 :
se.both <- spiec.easi(list(otu.fung.c, otu.bac.c), method='mb', nlambda=20,
                      lambda.min.ratio=1e-2, pulsar.params = list(rep.num = 50))

dtype <- c(rep(1,ntaxa(otu.fung.c)), rep(2,ntaxa(otu.bac.c)))
plot(adj2igraph(getRefit(se.both)), vertex.color=dtype+1, vertex.size=9)

saveRDS(se.both, './Output/surface_class_network.RDS')
beepr::beep()


#######################################
#Extract adjacency matrix. Indicates which pairs of ASVs are adjacent or not in the graph
se.both <- readRDS('./Output/surface_class_network.RDS')

getStability(se.both) #check stability of network - should be => 0.05

spieceasi.matrix <- symBeta(getOptBeta(se.both), mode = 'maxabs')
spieceasi.matrix <- as.matrix(spieceasi.matrix)
spieceasi.matrix <- readRDS('./Output/surface_spieceasi.matrix.RDS')

colnames(spieceasi.matrix) <- rownames(spieceasi.matrix) <- c(as.character(colnames(ps.bac.sub.f@otu_table)), as.character(colnames(ps.fung.sub.f@otu_table)))
asv.names <- colnames(spieceasi.matrix)

#Build weighted network. Edges in a weighted network represent the strength of association between ASVs
net <- graph.adjacency(spieceasi.matrix, mode = 'undirected', weighted = T, diag = F)
V(net)$name <- asv.names

#Convert edge weights into distances: larger weights = shorter distances -> then output a distance-based network
net.dist <- net
max(abs(E(net.dist)$weight))
weights.dist <- 1 - abs(E(net.dist)$weight)
E(net.dist)$weight <- weights.dist

#convert the weighted network to a separate absolute network
net.abs <- net
E(net.abs)$weight <- abs(E(net.abs)$weight)

##calculate centrality metrics and create a summary
#alpha centrality
net.alpha <- alpha.centrality(net)
#degree distribution
net.strength <- strength(net.abs)
#betweenness centrality
bet <- betweenness(net.dist, v = V(net.dist))

#summary of the metrics
summary_cent <- as.data.frame(net.alpha)
colnames(summary_cent) <- ("Alpha_centrality")
rownames(summary_cent) <- asv.names
summary_cent$Weighted_vertex_degree <- net.strength
summary_cent$Betweenness_centrality <- bet

metrics <- summary_cent

#cluster nodes into modules
wt <- cluster_louvain(net, weights = E(net.dist)$weight)
temp <- V(net)$name
temp <- as.data.frame(temp)
temp$louvain <- membership(wt)
V(net)$louvain <- temp$louvain

#Investigate notes that have been put into modules with <= 3 members and combine them
length(unique(temp$louvain))
summary_modules <- data.frame(table(temp$louvain))
colnames(summary_modules) <- c("louvain","n")
summary_modules
modules <- as.numeric(summary_modules$louvain[which(summary_modules$n>3)])

x <- max(modules) + 1
for (i in c(1:length(temp$temp))) {
  if (temp$louvain[i] %!in% modules) {
    temp$louvain[i] <- paste(x)
  }
}

modules <- temp
modules$louvain <- as.numeric(modules$louvain)
modules <- modules[order(modules$louvain),]
module.lookup <-
  data.frame("louvain" = unique(modules$louvain), "new_louvain" = c(1:length(unique(modules$louvain))))
new <- merge(modules, module.lookup)
modules <- new
modules <- modules[,2:3]
summary_modules <- data.frame(table(modules$new_louvain))
summary_modules
max(modules$new_louvain)

## Test whether centrality metrics of nodes differ between groups, when considered as a multivariate dataset
## Examine whether ASV's relative abundance affects its centrality metrics
# to include multiple metrics in the same model they must be z score transformed to all be on the same scale
# z score transformation
metrics.stand <- decostand(metrics, method = 'standardize')

# cannot be negative so we transform them all to make them positive
x <- abs(floor(min(metrics.stand)))
metrics.stand.abs <- metrics.stand + x
# Extract average abundance of each ASV
otu_fung <- as(otu_table(ps.fung.sub), 'matrix')
colnames(otu_fung) <- paste0('fung_', colnames(otu_fung)) 

otu_bact <- as(otu_table(ps.bac.sub), 'matrix')
otu.table <- merge(otu_bact, otu_fung, all= T)

av.abund <- as.data.frame(colMeans(otu.table))
colnames(av.abund) <- "average_abundance"
# Test whether abundance of an OTU significantly influences how important it is in the network
perm <- adonis2(metrics.stand.abs ~ av.abund$average_abundance)
perm # Abundance of OTU significantly influence how important it is in this network

# Evaluate differences in individual metrics between groups
cor.test(av.abund$average_abundance, metrics$Alpha_centrality, method = "pearson") # non-sig

cor.test(av.abund$average_abundance, metrics$Weighted_vertex_degree, method = "pearson") # sig

cor.test(av.abund$average_abundance, metrics$Betweenness_centrality, method = "pearson") # sig

# Test whether certain modules have higher centrality than other modules
modules.test <- as.data.frame(modules[which(modules$new_louvain != x),])
colnames(modules.test) <- c("OTU","louvain")
metrics.stand.abs.test <- metrics.stand.abs[which(modules$new_louvain  != x),]
metrics.test <- metrics[which(modules$new_louvain != x),]

#Test whether modules differ when considering all metrics in a single model
perm.1 <- adonis2(metrics.stand.abs.test ~ modules.test$louvain) #significant, the different metrics are significantly affected by OTU louvain
perm.1


# Prepare for gephi
spieceasi.matrix.m <- melt(spieceasi.matrix)

# name cols
colnames(spieceasi.matrix.m) <- c("source","target","weight")

# name of nodes
node.names <- unique(c(as.character(unique(spieceasi.matrix.m$source)), as.character(
  unique(spieceasi.matrix.m$target)
)))

# number them as an alphabetical node list, write to a csv
node.names <- as.data.frame(node.names)
node.names$node_number <- c(1:length(node.names$node.names))
node.names$node_number2 <- c(1:length(node.names$node.names))


colnames(node.names) <- c("Taxonomy","Label","Id")
row.names(node.names) <- node.names$Taxonomy
row.names(modules) <- modules$temp
modules <- modules[order(modules$temp),]
node.names <- node.names[order(node.names$Taxonomy),]
ps.arch.sub.f <- ps.bac.sub %>% subset_taxa(Kingdom == 'Archaea') %>% microbiomeutilities::format_to_besthit()
Kingdom <- ifelse(node.names$Taxonomy %in% colnames(ps.fung.sub.f@otu_table), "Fungi",
                  ifelse(node.names$Taxonomy %in% colnames(ps.arch.sub.f@otu_table), "Archaea", "Bacteria"))
node.names$Kingdom <- Kingdom
metrics$temp <- rownames(metrics)
metrics <- metrics[order(metrics$temp),]
row.names(node.names) == row.names(metrics)
row.names(node.names) == row.names(modules)
node.names.final <- cbind(node.names, metrics, modules)
node.names.final <- node.names.final[,-which(names(node.names.final) %in% c('temp'))]
write.table(node.names.final, "node.names.surface_class.csv", sep = ",", row.names = F)
write.table(node.names, "node.names.surface_class.csv", sep = ",", row.names = F)

# convert node names to numbers, write to a CSV
temp <- merge(x = spieceasi.matrix.m, y = node.names, by.x = "source", by.y = "Taxonomy")

# create the edge list
colnames(temp) <- c("source","target","weight","remove","source_number")
temp <- temp[,-4]
edge.list <- merge(x = temp, y = node.names, by.x = "target", by.y = "Taxonomy")
colnames(edge.list) <- c("source","target","weight","source.number","target.number")
edge.list <- edge.list[, c(3,4,6)]
colnames(edge.list) <- c("weight","source","target")
edge.list$Type <- "Undirected"
negative <- ifelse(edge.list$weight<0, "negative","positive")
edge.list$Negative <- negative
edge.list$weight <- abs(edge.list$weight)
edge.list <- edge.list[which(abs(edge.list$weight)>0),]
write.table(edge.list, "edge.list.surface_class.csv", sep = ",", row.names = F)

### edge list and node list can now be imported into Gephi ###

## Examining important nodes (ASVs)
node.names.final[order(node.names.final$Alpha_centrality, decreasing = T),][1,]
node.names.final[order(node.names.final$Betweenness_centrality, decreasing = T),][1,]


## Analysis of important fungi and bacteria
node.names.final <- read.table("./node.names_family_JD.csv", sep = ',', header = T)
node.names.impt <- node.names.final %>% filter(new_louvain != 10)
node.names.impt.fung <- node.names.impt %>% filter(Kingdom == 'Fungi') %>% .[order(.$Betweenness_centrality,decreasing = T),]
node.names.impt.bact <- node.names.impt %>% filter(Kingdom == 'Bacteria') %>% .[order(.$Betweenness_centrality,decreasing = T),]
  
library(ShortRead)
writeFasta(ps.fung@refseq['ASV1422'], './Output/refs/eurotiomycetes.fasta', mode = 'a')
writeFasta(ps.fung@refseq['ASV2547'], './Output/refs/herpotrichiellaceae.fasta', mode = 'a')
writeFasta(ps.bac@refseq['ASV1838'], './Output/refs/sandaracinaceae.fasta', mode = 'a')

writeFasta(ps.fung@refseq['ASV765'], './Output/refs/trichoderma.fasta', mode = 'a')
writeFasta(ps.bac@refseq['ASV220'], './Output/refs/napoli-4B-65.fasta', mode = 'a')
writeFasta(ps.fung@refseq['ASV920'], './Output/refs/hymenocallidicola.fasta', mode = 'a')

writeFasta(ps.bac@refseq['ASV362'], './Output/refs/bd2-11.fasta', mode = 'a')
writeFasta(ps.bac@refseq['ASV1695'], './Output/refs/anaerolineae.fasta', mode = 'a')
writeFasta(ps.bac@refseq['ASV3067'], './Output/refs/dehalococcoidia.fasta', mode = 'a')

