library(dplyr)
library(igraph)
library(SpiecEasi)
library(ggplot2)
library(emmeans)
source('./Zi_Pi_functions_2.R')

#### surface ####
edge.list_surface <- read.csv('./edge.list.surface_core.csv') %>% arrange(source)
node.names_surface <- read.csv('./node.names.surface_core.csv') %>% arrange(Id)

# create edgelist matrix
surface_el <- matrix(c(edge.list_surface$source, edge.list_surface$target), nc = 2)

# create igraph object - ensure that graph is undirected with directed = FALSE
surface_graph <- graph_from_edgelist(surface_el, directed = FALSE)
E(surface_graph)$weight <- edge.list_surface$weight # convert graph to weighted
surface_graph <- set.vertex.attribute(surface_graph, name = 'module', value = node.names_surface$new_louvain) # define modules
surface_graph <- set.vertex.attribute(surface_graph, name = 'name', value = node.names_surface$Id) # define names

# calculate zi and pi
node.names_surface$p <- part_coeff(surface_graph, get.vertex.attribute(surface_graph, 'module'), weighted = TRUE)
node.names_surface$z <- within_module_deg_z_score(surface_graph, get.vertex.attribute(surface_graph, 'module'), weighted = TRUE)

# plot with defined threshold for z > 2.5 as hub nodes and p > 0.62 as connector nodes (nodes with links to other modules)
surface_plot <- ggplot(data = node.names_surface, aes(x = p, y = z, colour = Kingdom)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('#619CFF','#F8766D','#00BA38')) +
  geom_hline(yintercept = 2.5, linetype = 'dashed') +
  geom_vline(xintercept = 0.62, linetype = 'dashed') +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab('Participation coefficient, P') +
  ylab('Within-module degree, z') +
  ylim(-3,5) +
  xlim(-0.01,1)

# identifying roles 
node.names_surface$roles <- assign_module_roles_4(node.names_surface %>% select(z,p))$roles

node.names_surface %>% group_by(roles, Kingdom) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

node.names_surface %>% group_by(Kingdom, roles) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(aes(x = Kingdom, y = z), data = node.names_surface) + 
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

ggplot(aes(x = Kingdom, y = p), data = node.names_surface) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

surface_p <- aov(p ~ Kingdom, data = node.names_surface)
summary(surface_p)
test(emmeans(surface_p, ~Kingdom))

#### subsurface ####
edge.list_subsurface <- read.csv('./edge.list.subsurface_core.csv') %>% arrange(source)
node.names_subsurface <- read.csv('./node.names.subsurface_core.csv') %>% arrange(Id)

# create edgelist matrix
subsurface_el <- matrix(c(edge.list_subsurface$source, edge.list_subsurface$target), nc = 2)

# create igraph object - ensure that graph is undirected with directed = FALSE
subsurface_graph <- graph_from_edgelist(subsurface_el, directed = FALSE)
E(subsurface_graph)$weight <- edge.list_subsurface$weight # convert graph to weighted
subsurface_graph <- set.vertex.attribute(subsurface_graph, name = 'module', value = node.names_subsurface$new_louvain) # define modules
subsurface_graph <- set.vertex.attribute(subsurface_graph, name = 'name', value = node.names_subsurface$Id) # define names

# calculate zi and pi
node.names_subsurface$p <- part_coeff(subsurface_graph, get.vertex.attribute(subsurface_graph, 'module'), weighted = TRUE)
node.names_subsurface$z <- within_module_deg_z_score(subsurface_graph, get.vertex.attribute(subsurface_graph, 'module'), weighted = TRUE)

# plot with defined threshold for z > 2.5 as hub nodes and p > 0.62 as connector nodes (nodes with links to other modules)
subsurface_plot <- ggplot(data = node.names_subsurface, aes(x = p, y = z, colour = Kingdom)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('#619CFF','#F8766D','#00BA38')) +
  geom_hline(yintercept = 2.5, linetype = 'dashed') +
  geom_vline(xintercept = 0.62, linetype = 'dashed') +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab('Participation coefficient, P') +
  ylab('Within-module degree, z') +
  ylim(-3,5) +
  xlim(-0.01,1)

# identifying roles 
node.names_subsurface$roles <- assign_module_roles_4(node.names_subsurface %>% select(z,p))$roles

node.names_subsurface %>% group_by(roles) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

node.names_subsurface %>% group_by(Kingdom, roles) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(aes(x = Kingdom, y = z), data = node.names_subsurface) + 
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

ggplot(aes(x = Kingdom, y = p), data = node.names_subsurface) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

subsurface_p <- aov(p ~ Kingdom, data = node.names_subsurface)
summary(subsurface_p)
test(emmeans(subsurface_p, ~Kingdom))

#### deep ####
edge.list_deep <- read.csv('./edge.list.deep_core.csv') %>% arrange(source)
node.names_deep <- read.csv('./node.names.deep_core.csv') %>% arrange(Id)

# create edgelist matrix
deep_el <- matrix(c(edge.list_deep$source, edge.list_deep$target), nc = 2)

# create igraph object - ensure that graph is undirected with directed = FALSE
deep_graph <- graph_from_edgelist(deep_el, directed = FALSE)
E(deep_graph)$weight <- edge.list_deep$weight # convert graph to weighted
deep_graph <- set.vertex.attribute(deep_graph, name = 'module', value = node.names_deep$new_louvain) # define modules
deep_graph <- set.vertex.attribute(deep_graph, name = 'name', value = node.names_deep$Id) # define names

# calculate zi and pi
node.names_deep$p <- part_coeff(deep_graph, get.vertex.attribute(deep_graph, 'module'), weighted = TRUE)
node.names_deep$z <- within_module_deg_z_score(deep_graph, get.vertex.attribute(deep_graph, 'module'), weighted = TRUE)

# plot with defined threshold for z > 2.5 as hub nodes and p > 0.62 as connector nodes (nodes with links to other modules)
deep_plot <- ggplot(data = node.names_deep, aes(x = p, y = z, colour = Kingdom)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('#619CFF','#F8766D','#00BA38')) +
  geom_hline(yintercept = 2.5, linetype = 'dashed') +
  geom_vline(xintercept = 0.62, linetype = 'dashed') +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab('Participation coefficient, P') +
  ylab('Within-module degree, z') +
  ylim(-3,5) +
  xlim(-0.01,1)

# identifying roles 
node.names_deep$roles <- assign_module_roles_4(node.names_deep %>% select(z,p))$roles

node.names_deep %>% group_by(roles, Kingdom) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

node.names_deep %>% group_by(Kingdom, roles) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n))

ggplot(aes(x = Kingdom, y = z), data = node.names_deep) + 
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

ggplot(aes(x = Kingdom, y = p), data = node.names_deep) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

deep_p <- aov(p ~ Kingdom, data = node.names_deep)
summary(deep_p)
test(emmeans(deep_p, ~Kingdom))


# collated data frame
node.names_combined <- data.frame('Names' = c(node.names_surface$Taxonomy, node.names_subsurface$Taxonomy, node.names_deep$Taxonomy),
                                  'Layer' = c(rep('Surface',nrow(node.names_surface)), rep('Subsurface', nrow(node.names_subsurface)), rep('Deep', nrow(node.names_deep))),
                                  'Kingdom' = c(node.names_surface$Kingdom, node.names_subsurface$Kingdom, node.names_deep$Kingdom),
                                  'Z' = c(node.names_surface$z, node.names_subsurface$z, node.names_deep$z),
                                  'P' = c(node.names_surface$p, node.names_subsurface$p, node.names_deep$p))

node.names_combined$Layer <- factor(node.names_combined$Layer, levels = c('Surface','Subsurface','Deep'))

ggplot(aes(x = Layer, y = P), data = node.names_combined) + 
  geom_jitter() +
  facet_grid(~Kingdom) +
  theme_bw()

ggplot(aes(x = Layer, y = Z), data = node.names_combined) + 
  geom_boxplot() +
  facet_grid(~Kingdom) +
  theme_bw()

comb_p <- vegan::adonis2(node.names_combined$P ~ Kingdom*Layer, data = node.names_combined, by = "terms", method = 'euclidean')
comb_z <- vegan::adonis2(node.names_combined$Z ~ Kingdom*Layer, data = node.names_combined, by = "terms", method = 'euclidean')

glm_p <- aov(node.names_combined$P ~ Kingdom*Layer, data = node.names_combined)
plot(glm_p)


library(pairwiseAdonis)
archaea_nodes <- node.names_combined %>% filter(Kingdom == 'Archaea')
archaea_nodes_pair <- pairwise.adonis(archaea_nodes$P, factors = archaea_nodes$Layer,
                sim.method = 'euclidean', p.adjust.m = 'fdr')

bacteria_nodes <- node.names_combined %>% filter(Kingdom == 'Bacteria')
bacteria_nodes_pair <- pairwise.adonis(bacteria_nodes$P, factors = bacteria_nodes$Kingdom,
                sim.method = 'euclidean', p.adjust.m = 'fdr')

fungi_nodes <- node.names_combined %>% filter(Kingdom == 'Fungi')
fungi_nodes_pair <- pairwise.adonis(fungi_nodes$P, factors = fungi_nodes$Kingdom,
                sim.method = 'euclidean', p.adjust.m = 'fdr')

fung_p <- kruskal.test(P ~ Layer, data = node.names_combined %>% filter(Kingdom == 'Fungi'))
pairwise.wilcox.test(node.names_combined$P, node.names_combined$Layer)
p_pairs <- emmeans::emmeans(comb_p, ~Kingdom:Layer)
test(p_pairs) %>% arrange(Kingdom)

ggplot(aes(x = Kingdom, y = P), data = node.names_combined) + 
  geom_boxplot() +
  facet_grid(~Layer) +
  theme_bw()

ggplot(aes(x = Kingdom, y = Z), data = node.names_combined) + 
  geom_boxplot() +
  facet_grid(~Layer) +
  theme_bw()
