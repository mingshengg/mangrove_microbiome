#### Network post analysis

library(emmeans)
library(ggplot2)
library(dplyr)

## centrality analysis - investigating if centrality metrics are significantly different across kingdoms ####
centrality <- readxl::read_xlsx('./combined_metrics.xlsx')
centrality$Layer <- factor(centrality$Layer, levels = c('surface','subsurface','deep'))

## eigenvalue
eigen_mod <- aov(eigen_centrality ~ Kingdom*Layer, data = centrality)
anova(eigen_mod)
plot(eigen_mod)
eigen_pairwise <- emmeans(eigen_mod, specs = pairwise~Kingdom|Layer, type = 'response')
eigen_pairwise$contrasts

write.csv(eigen_pairwise$contrasts, './Output/eigen_pairwise.csv')
write.csv(anova(eigen_mod), './Output/eigen_mod.csv')


ggplot(aes(x = Kingdom, y = eigen_centrality, fill = Kingdom), data = centrality) +
  geom_boxplot() + 
  theme_classic() + 
  facet_wrap(~Layer)

## degree
degree_mod <- aov(Degree ~ Kingdom*Layer, data = centrality)
summary(degree_mod)
plot(degree_mod)
degree_pairwise <- emmeans(degree_mod, specs = pairwise~Kingdom|Layer, type = 'response')
degree_pairwise$contrasts

write.csv(degree_pairwise$contrasts, './Output/degree_pairwise.csv')
write.csv(anova(degree_mod), './Output/degree_mod.csv')

ggplot(aes(x = Kingdom, y = Degree, fill = Kingdom), data = centrality) +
  geom_boxplot() + 
  theme_classic() + 
  facet_wrap(~Layer)

## closeness
closeness_mod <- aov((closness_centrality) ~ Kingdom*Layer, data = centrality)
summary(closeness_mod)
plot(closeness_mod)
closeness_pairwise <- emmeans(closeness_mod, specs = pairwise~Kingdom|Layer, type = 'response')
closeness_pairwise$contrasts

write.csv(closeness_pairwise$contrasts, './Output/closeness_pairwise.csv')
write.csv(anova(closeness_mod), './Output/closeness_mod.csv')


ggplot(aes(x = Kingdom, y = closness_centrality, fill = Kingdom), data = centrality) +
  geom_boxplot() + 
  theme_classic() + 
  facet_wrap(~Layer)

## betweenness
betweenness_mod <- aov((betweenness_centrality) ~ Kingdom*Layer, data = centrality)
summary(betweenness_mod)
plot(betweenness_mod)
betweenness_pairwise <- emmeans(betweenness_mod, specs = pairwise~Kingdom|Layer, type = 'response')
betweenness_pairwise$contrasts

write.csv(betweenness_pairwise$contrasts, './Output/betweenness_pairwise.csv')
write.csv(anova(betweenness_mod), './Output/betweenness_mod.csv')


ggplot(aes(x = Kingdom, y = (betweenness_centrality), fill = Kingdom), data = centrality) +
  geom_boxplot() + 
  theme_classic() + 
  facet_wrap(~Layer)

# positive vs negative edges ####
surface_edgelist<-read.csv('./edge.list.surface_core.csv')

surface_edgelist$source <- as.factor(surface_edgelist$source)

unique_node <- surface_edgelist %>% count(source) %>% select(source)
surface_edge_number <- data.frame(id = unique_node, positive = rep("0", ncol(unique_node)), negative = rep("0", ncol(unique_node)))

# for each node, calculate number of positive and negative edges
for (i in surface_edgelist$source){
  temp <- surface_edgelist %>% filter(source == i) %>% count(Negative)
  surface_edge_number[surface_edge_number$source == i,]$positive <- temp %>% filter(Negative == 'positive') %>% .$n
  surface_edge_number[surface_edge_number$source == i,]$negative <- ifelse(match('negative', temp$Negative), temp %>% filter(Negative == 'negative') %>% .$n, 0)
}

surface_edge_number['negative'][is.na(surface_edge_number['negative'])] <- 0

node.names <- read.csv('./node.names.surface_core.csv')
node.names <- node.names[order(node.names$Label),]

all(surface_edge_number$source == node.names$Label)

surface_edge_number$kingdom <- node.names$Kingdom

prop_edge <- surface_edge_number
surface_edge_number$positive <- as.numeric(surface_edge_number$positive)
surface_edge_number$negative <- as.numeric(surface_edge_number$negative)
prop_edge$positive <- surface_edge_number$positive/(surface_edge_number$positive + surface_edge_number$negative)
prop_edge$negative <- surface_edge_number$negative/(surface_edge_number$positive + surface_edge_number$negative)

library(ggplot2)

ggplot(data = prop_edge, aes(x = kingdom, y = positive)) +
  geom_boxplot()

temp <- aov(positive ~ kingdom, data = surface_edge_number)
summary(temp)

par <- emmeans::emmeans(temp, ~kingdom)
emmeans::test(par)
aov(par)

ggplot(data = surface_edge_number, aes(x = kingdom, y = positive)) +
  geom_boxplot()
