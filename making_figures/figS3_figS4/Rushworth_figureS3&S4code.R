

library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
source("../gf_sims.R")



# The code below generate the data loaded for making figure S1
#M0.vals <- lapply( seq(0.05,.95,.1), function(M0) {
#  print(M0)
#  my.tmp <- runGFsim(n.gen = 5e+05, discrim = 1, s   =     .75,
#                     r12 =     .0001, r23 =     .0001, 
#                     init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0, fA_1 = 1, fM_1 = M0, fF_1 = 0.01),
#                     prop.replaced0 = .1,   prop.replaced1 = .1,  n.unlinked     = 0,
#                     delta_hap_components = FALSE,
#                     return.reinforce.only = TRUE,decompS = FALSE) 
#  return(as_tibble(my.tmp$geno.time) %>% mutate(M0 = M0))
#})

#write_csv(bind_rows(M0.vals), path = "M0.csv")


read_csv("M0.csv") %>%
  ggplot(aes(x = gen, y = reinf_1, group=M0,color = factor(M0)))+
  geom_line()+
  scale_x_continuous(trans = "log10", expand = c(0,0))+
  geom_hline(yintercept = 0, alpha = .4)+
  theme_tufte(base_family = "Helvetica")+
  theme(axis.line = element_line(color = "black"), legend.position = "bottom")+
  labs(x = "Generation", y = "Reinforcement", color = expression(M[0]), title = "Results do not depend on initial\nfrequency of compatibility allele")+
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 8, name = "RdPu"))(13)[-(1:2)])


# The code below generate the data loaded for making figure S2
#disc.vals <- lapply( seq(0,1,.1), function(DISC) {
#  print(DISC)
#  my.tmp <- runGFsim(n.gen = 5e+05, discrim = DISC, s   =     .75,
#                     r12 =     .0001, r23 =     .0001, 
#                     init.freqs =  c(fA_0 = 0, fM_0 = 0, fF_0 = 0, fA_1 = 1, fM_1 = 1, fF_1 = 0.01),
#                     prop.replaced0 = .1,   prop.replaced1 = .1,  n.unlinked     = 0,
#                     delta_hap_components = FALSE,
#                     return.reinforce.only = TRUE,decompS = FALSE) 
#  return(as_tibble(my.tmp$geno.time) %>% mutate(discrim = DISC))
#})
#write_csv(bind_rows(disc.vals), path = "discrim.csv")


read_csv("discrim.csv") %>%
  filter(discrim >0)%>%
  ggplot(aes(x = gen, y = reinf_1, group=discrim,color = factor(discrim)))+
  geom_line(show.legend = FALSE)+
  scale_x_continuous(trans = "log10", expand = c(0,0))+
  geom_hline(yintercept = 0, alpha = .4)+
  theme_tufte(base_family = "Helvetica")+
  theme(axis.line = element_line(color = "black"))+
  geom_text(data = .  %>% group_by(discrim) %>% filter(reinf_1 == max(reinf_1))%>%  ungroup()%>% mutate(gen = gen*(3)), 
            aes(label = discrim),show.legend = FALSE)+
  labs(x = "Generation", y = "Reinforcement", color = expression(M[0]), title = "Results quantitatively depend\non barrier strength")+
  scale_color_manual(values = colorRampPalette(brewer.pal(n = 8, name = "YlGn"))(15)[5:15])+
  annotate(x = .6e4, y = .8, label = "Barrier\nstrength (I)", geom = "text", hjust = 0)
  