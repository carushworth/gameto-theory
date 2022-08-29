# Rushworth et al Figures 2 & 3
library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)


#######################
###### Figure 2. ######
#######################

# load and process data
load("exampleDynamics.Robj") # loads "output"
output <- output[!names(output) == "meanUs"]        # remove "meanUs" which only has meaning in th unlinked model
output <- output[!names(output) == "summary.stats"] # we did not use these summary stats

# Output is a list with elements: "geno.time", "dhaps", "marj", "params"   

   # geno.time: is matrix of genotype frequencies across populations over time. 
     # each row is a generation 
     # columns starting with X are genotypes in populations 
     # where the first three numbers reflect maternal haplotypes 
     # and the next three numbers reflect paternal haplotypes 
     # at the A, M, and F loci in that order
     # where 0 refers to a, m, or f, and 1 refers to A, M, or F.
     # The final number is population, where 0 is maize and 1 is teosinte

   # dhaps: is matrix of change in haplotype frequencies across the life cycle every generation. 
     # Each row is a generation 
     # Haplotypes show the A, M, and F loci in that order
     # where 0 refers to a, m, or f, and 1 refers to A, M, or F.
     # p0 and p1 refer to maize and teosinte populations, respectively.

   # marj: is a matrix summarizing linked and direct selection on the female incompatability 
     # allele each generation in each population, 
     # where 0 and one refer to maize and teosinte, respectively

   # params: is a vector containing parameter values used in this simulation

# thin data so as to not overwhelm R or the graphs
to.keep  <- which(output$geno.time$gen %% floor(nrow(output$geno.time) / 4000) ==1 |  c(1, rowSums(abs(output$geno.time[-1,-c(1:3)] - output$geno.time[-nrow(output$geno.time),-c(1:3)]  ))) > 5e-4)
output$geno.time <- slice(output$geno.time, to.keep)
output$meanUs    <- slice(output$meanUs, to.keep)
output$dhaps     <- slice(output$dhaps, to.keep)
output$marj      <- output$marj %>%
  mutate(gen = 1:nrow(output$marj)) %>%
  slice(to.keep)


## Standardize reinforcement by value at first generation
output$geno.time <- output$geno.time %>%      #
  mutate(teo_reinforce = (reinf_1 - reinf_1[1])  / (1- reinf_1[1])  )




#### tidying data for plotting

source("functions_forsummaries.R")

# tidying data
tidy.haps           <- tidyingHaps(output$geno.time) 
allele.freqs.and.ld <- findFreqs(tidy.haps)
tidy.allele.freqs   <- tidyingAlleleFreqs(allele.freqs.and.ld)
tidy.ld             <- tidyingLD(allele.freqs.and.ld)
tidy.freq.change    <- tidyingFreqChange(output$dhaps)
tidy.allele.change <- alleleChange(tidy.freq.change)
tidy.haps           <- tidy.haps  %>%  
  filter(!duplicated(paste(  pop,   haplo,round(log10(gen) ,digits = 2)   ,sep="")))%>%
  mutate(pop = factor(ifelse(pop == 1,"teosinte","maize"),levels = c("teosinte","maize")))
tidy.ld.r2             <- allele.freqs.and.ld %>%
  mutate(ld_MF = ld_MF / sqrt(freq_M * (1 -freq_M ) * freq_F * (1 - freq_F)),
         ld_AM = ld_AM / sqrt( freq_A* (1 - freq_A) * freq_M * (1 - freq_M)),
         ld_AF = ld_AF / sqrt( freq_A* (1 - freq_A) * freq_F * (1 - freq_F)))%>%
  tidyingLD()

tidy.ld <-  tidy.ld %>%
  filter(!duplicated(paste(  pop,   pair,round(log10(gen) ,digits = 2)   ,sep="")))%>%
  mutate(pop = factor(ifelse(pop == 1,"teosinte","maize"),levels = c("teosinte","maize")))
tidy.ld.r2 <-  tidy.ld.r2 %>%
  filter(!duplicated(paste(  pop,   pair,round(log10(gen) ,digits = 2)   ,sep="")))%>%
  mutate(pop = factor(ifelse(pop == 1,"teosinte","maize"),levels = c("teosinte","maize")))%>%
  filter(!is.na(ld), is.finite(ld), abs(ld) <1.000001) # remove numerically unstable results


# Defining 'phases'
phase1 <- 1
phase2 <- output$geno.time %>% 
  filter(teo_reinforce==max(teo_reinforce)) %>% 
  pull(gen)
phase3 <- tidy.ld %>% 
  filter(pop == "teosinte", pair == "AF")%>%
  filter(ld ==  max(ld))%>%
  pull(gen)

#annotation layer for plot
z <- annotate(geom = "rect",
              xmin = c(phase1,phase3), xmax = c(phase2,Inf), 
              ymin =c(-Inf,-Inf), ymax = c(Inf,Inf),
              fill =  c("lightblue","lightgrey"), alpha = c(.15,.2))
gen_scale <- scale_x_continuous(trans = "log10", breaks =(10^(0:6)), 
                                labels = c(expression(10^0),expression(10^1),
                                           expression(10^2),expression(10^3),
                                           expression(10^4),expression(10^5),
                                           expression(10^6)),
                                limits = c(1,5e5),name = "Generation")
#color schemes
allele.colors <- wes_palette("FantasticFox1")[c(2,3,1)]

# 2A: reinforcement_dynamics 
reinforcement_dynamics  <- output$geno.time%>%
  ggplot(aes(x = gen, y = teo_reinforce))+
  geom_line()+
  gen_scale  +
  theme_tufte(base_family = "Helvetica")+
  ylab("Reinforcement     ")+
  # ylab("Reinforcement     ")+
  scale_y_continuous(limits = c(-.1,1),expand = c(0,0))+
  geom_hline(yintercept = 0,color = "red",size = 1.5, alpha = .2)+
  annotate(geom = "text",x = c(1,190, 1000), y = c(.95,.95,.8), 
           label = c("Phase 1", "Phase 2","Phase 3"),
           hjust = 0, size = 3)+
  theme(axis.line = element_line(color = "black"))+
  z

# 2B: allele_freq_dynamics
allele_freq_dynamics <- tidy.allele.freqs  %>%  
  filter(!duplicated(paste(  pop,   allele,round(log10(gen) ,digits = 2)   ,sep="")))%>%
  mutate(pop = factor(ifelse(pop == 1,"teosinte","maize"),levels = c("teosinte","maize")))%>%
  ggplot(aes(x = gen, y  = freq, color = allele, lty = pop))+
  geom_line(show.legend = FALSE)              +
  gen_scale+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  theme_tufte(base_family = "Helvetica")+
  ylab("Allele freq")+ 
  scale_color_manual(values = allele.colors)+
  annotate(x = c(1,7, 500, 1.5e3,40000,70000), 
           y = c(.9,.6,.6, .1,.9,.1), 
           geom="text",
           color = c(allele.colors[3],allele.colors[2],allele.colors[3],
                     allele.colors[2],"gold","gold"),
           label = c(expression(M[teosinte]),expression(F[teosinte]) ,expression(M[maize]) ,
                     expression(F[maize]),expression(A[teosinte]) ,expression(A[maize])),
           hjust = 0, size = 3  )+
  theme(axis.line = element_line(color = "black"))+
  z

# 2C: teo_hap_dynamics
teo_hap_dynamics <- tidy.haps %>%
  filter(pop != "maize")%>%
  ggplot(aes(x = gen, y = freq, color = pollen_style, lty = factor(pop), group = paste(pop,pollen_style,local_adapt)))     +
  geom_line(show.legend = FALSE)              +
  gen_scale+
  scale_y_continuous(expand = c(0,0),limits = c(0,1))+
  theme_tufte(base_family = "Helvetica")+
  ylab("Hap freq: teosinte        ")+ 
  scale_color_manual(values = c("#78c679","#ffffcc","#41b6c4", "#225ea8"))+
  annotate(x = c(10,2000,2,2,10,10,2,10,4), 
           y = c(.9,.9,.35,.225,.35,.225,.1,.1,.5), 
           geom="text",color = c("#41b6c4", "#225ea8","#78c679","gold","#78c679","gold","#41b6c4", "#225ea8","grey"), 
           label = c("AMf","AMF","Amf","AmF","amf","amF","aMf","aMF","Rare haps"),
           size = c(3,3,2,2,2,2,2,2,2))+
  theme(axis.line = element_line(color = "black"))+
  z




# 2D: teo_ld_dynamics
teo_ld_dynamics <- tidy.ld.r2  %>%
  filter(pop != "maize")%>%
  ggplot(aes(x = gen, y = ld, color = pair, group = paste(pair,pop),lty = pop))     +
  geom_line(show.legend = FALSE) +   
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  annotate(x = c(1, 2e5, 1), 
           y = c(0.9, .085, .085), 
           geom="text",
           color = (wes_palette("Darjeeling1")[c(2,1,3)]),
           label = c("AM","AF","MF"),
           hjust = 0, size = 3  )+
  ylab(expression(paste("LD: teosinte (",r^2,")")))+
  theme(axis.line = element_line(color = "black"))  +
  z


# 2E: maize_hap_dynamics
maize_hap_dynamics <- tidy.haps%>%
  filter(pop == "maize")%>%
  ggplot(aes(x = gen, y = freq, color = pollen_style, group = factor(paste(pollen_style,local_adapt))))     +
  geom_line(show.legend = FALSE)+
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  ylab("Hap freq: maize    ")+ 
  scale_color_manual(values = c("#78c679","#ffffcc","#41b6c4", "#225ea8"))+
  annotate(x = c(10,2000,10000),
           y = c(.9,.9,.1), 
           geom="text",
           color = c("#78c679","#41b6c4", "#225ea8"), 
           label = c("amf","aMf","aMF"), size = 3)+
  annotate(x = c(2,2,10,2,10,4), 
           y = c(.35,.225,.225,.1,.1,.5), 
           geom="text",color = c("#78c679","gold","gold","#41b6c4", "#225ea8","grey"), 
           label = c("Amf","AmF","amF","AMf","AMF","Rare haps"),
           size = 2)+
  theme(axis.line = element_line(color = "black"))+
  z

# 2F: maize_ld_dynamics
maize_ld_dynamics  <- tidy.ld.r2  %>%  
  filter(pop == "maize")%>%
  ggplot(aes(x = gen, y = ld, color = pair))     +
  geom_line(lty = 1, show.legend = FALSE) +               
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  annotate(x = c(1, 2e5, 1), 
           y = c(.9, .085, .15), 
           geom="text",
           color = (wes_palette("Darjeeling1")[c(2,1,3)]),
           label = c("AM","AF","MF"),
           hjust = 0, size = 3  )+
  ylab("LD: maize    ")+
  ylab(expression(paste("LD: maize   (",r^2,")")))+
  theme(axis.line = element_line(color = "black"))+
  z

plot_grid(reinforcement_dynamics,
          allele_freq_dynamics,
          teo_hap_dynamics,
          teo_ld_dynamics,
          maize_hap_dynamics,
          maize_ld_dynamics,
          ncol=2,
          labels = c("A","B","C","D","E","F"))






#######################
###### Figure 3. ######
#######################




tidy.allele.change <- tidy.allele.change  %>%
  filter(gen %in% unique(output$geno.time$gen))%>%
  mutate(pop = factor(ifelse(pop == 1,"teosinte","maize"),levels = c("teosinte","maize")))

allele_freq_change_scale <- scale_y_continuous(limits = c(-0.05,0.05),
                                               breaks = seq(-0.05,0.05,0.025),
                                               labels = c("-.050","-.025",".000",".025",".050"))

###### Fig 3
allele_change_mig <- tidy.allele.change %>%  
  rename(Subspecies= "pop", Allele = "allele")%>%
  filter(time == "migrationPollen")%>%
  ggplot(aes(x = gen, y  = delta_p, color = Allele, lty = Subspecies))+
  geom_line()              +
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  ylab(expression(paste(Delta~p[migration],"     " )))+ 
  geom_hline(yintercept = 0,color = "red",size = 1.5, alpha = .2)+
  scale_color_manual(values = allele.colors)+
  theme(axis.line = element_line(color = "black"))+
  z+
  annotate(x = c(650,650, 3*10^5, 3*10^5, 3*10^5,  3*10^5),
           y = c(-0.02,.02,  .02,     -0.02, .04  ,  -.04),
           geom = "text",
           size = 3,
           color = c( allele.colors[3], allele.colors[3], allele.colors[2], allele.colors[2],"gold","gold"),
           label = c("M","M","F","F","A","A"))+
  allele_freq_change_scale 



allele_change_mating <- tidy.allele.change %>%  
  filter(time == "matingPollen")%>%
  ggplot(aes(x = gen, y  = delta_p, color = allele, lty = pop))+
  geom_line(show.legend = FALSE)              +
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  ylab(expression(paste(Delta~p[fertilization],"     " )))+ 
  geom_hline(yintercept = 0,color = "red",size = 1.5, alpha = .2)+
  scale_color_manual(values =  allele.colors)+
  theme(axis.line = element_line(color = "black"))+
  z+
  allele_freq_change_scale +
  annotate(x = c(550,700, 600),
           y = c(0.05, 0.03, 0.04),
           geom = "text",
           size = 3,
           color = c( allele.colors[3], allele.colors[2],"gold"),
           label = c("M","F","A"))


allele_change_selection <- tidy.allele.change %>%  
  filter(time == "sel")%>%
  ggplot(aes(x = gen, y  = delta_p, color = allele, lty = pop))+
  geom_line(show.legend = FALSE)              +
  gen_scale+
  theme_tufte(base_family = "Helvetica")+
  ylab(expression(paste(Delta~p[ selection],"     " )))+ 
  geom_hline(yintercept = 0,color = "red",size = 1.5, alpha = .2)+
  scale_color_manual(values = allele.colors)+
  theme(axis.line = element_line(color = "black"))+
  z+  
  allele_freq_change_scale +
  annotate(x = c(10,        10, 3*10^5, 3*10^5, 3*10^5,  3*10^5),
           y = c(-0.0325,.0325,  .015,     -0.015, .0325,  -.0325),
           geom = "text",
           size = 3,
           color = c(allele.colors[3],allele.colors[3],allele.colors[2],allele.colors[2],"gold","gold"),
           label = c("M","M","F","F","A","A"))


marj_wdat <- pivot_longer(output$marj, cols = -1, values_to = "marginal s", 
                          names_sep = "_", names_to = c("allele","pop")) %>%
  mutate(Population = factor(case_when(pop == 0 ~ paste("- - - ","  maize", sep=""),pop == 1~paste("\u2013\u2013\u2013","  teosinte", sep=""))),
         Population = fct_rev(Population))%>%
  filter(gen %in% unique(output$geno.time$gen))


marj_wdat_fig <- ggplot(marj_wdat %>% filter(!is.na(`marginal s`)),
                        aes(x = gen, y = `marginal s`, lty = Population, group = paste(allele, Population),  color = allele))+
  geom_line(show.legend = FALSE)+
  facet_wrap(~Population, ncol = 1, scale = "free_y")+
  gen_scale+  
  #scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  theme_tufte(base_family = "Helvetica")+
  scale_color_manual(values = c("black","grey"))+
  geom_text(data = .%>%
              group_by(Population,pop,allele)%>% 
              tally()%>% 
              ungroup()%>%
              mutate(label = c("On F",
                               "Linked",
                               "On F",
                               "Linked"),
                     gen = c(12, 80000, 12,80000), 
                     'marginal s' = c(.055,.075,-.3,-.4)),
            aes(label = label), size = 3, show.legend = FALSE, color = "black")+
  scale_y_continuous(n.breaks = 3)+
  theme(axis.line = element_line(color = "black"))+
  annotate(geom = "rect",xmin  = phase1, xmax = phase2,  
           ymin  = c(-Inf), ymax = c(Inf),fill  = c("lightblue"), alpha = c(.15))+
  annotate(geom = "rect",xmin  = phase3, xmax = c(Inf),  
           ymin  = c(-Inf), ymax = c(Inf),fill  = c("lightgrey"), alpha = c(.15))


top <-  plot_grid(allele_change_mig + theme(legend.position = "none"),
                  allele_change_mating,
                  allele_change_selection,
                  ncol=3,
                  labels = c("A","B","C"))
left <- plot_grid(top, 
                  get_legend(allele_change_mig+theme(legend.position = "bottom")),
                  ncol = 1,
                  rel_heights = c(1,.2))

plot_grid(left,marj_wdat_fig, ncol = 2,rel_widths = c(3,1), labels = c("","D"))
