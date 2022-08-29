# Code for Rushworth et al Figure 5

library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
#######################
###### Figure 5. ######
#######################
  
dur.color    <- colorRampPalette(brewer.pal(9, 'BrBG'))(1000)
amount_scale <- scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, 'YlGnBu'))(1000),
                                     limit = c(0,1))
time_scale   <- scale_fill_gradientn(colours = dur.color,
                                     trans = "log10",na.value = "tan4",
                                     limit = c(10,40000),
                                     breaks = c(10,100,1000,10000),
                                     labels = c(expression(10^{1}),expression(10^{2}),
                                                expression(10^{3}),expression(10^{4}))) 

unlinked_amount <- bind_rows(read_csv("nunlinked1_sum_fix.csv"),
                             read_csv("nunlinked2_sum_fix.csv"),
                             read_csv("nunlinked3_sum_fix.csv"),
                             read_csv("nunlinked4_sum_fix.csv"),
                             read_csv("nunlinked5_sum_fix.csv"))%>% 
  mutate("n loci" = number_loci+1)%>%
  filter(m0>0)%>%
  ggplot(aes(x = s, y = m0, fill = max.reinf_1))+
  geom_tile()+
  scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, 'YlGnBu'))(1000),
                       limit = c(0,1), breaks  = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
  labs()+
  facet_wrap(~`n loci`, ncol = 3, label = "label_both") +
  theme_tufte(base_family = "Helvetica")+
  theme(legend.position = c(0.85, 0.2), legend.key.width  = unit(.5,"cm"),legend.direction = "horizontal",
        strip.text = element_text(color = "black"),strip.background = element_rect(fill = "lightgrey",colour = "white")
  )+
  labs(fill = "Amount", y = "migration (m)", x = "selection (s)")+
  guides(fill= guide_colourbar(title.position="top", title.hjust = 0.5))+
  scale_x_continuous(breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))


unlinked_dur <- bind_rows(read_csv("nunlinked1_sum_fix.csv"),
                          read_csv("nunlinked2_sum_fix.csv"),
                          read_csv("nunlinked3_sum_fix.csv"),
                          read_csv("nunlinked4_sum_fix.csv"),
                          read_csv("nunlinked5_sum_fix.csv"))%>% 
  filter(m0>0)%>%
  mutate("n loci" = number_loci+1)%>%
  ggplot(aes(x = s, y = m0, fill = reinf.time))+
  geom_tile()+
  time_scale+
  labs()+
  facet_wrap(~`n loci`, ncol = 3, label = "label_both") +
  theme_tufte(base_family = "Helvetica")+
  theme(legend.position = c(0.85, 0.2), legend.key.width  = unit(.5,"cm"),legend.direction = "horizontal",
        strip.text = element_text(color = "black"),strip.background = element_rect(fill = "lightgrey",colour = "white")
  )+
  labs(fill = "Duration", y = "migration (g)", x = "selection (s)")+
  guides(fill= guide_colourbar(title.position="top", title.hjust = 0.5))+
  #  guides(fill= guide_colourbar(title.position="left", title.hjust = 0.5,
  #                               title.theme = element_text(angle = 90)))+
  scale_x_continuous(breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  geom_tile(data = . %>%
              filter(!is.finite(reinf.time)),
            fill = "green4")

plot_grid(unlinked_amount,unlinked_dur, ncol = 1,  labels = c("A","B"))








  