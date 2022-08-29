library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
library(conflicted)
conflict_prefer("filter", "dplyr")
#source("/Users/catherinerushworth1/projects/gftheory_sandbox/functions_forsummaries.R")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[leg]
  return(legend)
}



### Supp Figure for two gfa





dur.color    <- colorRampPalette(brewer.pal(9, 'BrBG'))(1000)
amount_scale <- scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, 'YlGnBu'))(1000),
                                     limit = c(0,1))
time_scale   <- scale_fill_gradientn(colours = dur.color,
                                     trans = "log10",na.value = "tan4",
                                     limit = c(10,40000),
                                     breaks = c(10,100,1000,10000),
                                     labels = c(expression(10^{1}),expression(10^{2}),
                                                expression(10^{3}),expression(10^{4}))) 

time_scale2   <- scale_fill_gradientn(colours = dur.color,
                                      trans = "log10",na.value = "saddlebrown",
                                      limit = c(10,40000),
                                      oob = scales::squish_infinite(40000)) #?






### Load files
twogf_sum <- read_csv("~/Downloads/twogf_complete_summary.csv") %>%
  filter(m0 != 0)



twogf_amount <- ggplot(twogf_sum, aes(x = s, y = m0, fill = max.reinf_1)) +
  geom_tile()+
  amount_scale+
  theme_tufte(base_family="Helvetica") +
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5),
        legend.position = "bottom")+
  labs(x = "selection (s)", y = "migration (g)", fill = "amount")
  
twogf_dur <- ggplot(twogf_sum ,
                    aes(x = s, y = m0, fill = reinf.time))+
  geom_tile()+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = twogf_sum  %>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "#1d5a4c")+
  labs(x = "selection (s)", y = "migration (g)", fill = "duration")+
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5),
        legend.position =  "bottom")


plot_grid(twogf_amount, twogf_dur, labels = c('A', 'B'))
