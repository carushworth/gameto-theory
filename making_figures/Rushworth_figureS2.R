# Code for Rushworth et al Figure S2

library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
source("functions_forsummaries.R")
  library(gridExtra)
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  
##########################
##########################
##########################
  
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
  
  
  
  
  
  ramf_sum <- read_csv("finescaler_sum.csv")
  
  tmp <- ramf_sum %>% 
    mutate(z =  as.numeric(as.factor(r12)) ) %>% 
    filter(r12 %in% c(.5,.1,.01,.001,.0001,.00001,.000001))%>% 
    group_by(r12,z) %>%
    tally()%>% 
    ungroup()%>%
    mutate(breaks = z, cM = 100 * r12)
  
  
  scale_x_ramf <- scale_x_continuous(breaks = pull(tmp,breaks)[-length(pull(tmp,breaks))], 
                                     labels = c(expression(10^{-1}),expression(10^{-2}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                     expand = c(0,0), name = expression(italic(r[` AM        `])),
                                     lim = c(14,49))
  
  scale_y_ramf <- scale_y_continuous(breaks = pull(tmp,breaks)[-length(pull(tmp,breaks))], 
                                     labels = c(expression(10^{-1}),expression(10^{-2}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                     expand = c(0,0), name = expression(italic(r[` MF        `])),
                                     lim = c(14,49))
  
  
  
  

  ramf_Supp_amount  <-ramf_sum %>% 
    filter(s != 1, m0 == 0.01)%>%
    mutate(r12 = as.numeric(as.factor(r12)),
           r23 = as.numeric(as.factor(r23))) %>%
    ggplot(aes(x = r12, y = r23, fill = max.reinf_1))+
    geom_tile()+
    amount_scale+
    theme_light(base_family="Helvetica")+
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.height = unit(.3, "cm"),
          legend.margin=margin(t=7,b=3))+
    labs(fill = "Amount\n\n", y = "migration (m)", x = "selection (s)")+
    scale_x_ramf+
    scale_y_ramf+
    facet_wrap(~s,labeller = "label_both")
  
  ramf_Supp_dur  <- ramf_sum %>% 
    filter(s != 1, m0 == 0.01)%>%
    mutate(r12 = as.numeric(as.factor(r12)),
           r23 = as.numeric(as.factor(r23))) %>%
    ggplot(aes(x = r12, y = r23, fill  = reinf.time))+
    geom_tile()+
    time_scale+
    theme_light(base_family="Helvetica")+
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.height = unit(.3, "cm"),
          legend.margin=margin(t=7,b=3))+
    scale_x_ramf+
    scale_y_ramf+
    facet_wrap(~s,labeller = "label_both")+
    geom_tile(data = . %>%
                filter(!is.finite(reinf.time)),
              fill = "green4")+
    labs(fill = "Duration\n\n", y = "migration (m)", x = "selection (s)")+
    annotate(x = .08, y = 0.2, color = "white", label = "< 10\ngen", geom = "text", size = 2.5)
  
  
  
  
  
myplot <-  plot_grid(ramf_Supp_amount,ramf_Supp_dur,  labels = c("A","B"))
  ggsave(filename = "figS4.pdf",plot = myplot ,
         device = "pdf", width = 10.2, height =  3.37)
  
