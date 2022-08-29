# Code for Rushworth et al Figure S1

library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
source("../scripts/functions_forsummaries.R")
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
  


  
moreasymS2_sum <- read_csv("moreasyms_sum.csv")    
  
  supp_asymS_amount <- moreasymS2_sum %>%
    filter(m0 !=0 & m0 !=.3)%>%
    rename(m = m0)%>%
    select(-m1)%>%
    ggplot(aes(x = s0, y = s1, fill = max.reinf_1))+
    geom_tile()+
    amount_scale+
    theme_light(base_family="Helvetica")+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    labs(x =  bquote(s [~maize]), y =  bquote(s[~teo]), fill = "Amount\n\n")+
    facet_wrap(~m, labeller = "label_both")+
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.height = unit(.3, "cm"),
          legend.margin=margin(t=7,b=3))
  
  
  
  
  
  supp_asymS_dur <- moreasymS2_sum %>%
    filter(m0 !=0 & m0 !=.3)%>%
    rename(m = m0)%>%
    select(-m1)%>%
    ggplot(aes(x = s0, y = s1, fill = reinf.time))+
    geom_tile()+
    time_scale+
    theme_light(base_family="Helvetica")+
    theme(legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.height = unit(.3, "cm"),
          legend.margin=margin(t=7,b=3))+
    scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    labs(x =  bquote(s [~maize]), y =  bquote(s[~teo]),fill = "Duration\n\n")+
    geom_tile(data = . %>%
                filter(!is.finite(reinf.time)),
              fill = "green4")+
    facet_wrap(~m,labeller = "label_both")+
    geom_tile(data = . %>%
                filter(!is.finite(reinf.time)),
              fill = "green4")
  
  
  plot_grid(supp_asymS_amount,  supp_asymS_dur,labels = c("A","B"))
