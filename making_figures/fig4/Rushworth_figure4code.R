# Code for Rushworth et al Figure 4

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
###### Figure 4. ######
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
  
  time_scale2   <- scale_fill_gradientn(colours = dur.color,
                                       trans = "log10",na.value = "saddlebrown",
                                       limit = c(10,40000),
                                       oob = scales::squish_infinite(40000)) #?


  
### FIgure 4A&E symmetric migration and selection 
  
smint_sum  <- read_csv("smint_sum.csv")  
  
  mig_sel_amount <- ggplot(smint_sum, aes(x = s, y = m0, fill = max.reinf_1))+
    geom_tile()+
    amount_scale+
    theme_tufte(base_family="Helvetica")+ 
    theme(legend.position = "top", legend.key.width = unit(.55, "cm"), legend.key.height = unit(.3, "cm"),
          legend.margin=margin(t=7,b=3))+
    scale_x_continuous(expand = c(0,0),  limits = c(0,1),breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,0.35))+
    labs(fill = "Amount\n\n", y = "migration (g)       ", x = "selection (s)")

  
  mig_sel_dur <- ggplot(smint_sum, aes(x = s, y = m0, fill = reinf.time ))+
    geom_tile()+
    time_scale+
    theme_tufte(base_family="Helvetica")+
    scale_x_continuous(expand = c(0,0), limits = c(0,1),breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,0.35))+    #  theme(legend.key.height = unit(.5, "cm"),
    labs(fill = "Duration\n\n", y = "migration (g)       ", x = "selection (s)")+
    annotate(x = .08, y = 0.2, color = "white", label = "< 10\ngen", geom = "text", size = 2.5)
  
  

  
  
  
  
  
  
  
### FIgure 4B&F asymmetric selection  
  
asyms_sum <- read_csv("asyms_sum.csv")   

asymS_amount <- asyms_sum %>%
  filter(m0==0.03)%>%
  ggplot(aes(x = s0, y = s1, fill = max.reinf_1))+
  geom_tile(show.legend = FALSE)+
  amount_scale+
  theme_tufte(base_family="Helvetica")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
  labs(x =  bquote(s [~maize]), y =  bquote(s[~teo]))

  
asymS_dur <- asyms_sum %>%
  filter(m0==0.03)%>%
  ggplot(aes(x = s0, y = s1, fill = reinf.time))+
  geom_tile(show.legend = FALSE)+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,1,.25), labels = c("0","0.25","0.5","0.75","1"))+
  labs(x =  bquote(s [~maize]), y =  bquote(s[~teo]))+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = asyms_sum %>%
              filter(m0==0.03)%>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "green4")+
  annotate(x = .035,y=.07,color = "white", label  = "<10",geom ="text", size= 2)+
  annotate(x = .991,y=.9, color = "white", label  =  expression(infinity) ,geom ="text", size= 2.5 )
  
  
  
  
  
  
### Figure 4C&G asymmetric migration

asymm_sum  <- read_csv("asymm_sum.csv") 

asymM_amount <- asymm_sum %>% filter(s != 0 & s!=1)%>%
  ggplot(aes(x = m1, y = m0))+
  geom_tile(aes(fill = max.reinf_1),show.legend = FALSE)+
  facet_wrap(~s, labeller = "label_both")+
  amount_scale+
  scale_x_continuous(breaks = c(0,0.2,.4),labels = c("0",".2",".4"))+
  scale_y_continuous(breaks = c(0,0.2,.4))+ 
  theme_tufte()+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  geom_text(data = . %>% group_by(s) %>% tally() %>%mutate(m1 = .25,m0 = .45, sel = paste("s",s,sep = " = "), col = s <0.7) , 
            aes(label = sel,color = col), size = 2.2,show.legend = FALSE)+
  scale_color_manual(values = c("white","black"))+
  labs(x =  bquote(g [~maize %->% teo]), y =  bquote(g[~teo %->% maize]~"     "))

asymM_dur <-asymm_sum %>% filter(s != 0 & s!=1)%>%
  ggplot(aes(x = m1, y = m0))+
  geom_tile(aes(fill = reinf.time),show.legend = FALSE)+
  facet_wrap(~s, labeller = "label_both")+
  time_scale+
  scale_x_continuous(breaks = c(0,0.2,.4),labels = c("0",".2",".4"))+scale_y_continuous(breaks = c(0,0.2,.4))+ 
  theme_tufte(base_family="Helvetica") +
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  geom_text(data = . %>% group_by(s) %>% tally() %>%mutate(m1 = .25,m0 = .45, sel = paste("s",s,sep = " = "), col = s >0.3) , 
            aes(label = sel,color =col), size = 2.2,show.legend = FALSE)+
  scale_color_manual(values = c("white","black"))+
  labs(x =  bquote(g [~maize %->% teo]), y =  bquote(g[~teo %->% maize]~"     "))




  
  
### Figure 4D&H recombination    

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
#scale_x_ramf <- scale_x_continuous(breaks = pull(tmp,breaks)[-length(pull(tmp,breaks))], 
#                                labels = c(expression(10^{-1}),expression(10^{-2}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
#                                expand = c(0,0), name = expression(italic(r[` AM        `])),
#                                lim = c(14,49))
#scale_y_ramf <- scale_y_continuous(breaks = pull(tmp,breaks)[-length(pull(tmp,breaks))], 
#                                labels = rev(c(expression(10^{-6}),expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1}))), 
#                                expand = c(0,0),name = expression(italic(r[` MF        `])),
#                                lim = c(14,49))

ramf_plot_format_sub <-  ramf_sum %>% 
  filter(s == 0.8, m0 == 0.01)%>%
  mutate(r12 = as.numeric(as.factor(r12)),
         r23 = as.numeric(as.factor(r23)))



ramf_amount <- ggplot(ramf_plot_format_sub ,
                      aes(x = r12, y = r23, fill = max.reinf_1))+
  geom_tile(show.legend = FALSE)+
  amount_scale+
  theme_tufte(base_family="Helvetica")+
  scale_x_ramf+
  scale_y_ramf

ramf_dur <- ggplot(ramf_plot_format_sub ,
                   aes(x = r12, y = r23,fill = reinf.time))+
  geom_tile(show.legend = FALSE)+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  scale_x_ramf+
  scale_y_ramf+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = ramf_plot_format_sub  %>%
              filter(m0==0.01)%>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "green4")





####  Figure 4  

plot_grid(mig_sel_amount+ylim(c(0,.3)) + theme(legend.position="none"), 
          asymS_amount,
          asymM_amount,
          ramf_amount,
          get_legend(mig_sel_amount +labs(fill = "Amount")+ 
                       theme(legend.position= "right",
                             legend.justification = "top", 
                             legend.key.width = unit(.45, "cm"),
                             legend.key.height = unit(.45, "cm"),
                             legend.margin =  margin(t=0,b=2))),
          mig_sel_dur +ylim(c(0,.3))+ theme(legend.position="none"),
          asymS_dur,
          asymM_dur,
          ramf_dur,
          get_legend(mig_sel_dur +labs(fill = "Time")+
                       theme(legend.position= "right",
                             legend.justification = "top", 
                             legend.key.width = unit(.45, "cm"),
                             legend.key.height = unit(.45, "cm"),
                             legend.margin =  margin(t=0,b=2))),
          ncol=5, rel_widths = c(1.8,1.8,1.8,1.8,0.6),
          labels = c("A","B","C","D","","E","F","G","H",""))

######### Figure S5











  