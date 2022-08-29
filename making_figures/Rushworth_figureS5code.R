library(tidyverse)
library(cowplot)
library(ggthemes)
library(RColorBrewer)
library(wesanderson)
library(gridExtra)
#source("/Users/catherinerushworth1/projects/gftheory_sandbox/functions_forsummaries.R")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[leg]
  return(legend)
}



### Supp Figure for marker order


### make summary files using skip bigass



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
ramf_sum <- read_csv("AMFfinescaler_sum.csv")
rmfa_sum <- read_csv("MFAfinescaler_sum.csv")
rmaf_sum <- read_csv("MAFfinescaler_sum.csv")



### AMF


tmp <- ramf_sum %>% 
  mutate(z =  as.numeric(as.factor(r12)) ) %>% 
  filter(r12 %in% c(10^(-5:-1))) %>% 
  group_by(r12,z) %>%
  tally()%>% 
  ungroup()%>%
  mutate(breaks = z, cM = 100 * r12)




scale_x_ramf <- scale_x_continuous(breaks = pull(tmp,breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` AM        `])),
                                   lim = c(3,17))

scale_y_ramf <- scale_y_continuous(breaks = pull(tmp,breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` MF        `])),
                                   lim = c(3,17))

ramf_plot_format_sub <-  ramf_sum %>% 
  filter(s == 0.8, m0 == 0.01)%>%
  mutate(r12 = as.numeric(as.factor(r12)),
         r23 = as.numeric(as.factor(r23)))


ramf_amount <- ggplot(ramf_plot_format_sub,
                      aes(x = r12, y = r23, fill = max.reinf_1))+
  geom_tile(show.legend = FALSE)+
  amount_scale+
  theme_tufte(base_family="Helvetica") +
  labs(title="AMF")+
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5))+
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
              pull(color), fill = "#1d5a4c")




### MFA

tmp_rmfa <- rmfa_sum %>% 
  mutate(z =  as.numeric(as.factor(r12)) ) %>% 
  filter(r12 %in% c(10^(-5:-1))) %>% 
  group_by(r12,z) %>%
  tally()%>% 
  ungroup()%>%
  mutate(breaks = z, cM = 100 * r12)


scale_x_rmfa <- scale_x_continuous(breaks = pull(tmp_rmfa, breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` MF        `])),
                                   lim = c(3,17))

scale_y_rmfa <- scale_y_continuous(breaks = pull(tmp_rmfa,breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` FA        `])),
                                   lim = c(3,17))


rmfa_plot_format_sub <-  rmfa_sum %>% 
  filter(s == 0.8, m0 == 0.01)%>%
  mutate(r12 = as.numeric(as.factor(r12)),
         r23 = as.numeric(as.factor(r23)))


rmfa_amount <- ggplot(rmfa_plot_format_sub ,
                      aes(x = r12, y = r23, fill = max.reinf_1))+
  geom_tile(show.legend = FALSE)+
  amount_scale+
  theme_tufte(base_family="Helvetica")+
  labs(title="MFA")+
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5))+
  scale_x_rmfa+
  scale_y_rmfa

rmfa_dur <- ggplot(rmfa_plot_format_sub ,
                   aes(x = r12, y = r23,fill = reinf.time))+
  geom_tile(show.legend = FALSE)+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  scale_x_rmfa+
  scale_y_rmfa+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = rmfa_plot_format_sub  %>%
              filter(m0==0.01)%>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "#1d5a4c")



### MAF

tmp_rmaf <- rmaf_sum %>% 
  mutate(z =  as.numeric(as.factor(r12)) ) %>% 
  filter(r12 %in% c(10^(-5:-1))) %>% 
  group_by(r12,z) %>%
  tally()%>% 
  ungroup()%>%
  mutate(breaks = z, cM = 100 * r12)



scale_x_rmaf <- scale_x_continuous(breaks = pull(tmp_rmaf, breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` AM        `])),
                                   lim = c(3,17))

scale_y_rmaf <- scale_y_continuous(breaks = pull(tmp_rmaf,breaks), 
                                   labels = c(expression(10^{-5}),expression(10^{-4}),expression(10^{-3}),expression(10^{-2}),expression(10^{-1})),
                                   expand = c(0,0), name = expression(italic(r[` FA        `])),
                                   lim = c(3,17))



rmaf_plot_format_sub <-  rmaf_sum %>% 
  filter(s == 0.8, m0 == 0.01)%>%
  mutate(r12 = as.numeric(as.factor(r12)),
         r23 = as.numeric(as.factor(r23)))

g2 <- function(a.gplot){
  if (!gtable::is.gtable(a.gplot))
    a.gplot <- ggplotGrob(a.gplot)
  gtable::gtable_filter(a.gplot, 'guide-box', fixed=TRUE)
}

rmaf_amount <- ggplot(rmaf_plot_format_sub ,
                      aes(x = r23, y = r12, fill = max.reinf_1))+
  geom_tile()+
  amount_scale+
  theme_tufte(base_family="Helvetica") +
  labs(title="MAF",fill="Amount")+
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5))+
        #plot.margin=unit(c(0.2,-0.5,0.2,0.2),"cm"))+
  # theme(legend.position = "right", legend.key.width = unit(.3, "cm"), legend.key.height = unit(.7, "cm"),
  #       legend.margin=margin(t=3,r=0,b=7,l=0), #legend.title = element_text(size = 8),
  #       legend.box.margin = margin(6, 7, 6, 0))+
  scale_x_rmaf+
  scale_y_rmaf 




rmaf_dur <- ggplot(rmaf_plot_format_sub ,
                   aes(x = r23, y = r12,fill = reinf.time))+
  geom_tile()+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  labs(fill="Time")+
  theme(legend.position = "right", legend.key.width = unit(.3, "cm"), legend.key.height = unit(.7, "cm"),
        legend.margin=margin(t=3,b=7))+
  scale_x_rmaf+
  scale_y_rmaf+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = rmaf_plot_format_sub  %>%
              filter(m0==0.01)%>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "#1d5a4c")



leg1 <- cowplot::get_legend(rmaf_amount)
leg2 <- cowplot::get_legend(rmaf_dur)
rmaf_amount <- rmaf_amount + theme(legend.position="none")
rmaf_dur <- rmaf_dur + theme(legend.position="none")
#ggplot(leg)


### Assemble with separate legends
testrun <- plot_grid(ramf_amount,
                     rmfa_amount,
                     rmaf_amount,
                     leg1,
                     ramf_dur,
                     rmfa_dur,
                     rmaf_dur,
                     leg2,
                     # get_legend(rmaf_dur +labs(fill = "Time")+
                     #              theme(legend.position= "right",
                     #                    legend.justification = "top", 
                     #                    legend.key.width = unit(.45, "cm"),
                     #                    legend.key.height = unit(.45, "cm"),
                     #                    legend.margin =  margin(t=0,b=2))),
                     ncol=4, rel_widths = c(1.8,1.8,1.8,0.6), rel_heights = c(1.8,1.6),
                     labels = c("A","B","C","","D","E","F",""))



### Assemble

quartz()
testrun <- plot_grid(ramf_amount,
          rmfa_amount,
          rmaf_amount,
          ramf_dur,
          rmfa_dur,
          rmaf_dur,
          ncol=3, rel_widths = c(1.8,1.8,2.1), rel_heights = c(1.8,1.6),
          labels = c("A","B","C","D","E","F"))



MO <- tempfile("markerorder_amtdur", fileext = ".pdf")
#embed_fonts("markerorder_amtdur.pdf", outfile="font_markerorder_amtdur.pdf")
save_plot(testrun, filename="markerorder_amtdur_2.pdf", ncol=2, base_asp=1)







### no legend

rmaf_amount <- ggplot(rmaf_plot_format_sub ,
                      aes(x = r23, y = r12, fill = max.reinf_1))+
  geom_tile(show.legend = FALSE)+
  amount_scale+
  theme_tufte(base_family="Helvetica") +
  labs(title="MAF")+
  theme(plot.title = element_text(family="Palatino",face="italic",size=16,hjust=0.5))+
  #plot.margin=unit(c(0.2,-0.5,0.2,0.2),"cm"))+
  # theme(legend.position = "right", legend.key.width = unit(.3, "cm"), legend.key.height = unit(.7, "cm"),
  #       legend.margin=margin(t=3,r=0,b=7,l=0), #legend.title = element_text(size = 8),
  #       legend.box.margin = margin(6, 7, 6, 0))+
  scale_x_rmaf+
  scale_y_rmaf 

rmaf_dur <- ggplot(rmaf_plot_format_sub ,
                   aes(x = r23, y = r12,fill = reinf.time))+
  geom_tile(show.legend=FALSE)+
  time_scale+
  theme_tufte(base_family="Helvetica")+
  labs(fill="Time")+
  theme(legend.position = "right", legend.key.width = unit(.3, "cm"), legend.key.height = unit(.7, "cm"),
        legend.margin=margin(t=3,b=7))+
  scale_x_rmaf+
  scale_y_rmaf+
  geom_tile(aes(fill = is.infinite(reinf.time)), 
            alpha = rmaf_plot_format_sub  %>%
              filter(m0==0.01)%>%
              mutate(color = ifelse(is.infinite(reinf.time),1,0))%>% 
              pull(color), fill = "#1d5a4c")
