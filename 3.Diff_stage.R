
library(ggplot2)
library(tidyverse)
library(ggpubr)


Expression_data <- read.csv("Expression_data.csv",row.names=1)
clinical_data <- read.csv("Clinical_data.csv")
RORA_ex_data <- Expression_data["RORA", ] %>% t() %>% as.data.frame() %>% rownames_to_column()
RORA_ex_data <- merge(RORA_ex_data, clinical_data, by.x="rowname", by.y = "Sample", all.x =TRUE)
RORA_ex_data <- RORA_ex_data %>% select(RORA, Stage) %>% filter(Stage != "<NA>") %>% filter(Stage != "—") %>% 
mutate(Group = case_when(Stage %in% c("I", "IA", "IB", "II", "IIA", "IIB", "IIC") ~ "Low stage", 
                            Stage %in% c("III", "IIIA", "IIIB", "IIIC", "IIID", "IV") ~ "High stage",TRUE ~ NA))


p <- ggplot(data=RORA_ex_data)+ 
  geom_boxplot(mapping=aes(x=Group,y=RORA,colour = Group ,fill = Group), 
               alpha = 0.8,notch=TRUE,
               size=0.9,
               width = 0.6,outlier.colour = NA)+ 
  scale_color_manual(limits=c("Low stage","High stage"), 
                     values=c("#38434e","#38434e"))+
  scale_fill_manual(values=c("#748A9F", "#748A9F"))+ 
  scale_x_discrete(limits=c("Low stage","High stage"))+
  geom_signif(mapping=aes(x=Group,y=RORA), 
              comparisons = list(c("Low stage", "High stage") 
              ),
              map_signif_level=T, 
              tip_length=c(0,0,0,0,0,0,0,0,0,0,0,0), 
              y_position = c(6.5), 
              size=1, 
              textsize = 4, 
              test = "wilcox.test")+ 
  theme_classic(  
    base_line_size = 1 
  )+
  labs(title="Our inhouse data (n=177)",x="",y="Expression")+
    rremove("axis.text")+
  theme(plot.title = element_text(size = 15,
                                  colour = "black",
                                  hjust = 0.5),
        axis.title.y = element_text(size = 15, 
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_text(color="black", 
                                    size=15, 
                                    face="bold"),
        legend.text = element_text(color="black", 
                                   size = 10, 
                                   face = "plain"), 
        axis.text.y = element_text(size = 13,  
                                   color = "black",
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0)
  )
ggsave(filename="box_plot.pdf", plot=p)
