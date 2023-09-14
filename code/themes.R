library(ggplot2)         #plotting with ggplot2
library(sysfonts)        #adding fonts 


theme_set(theme_minimal(base_size = 12, base_family = "Satoshi") +
            theme(plot.title = element_text(hjust = 0.5,size = 25, face = "bold"), 
                  legend.title.align=0.5,
                  legend.position = "bottom",
                  #panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
                  axis.ticks.x = element_line(colour = "#333333"),
                  axis.ticks.y = element_line(colour = "#333333"),
                  axis.ticks.length =  unit(0.26, "cm"),
                  panel.border = element_rect(colour="#333333", fill=NA), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.text=element_text(size=20),
                  axis.title=element_text(size=23),
                  legend.text = element_text(size=20),
                  legend.title = element_text(size=20)
            )) 



theme_set(theme_minimal(base_size = 12, base_family = "Roboto") +
            theme(plot.title = element_text(hjust = 0.5,size = 25, face = "bold"), 
                  legend.title.align=0.5,
                  legend.position = "bottom",
                  #panel.grid.major  = element_line(colour = "#F0F0F2", size=0.5),
                  axis.ticks.x = element_line(colour = "#333333"),
                  axis.ticks.y = element_line(colour = "#333333"),
                  axis.ticks.length =  unit(0.26, "cm"),
                  panel.border = element_rect(colour="#333333", fill=NA), 
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.text=element_text(size=20),
                  axis.title=element_text(size=23),
                  legend.text = element_text(size=20),
                  legend.title = element_text(size=20)
            )) 


theme_pres <- function(
    base_size =12,
    base_line_size = base_size/24,
    base_rect_size = base_size/24,
    base_family = "Roboto"
){
  theme_bw(
    base_size = base_size,
    base_family = base_family,
    base_line_size = base_line_size,
    base_rect_size = base_rect_size
  ) %+replace%
    theme(
      axis.text.x = element_text(size = 15, color = "grey30", margin=margin(t=5,b=10)),
      axis.text.y = element_text(size = 15, color = "grey30", margin=margin(t=0, r = 5, b = 0, l = 10)),
      axis.ticks = element_line(color = "grey91", size = .5),
      axis.ticks.length.x = unit(.5, "lines"),
      axis.ticks.length.y = unit(.7, "lines"),
      axis.title=element_text(size=18),
      legend.text=element_text(size=15),
      legend.title=element_text(size=18, hjust = 0.5),
      panel.grid.major = element_line(color=NA),
      plot.margin = margin(20, 40, 20, 40),
      plot.background = element_rect(fill = "grey98", color = "grey98"),
      panel.background = element_rect(fill = "grey98", color = "grey98"),
      legend.key= element_rect(fill = "grey98", colour = "grey98"),
      legend.background= element_rect(fill = "grey98", colour = "grey98"),
      plot.title = element_text(color = "grey10", size = 20, face = "bold",
                                margin = margin(t = 25,b=10), hjust=0.5),
      plot.subtitle = element_markdown(color = "grey30", size = 12, 
                                       lineheight = 1.35, hjust=0.5), 
      #margin = margin(t = 15, b = 40)), # i have no idea what is it doing 
      plot.caption = element_text(color = "grey30", size = 10,
                                  lineheight = 1.2, hjust = 0, 
                                  margin = margin(t = 80)),
      #legend.position = "none"
  complete=TRUE
    )
  }

