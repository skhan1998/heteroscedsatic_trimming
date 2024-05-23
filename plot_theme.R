theme_custom <- function(){ 
  #font <- "Georgia"   #assign font family up front
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(-8,0,0,0),      
      axis.title = element_text(size=9.5),
      axis.text = element_text(size=7.5),
      legend.title = element_text(size = 9.5),
      legend.text = element_text(size=7.5),
      strip.text = element_text(size = 9.5, margin = margin(t = 2, b = 2, unit = "pt")))
}


color_map <- c(rgb(31, 119, 180, maxColorValue = 256),
               rgb(255, 127, 14, maxColorValue = 256),
               rgb(44, 160, 44, maxColorValue = 256),
               rgb(214, 39, 40, maxColorValue = 256),
               rgb(148, 103, 189, maxColorValue = 256),
               rgb(140, 86, 75, maxColorValue = 256),               
               rgb(227, 119, 194, maxColorValue = 256))