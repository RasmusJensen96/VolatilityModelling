plot.FitModel <- function(obj){
    require(ggplot2)
  attach(obj)
  Filtered <- Filtered
  Obs    <- Observations
  Dates  <- Dates
    p =  ggplot() + geom_line(aes(x = Dates, y = Filtered, colour = obj$Targetfilter)) +
      geom_point(aes(x = Dates, y = Obs,colour = "Observed returns"),size = 0.5) +
      scale_color_manual(values = c(rgb(0,0.4470,0.7410),rgb(0.6350, 0.0780,0.1840))) +
      labs(colour=" ",
           x="Dates",
           y="Returns")

 if (exists("Name")){
   p = p + ggtitle(paste0(Model,"-",Dist, " fit"), paste0("Series: ", Name)) + theme_bw() + theme(legend.position=c(0.80, 1.04),
                                                                                  legend.background = element_blank(),
                                                                                  legend.direction='horizontal')
   } else {
  p =  p + ggtitle(paste0(Model,"-",Dist, " fit")) + theme_bw() + theme(legend.position=c(0.80, 1.04),
                                                         legend.background = element_blank(),
                                                         legend.direction='horizontal')
  }
detach(obj)
 p
}
