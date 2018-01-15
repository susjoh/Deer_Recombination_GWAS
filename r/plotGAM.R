
plotGAM <- function(fit, grep.value = NULL, make.plot = TRUE, return.data = FALSE, plottitle = ""){
  
  par(mfrow = c(3,3))
  x <- plot(fit)
  dev.off()
  
  gamplot <- NULL
  
  for(i in 1:length(x)) gamplot <- rbind(gamplot, data.frame(x = x[[i]]$x, y = x[[i]]$fit[,1], se = x[[i]]$se, ylab = x[[i]]$ylab))
  
  if(!is.null(grep.value)){
    
    gamplot <- gamplot[grep(grep.value, gamplot$ylab),]
    
  }
  
  if(make.plot == TRUE){
    print(ggplot(gamplot, aes(x, y, colour = ylab)) +
            geom_line(size = 1) +
            geom_line(aes(x, y+se), linetype = "dotted", size = 1) +
            geom_line(aes(x, y-se), linetype = "dotted", size = 1) +
            scale_colour_brewer(palette = "Set1") +
            theme(axis.text.x  = element_text (size = 12),
                  axis.text.y  = element_text (size = 12),
                  strip.text.x = element_text (size = 12),
                  axis.title.y = element_text (size = 14, angle = 90),
                  axis.title.x = element_text (size = 14),
                  strip.background = element_blank(),
                  legend.position = "top",
                  legend.direction = "vertical") +
            ggtitle(plottitle)
          )
  }
  
  if(return.data == TRUE) return(gamplot)
}




plotGAMalt <- function(fit, grep.value = NULL, make.plot = TRUE, return.data = FALSE, plottitle = ""){
  
  par(mfrow = c(3,3))
  x <- plot(fit)
  dev.off()
  
  gamplot <- NULL
  
  for(i in 1:length(x)) gamplot <- rbind(gamplot, data.frame(x = x[[i]]$x, y = x[[i]]$fit[,1], se = x[[i]]$se, ylab = x[[i]]$ylab))
  
  if(!is.null(grep.value)){
    
    gamplot <- gamplot[grep(grep.value, gamplot$ylab),]
    
  }
  
  if(make.plot == TRUE){
    ggplot(gamplot, aes(x, y, colour = ylab)) +
            geom_line() +
            geom_line(aes(x, y+se), linetype = "dotted") +
            geom_line(aes(x, y-se), linetype = "dotted") +
            scale_colour_brewer(palette = "Set1") +
            theme(axis.text.x  = element_text (size = 12),
                  axis.text.y  = element_text (size = 12),
                  strip.text.x = element_text (size = 12),
                  axis.title.y = element_text (size = 14, angle = 90),
                  axis.title.x = element_text (size = 14),
                  strip.background = element_blank(),
                  legend.position = "top",
                  legend.direction = "vertical") +
            ggtitle(plottitle)
    
  }
  
  if(return.data == TRUE) return(gamplot)
}