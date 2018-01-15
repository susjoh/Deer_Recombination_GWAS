countIF <- function(x){
  
  x <- data.frame(Var1 = x)
  x$Order <- 1:nrow(x)
  
  y <- data.frame(table(x$Var1))
  
  x <- merge(x, y, all.x = T)
  
  x <- x[with(x, order(Order)), ]
    
  x$Freq
  
}

