val1 <- data.frame(matrix(c(1,2,3,4), ncol=2))

eval(parse(text=paste("val1$X2")))

exp <- parse(text = c("
  x <- 4
  x
  5
"))


asdasod
