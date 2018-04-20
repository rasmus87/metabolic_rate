# Example showing the 

set.seed(42)
n <- 100
x <- rlnorm(n)
hist(x)
hist(log10(x))
y <- (x + x*rnorm(n)/10)^0.75

y2 <- y/x
y3 <- y/x^0.75

logx <- log10(x)
logy <- log10(y)
logy2 <- log10(y2)
logy3 <- log10(y3)

ggplot(data=data.frame(logx, logy, logy2, logy3), aes(x = logx)) + 
  geom_point(aes(y = logy), col = "red") + 
  geom_smooth(method = "lm", aes(y = logy), col = "red") +
  geom_point(aes(y = logy2), col = "green") +
  geom_smooth(method = "lm", aes(y = logy2), col = "green") +
  geom_point(aes(y = logy3), col = "blue") +
  geom_smooth(method = "lm", aes(y = logy3), col = "blue")

lm(logy ~ logx) %>% summary
lm(logy2 ~ logx) %>% summary
lm(logy3 ~ logx) %>% summary

