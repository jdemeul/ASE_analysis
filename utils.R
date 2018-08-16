## wrapper for pbetabinomial ~ binom.test
require(VGAM)

betabinom.test <- function (q, size, prob, rho, alternative = c("two.sided", "less", "greater")) 
{
  PVAL <- switch(alternative, less = pbetabinom(q = q, size = size, prob = prob, rho = rho), greater = 1 - pbetabinom(q = q - 1, size = size, prob = prob, rho = rho),
                 two.sided = { if (prob == 0) (q == 0) else if (prob == 1) (q == size) else {
                   relErr <- 1 + 1e-07
                   d <- dbetabinom(x = q, size = size, prob = prob, rho = rho)
                   m <- size * prob
                   if (q == m) 1 else if (q < m) {
                     i <- seq.int(from = ceiling(m), to = size)
                     y <- sum(dbetabinom(x = i, size = size, prob = prob, rho = rho) <= d * relErr)
                     p1 <- 1 - pbetabinom(q = size - y, size = size, prob = prob, rho = rho)
                     if (p1 < 0)
                       pbetabinom(q = q, size = size, prob = prob, rho = rho)
                     else
                       p1 + pbetabinom(q = q, size = size, prob = prob, rho = rho) 
                   } else {
                     i <- seq.int(from = 0, to = floor(m))
                     y <- sum(dbetabinom(x = i, size = size, prob = prob, rho = rho) <= d * relErr)
                     p1 <- 1 - pbetabinom(q = q - 1, size = size, prob = prob, rho = rho)
                     if (p1 < 0)
                       pbetabinom(q = y - 1, size = size, prob = prob, rho = rho)
                     else 
                       p1 + pbetabinom(q = y - 1, size = size, prob = prob, rho = rho)
                     
                   }
                 }
                 })
}


## wrapper for pbetabinomial ~ binom.test
betabinom.test.ab <- function (q, size, shape1, shape2, alternative = c("two.sided", "less", "greater")) 
{
  PVAL <- switch(alternative, less = pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2), greater = 1 - pbetabinom.ab(q = q - 1, size = size, shape1 = shape1, shape2 = shape2),
                 two.sided = { if (shape1 == 0) (q == 0) else if (shape2 == 0) (q == size) else {
                   relErr <- 1 + 1e-07
                   d <- dbetabinom.ab(x = q, size = size, shape1 = shape1, shape2 = shape2)
                   m <- size * shape1 / (shape1 + shape2)
                   if (q == m) 1 else if (q < m) {
                     i <- seq.int(from = ceiling(m), to = size)
                     y <- sum(dbetabinom.ab(x = i, size = size, shape1 = shape1, shape2 = shape2) <= d * relErr)
                     p1 <- 1 - pbetabinom.ab(q = size - y, size = size, shape1 = shape1, shape2 = shape2)
                     if (p1 < 0)
                       pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2)
                     else
                       p1 + pbetabinom.ab(q = q, size = size, shape1 = shape1, shape2 = shape2) 
                   } else {
                     i <- seq.int(from = 0, to = floor(m))
                     y <- sum(dbetabinom.ab(x = i, size = size, shape1 = shape1, shape2 = shape2) <= d * relErr)
                     p1 <- 1 - pbetabinom.ab(q = q - 1, size = size, shape1 = shape1, shape2 = shape2)
                     if (p1 < 0)
                       pbetabinom.ab(q = y - 1, size = size, shape1 = shape1, shape2 = shape2)
                     else 
                       p1 + pbetabinom.ab(q = y - 1, size = size, shape1 = shape1, shape2 = shape2)
                     
                   }
                 }
                 })
}



ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=".",cex=2, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}


ggqq <- function(pvector, title = "Quantile-quantile plot of p-values", spartan = T) {
  o <- -log10(sort(pvector,decreasing=F))
  e <- -log10( 1:length(o)/length(o) )
  # could use base graphics
  # plot(e,o,pch=19,cex=0.25, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)),ylim=c(0,max(e)))
  # lines(e,e,col="red")
  plot <- qplot(e,o) + geom_abline(intercept = 0,slope = 1, col = "red")
  # plot <- plot + labs(title = title)
  plot <- plot + scale_x_continuous(name = expression(Expected~~-log[10](italic(p))), limits = c(0, max(e)))
  plot <- plot + scale_y_continuous(name = expression(Observed~~-log[10](italic(p))), limits = c(0, max(o)))
  if (spartan) plot <- plot + theme_minimal()
  return(plot)
}
