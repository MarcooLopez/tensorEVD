library(ggplot2)

# 1. Simulating small matrices
m = 40; n = 75
p = 50; q = 40
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

res1 <- microbenchmark::microbenchmark(
        'kronecker\n(base)'            = kronecker(A,B),
        'kronecker.prod\n(fastmatrix)' = fastmatrix::kronecker.prod(A,B),
        'Kronecker\n(tensorEVD)'       = tensorEVD::Kronecker(A,B),
       times = 5)

tt1 <- paste0("'Kronecker('*A[",m,"*'x'*",n,"]*', '*B[",p,"*'x'*",q,"]*')'==K[",m*p,"*'x'*",n*q,"]")

# 2. Simulating large matrices
m = 100; n = 150
p = 100; q = 100
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

res2 <- microbenchmark::microbenchmark(
       'kronecker\n(base)'            = kronecker(A,B),
       'kronecker.prod\n(fastmatrix)' = fastmatrix::kronecker.prod(A,B),
       'Kronecker\n(tensorEVD)'       = tensorEVD::Kronecker(A,B),
      times = 5)

tt2 <- paste0("'Kronecker('*A[",m,"*'x'*",n,"]*', '*B[",p,"*'x'*",q,"]*')'==K[",m*p,"*'x'*",n*q,"]")

# 3. Making a barplot
res <- rbind(data.frame(group=tt1, expr=res1$expr, time=res1$time/1E9),
             data.frame(group=tt2, expr=res2$expr, time=res2$time/1E9))
dat <- aggregate(time~expr+group, data=res, FUN=mean)
dat$sd <- aggregate(time~expr+group, data=res, FUN=sd)$time

pp <- ggplot(dat, aes(expr,time)) + theme_bw() +
      geom_bar(stat="identity", fill='lightyellow3', width=0.65) +
      geom_errorbar(aes(ymin=time-sd, ymax=time+sd), width=0.2) +
      labs(x=NULL, y="Time (seconds)") +
      facet_wrap(~group, scales="free", labeller=label_parsed)
print(pp)