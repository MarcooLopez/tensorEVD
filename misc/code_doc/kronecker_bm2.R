library(ggplot2)

# 1. Simulating A and B matrices
m = 100; n = 50
p = 50; q = 100
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

dm <- c(m*p, n*q)   # dimension of the Kronecker

# 2. Subsetting a matrix with 10% of rows/columns
rows1 <- sample(seq(dm[1]), 0.10*dm[1], replace=TRUE)
cols1 <- sample(seq(dm[2]), 0.10*dm[2], replace=TRUE)

# 3. Subsetting a matrix with 200% of rows/columns
rows2 <- rep(seq(dm[1]), 2)
cols2 <- rep(seq(dm[2]), 2)

res <- microbenchmark::microbenchmark(
      'Kronecker(A,B,rows,cols)_Small submatrix'  = tensorEVD::Kronecker(A,B,rows=rows1,cols=cols1),
      'Kronecker(A,B)[rows,cols]_Small submatrix' = tensorEVD::Kronecker(A,B)[rows1,cols1],
      'Kronecker(A,B,rows,cols)_Large submatrix'  = tensorEVD::Kronecker(A,B,rows=rows2,cols=cols2),
      'Kronecker(A,B)[rows,cols]_Large submatrix' = tensorEVD::Kronecker(A,B)[rows2,cols2],
     times = 10)

# 4. Making a barplot
res <- data.frame(expr=res$expr, time=res$time/1E9)
res <- data.frame(res, do.call(rbind,strsplit(as.character(res$expr),"_")))
dat <- aggregate(time~X1+X2, data=res, FUN=mean)
dat$sd <- aggregate(time~X1+X2, data=res, FUN=sd)$time
dat$X2 <- gsub("Small",paste0(length(rows1),"x",length(cols1)),dat$X2)
dat$X2 <- gsub("Large",paste0(length(rows2),"x",length(cols2)),dat$X2)

title <- paste0("'Subsetting of Kronecker('*A[",m,"*'x'*",n,"]*', '*B[",p,"*'x'*",
                q,"]*')'==K[",m*p,"*'x'*",n*q,"]")

pp <- ggplot(dat, aes(X1,time)) + theme_bw() +
      geom_bar(stat="identity", fill='lightyellow3', width=0.6) +
      geom_errorbar(aes(ymin=time-sd, ymax=time+sd), width=0.1) +
      labs(x=NULL, y="Time (seconds)", title=parse(text=title)) +
      facet_wrap(~X2, scales="free") +
      theme(plot.title=element_text(hjust=0.5))
print(pp)
