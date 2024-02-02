library(microbenchmark)
library(tensorEVD)
library(ggplot2)

# 1. Simulating A and B matrices
m = 100; n = 50
p = 50; q = 100
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

dm <- c(nrow(A)*nrow(B), ncol(A)*ncol(B))   # dimension of the Kronecker

# 2. Subsetting an small matrix with 20% of rows/columns
rows1 <- sample(seq(dm[1]), 0.20*dm[1])
cols1 <- sample(seq(dm[2]), 0.20*dm[2])

# 3. Subsetting a large matrix with 200% of rows/columns
rows2 <- sample(seq(dm[1]), 2*dm[1], replace=TRUE)
cols2 <- sample(seq(dm[2]), 2*dm[2], replace=TRUE)

res <- microbenchmark(
       'Kronecker(A,B,rows,cols)_Small submatrix'  = Kronecker(A,B,rows=rows1,cols=cols1),
       'Kronecker(A,B)[rows,cols]_Small submatrix' = Kronecker(A,B)[rows1,cols1],
       'Kronecker(A,B,rows,cols)_Large submatrix'  = Kronecker(A,B,rows=rows2,cols=cols2),
       'Kronecker(A,B)[rows,cols]_Large submatrix' = Kronecker(A,B)[rows2,cols2],
     times = 30)

# 4. Making a barplot
res <- data.frame(expr=res$expr, time=res$time/1E9)
res <- data.frame(res, do.call(rbind,strsplit(as.character(res$expr),"_")))
dat <- aggregate(time~X1+X2, data=res, FUN=median)
dat$X2 <- paste0(dat$X2,": ",ifelse(dat$X2=="Small submatrix",
                                    paste0(length(rows1),"x",length(cols1)),
                                    paste0(length(rows2),"x",length(cols2))))

title <- paste0("'Subsetting of Kronecker('*A[",m,"*'x'*",n,"]*', '*B[",p,"*'x'*",
                q,"]*')'==K[",m*p,"*'x'*",n*q,"]")

pp <- ggplot(dat, aes(X1,time)) + theme_bw() +
      geom_bar(stat="identity", fill='lightyellow3', width=0.6) +
      labs(x=NULL, y="Time (seconds)", title=parse(text=title)) +
      facet_wrap(~X2, scales="free") +
      theme(plot.title=element_text(hjust=0.5))
print(pp)
