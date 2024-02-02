library(microbenchmark)
library(tensorEVD)
library(fastmatrix)
library(ggplot2)

# 1. Simulating incompatible matrices
m = 50; n = 60
p = 40; q = 50
A <- matrix(rnorm(m*n), ncol=n)
B <- matrix(rnorm(p*q), ncol=q)

# 2. Making Hadamard for a small matrix
dm <- c(2000, 2000)  # dimension of the Hadamard

rowsA <- sample(seq(nrow(A)), dm[1], replace=TRUE)
colsA <- sample(seq(ncol(A)), dm[2], replace=TRUE)

rowsB <- sample(seq(nrow(B)), dm[1], replace=TRUE)
colsB <- sample(seq(ncol(B)), dm[2], replace=TRUE)

res1 <- microbenchmark(
         'A[rowsA,colsA]*B[rowsB,colsB]' = A[rowsA,colsA]*B[rowsB,colsB],
         'Hadamard(A,B,rowsA,rowsB,colsA,colsB)' = Hadamard(A,B,rowsA,rowsB,colsA,colsB),
       times = 5)

tt1 <- paste0("Small Hadamard: ",length(rowsA)," x ",length(colsA))

# 3. Making Hadamard for a large matrix
dm <- c(12000, 12000)  # dimension of the Hadamard

rowsA <- sample(seq(nrow(A)), dm[1], replace=TRUE)
colsA <- sample(seq(ncol(A)), dm[2], replace=TRUE)

rowsB <- sample(seq(nrow(B)), dm[1], replace=TRUE)
colsB <- sample(seq(ncol(B)), dm[2], replace=TRUE)

res2 <- microbenchmark(
         'A[rowsA,colsA]*B[rowsB,colsB]' = A[rowsA,colsA]*B[rowsB,colsB],
         'Hadamard(A,B,rowsA,rowsB,colsA,colsB)' = Hadamard(A,B,rowsA,rowsB,colsA,colsB),
       times = 5)

tt2 <- paste0("Large Hadamard: ",length(rowsA)," x ",length(colsA))

# 3. Making a barplot
res <- rbind(data.frame(group=tt1, expr=res1$expr, time=res1$time/1E9),
             data.frame(group=tt2, expr=res2$expr, time=res2$time/1E9))
dat <- aggregate(time~expr+group, data=res, FUN=median)

title <- bquote('Hadamard between submatrices of '*A[.(m)*'x'*.(n)]*' and '*B[.(p)*'x'*.(q)])

pp <- ggplot(dat, aes(expr,time,fill=expr)) + theme_bw() +
      geom_bar(stat="identity", width=0.6) +
      labs(x=NULL, y="Time (seconds)", title=title, fill=NULL) +
      facet_wrap(~group, scales="free") +
      theme(plot.title=element_text(hjust=0.5),
            axis.text.x=element_blank(), legend.position="bottom",
            legend.spacing.x=unit(20,'pt'),
            legend.text=element_text(margin=margin(l=-18)))
print(pp)
