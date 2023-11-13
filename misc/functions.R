
#=========================================================
melt_data <- function(data, id, measure,
                      variable.name = "variable",
                      value.name = measure){

  stopifnot(length(measure)==length(value.name))

  names0 <- colnames(data)[!colnames(data) %in% id]

  levels0 <- lapply(measure, function(x){
    gsub(paste0("^",x), "", grep(paste0("^",x),names0,value=T))
  })
  var.levels <- unique(unlist(levels0))

  # Checkpoint
  tmp <- paste(paste0("'",var.levels,"'"),collapse=",")
  for(k in 1:length(levels0)){
    flag <- length(var.levels) == length(levels0[[k]])
    if(!(flag & all(levels0[[k]]%in%var.levels) & all(var.levels%in%levels0[[k]]))){
      stop("Variables for measure '",measure[k],"' must be equal to ",tmp)
    }
  }

  if(requireNamespace("reshape2")){
    out <- lapply(seq_along(measure), function(k){
      tmp <- grep(paste0("^",measure[k]),names0,value=T)
      dt <- data[,c(id, tmp)]
      colnames(dt) <- c(id, gsub(paste0("^",measure[k]), "", tmp))
      dt <- dt[,c(id,var.levels)]
      reshape2::melt(dt, id=id, variable.name=variable.name,value.name=value.name[k])
    })
  }
  #lapply(out, head)
  tmp <- do.call(cbind,lapply(out, function(dt){
     dt[,!colnames(dt) %in% c(id,variable.name),drop=F]
  }))

  data.frame(out[[1]][,c(id,variable.name),drop=F],tmp)
}

#=========================================================
# Capitalize a string:  "hello word" -> "Hello word"
#=========================================================
capitalize <- function(string){
  substr(string,1,1) <- toupper(substr(string,1,1))
  string
}

#=========================================================
# Get the total variance of the term Xb from MCMC
# samples using a matrix n x p matrix B=[b_1',...,b_n']'
# containing at each row, the n samples for each predictor
#=========================================================
get_variance <- function(X, B){
  stopifnot(ncol(X) == ncol(B))
  res <- rep(NA, nrow(B))
  tmp <- tcrossprod(X, B)
  for(k in 1:nrow(B)){
    res[k] <- var(tmp[,k])
  }
  c(mean=mean(res), SD=sd(res))
}

#=========================================================
# This function is a modification from the one in BGLR
# with the difference that reads incomplete data.
# If the last row (say n) is incomplete, it fixes the
# output matrix to be (n-1) x p. Otherwise n x p
# where p is the number of predictors
#=========================================================
readBinMat <- function(filename, byrow = TRUE, storageMode = "double")
{
    if(!storageMode %in% c("single", "double")) {
        stop("storageMode can either be 'single' or 'double' (default)")
    }
    size <- ifelse(storageMode=="single", 4, 8)
    fileIn <- gzfile(filename, open="rb")
    n <- readBin(fileIn, n=1, what=numeric(), size=size)
    p <- readBin(fileIn, n=1, what=numeric(), size=size)
    tmp <- readBin(fileIn, n=(n*p), what=numeric(), size=size)
    if(length(tmp) < (n*p)){
      n2 <- floor(length(tmp)/p)
      tmp <- tmp[seq(n2*p)]
      X <- matrix(data=tmp, nrow=n2, ncol=p, byrow=byrow)
      message("Only ",n2," complete rows (of ",n,") were read")
    }else{
      X <- matrix(data=tmp, nrow=n, ncol=p, byrow=byrow)
    }
    close(fileIn)
    return(X)
}

#=========================================================
# Function to get within year_loc correlation
#=========================================================
get_corr <- function(dat, y_name = "y",
                     yHat_name = "yHat",
                     by = "year_loc"){
  res <- do.call(rbind,lapply(split(dat, dat[,by]),function(x){
         x$y <- x[,y_name]
         x$yHat <- x[,yHat_name]
         cond1 <- (!is.na(x$y))
         cond2 <- (!is.na(x$yHat))
         nRecords <- sum(cond1&cond2)
         correlation <- cor(x$y,x$yHat,use='pairwise.complete')
         data.frame(x[1,c(by),drop=F], correlation,
                    nRecords,
                    SE=sqrt((1-correlation^2)/(nRecords-2))
         )}))
  rownames(res) <- NULL
  res
}

#=========================================================
# ggplot2-like colour scale in HCL space
#=========================================================
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#=========================================================
# Bar, line plots
#=========================================================
make_plot <- function(data, x, y, group, SD = NULL, nSD = 1,
                      by = NULL, by.label = by,
                      facet = NULL, type = c("bar","line"),
                      group.label = capitalize(group),
                      xlab = x, ylab = y, ylim = c(NA,NA),
                      pymin = 0.4, hline = NULL, hline.color = "blue",
                      text = NULL, text.color = "gray20", text.size = 2.8,
                      text.just = 0, text.pos = 0.05, text.scales = "fixed",
                      bar.width = 1, line.width = 0.5,
                      sec.axis = FALSE, ylab2 = NULL,
                      group.color = NULL, by.color = NULL,
                      errorbar.size = line.width, errorbar.width = 0.35,
                      by.alpha = 0.4, expand.y = c(0.03,0.03),
                      gap.by = 0.3, gap.x = 0.15,
                      rect.by.color = "gray60", rect.by.height = Inf,
                      ylabels = NULL, scales = "fixed")
{
  type <- match.arg(type)
  scales <- match.arg(scales, choices=c("fixed","free_x","free_y","free"))
  text.scales <- match.arg(text.scales, choices=c("fixed","free"))

  names0 <- c(facet, x, by, group)
  fm <- paste0(y,' ~ ',paste(names0,collapse="+"))
  DF <- aggregate(formula(fm), data=data, FUN=mean)
  colnames(DF)[colnames(DF)==y] <- "mean"
  DF$length <- aggregate(formula(fm), data=data, FUN=length)[,y]
  if(is.null(SD)){ # Take the SD of the y
    DF$sd <- nSD*aggregate(formula(fm), data=data, FUN=sd)[,y]
  }else{ # Take the mean of SD
    fm <- paste0(SD,' ~ ',paste(names0,collapse="+"))
    DF$sd <- nSD*aggregate(formula(fm), data=data, FUN=mean)[,SD]
  }

  DF$x <- DF[,x]
  DF$group <- DF[,group]
  if(is.null(by)){
    DF$by <- factor("by")
  }else{
    DF$by <- DF[,by]
  }
  if(is.null(facet)){
    DF$facet <- factor("facet")
  }else{
    DF$facet <- DF[,facet]
  }
  tmp <- setdiff(names0, c("x","group","facet","by"))
  if(length(tmp)>0){
    DF <- DF[,!colnames(DF) %in% tmp]
  }

  if(!is.factor(DF$x)){
    DF$x <- factor(DF$x)
  }
  if(!is.factor(DF$group)){
    DF$group <- factor(DF$group)
  }
  if(!is.factor(DF$facet)){
    DF$facet <- factor(DF$facet)
  }
  if(!is.factor(DF$by)){
    DF$by <- factor(DF$by)
  }

  bw <- (1-gap.x)/nlevels(DF$group)
  DF$x_by <- DF$by:DF$x
  DF$x0 <- as.numeric(DF$x_by) + gap.by*(as.numeric(DF$by)-1)
  DF$xm <- DF$x0
  if(type=="bar"){
    DF$xm <- (DF$x0-0.5) + (gap.x/2) + bw*(as.numeric(DF$group)-1) + bw/2
  }
  DF$x1 <- DF$xm - bar.width*bw/2
  DF$x2 <- DF$xm + bar.width*bw/2

  # Limits for y
  DF <- do.call(rbind,lapply(split(DF, DF$facet),function(df){
    rr <- c(min(df$mean-df$sd), max(df$mean+df$sd))
    tmp <- rr[1] - pymin*diff(rr)
    df$y1 <- ifelse(is.na(ylim[1]), ifelse(tmp<0,0,tmp), ylim[1])
    df$ymax <- ifelse(is.na(ylim[2]), rr[2], ylim[2])
    df$ry <- diff(c(df$y1[1],rr[2]))
    df
  }))
  rownames(DF) <- NULL

  errorbar.width <- errorbar.width*bw
  if(rect.by.height < 0){
    expand.y[1] <- 0
  }

  #  Data for rectangles if separated by 'by' parameter
  ry <- diff(c(min(DF$y1), max(DF$mean+DF$sd)))
  datby <- do.call(rbind,lapply(split(DF,paste(DF$facet,DF$by)),function(df){
    df <- df[df$x==levels(DF$x)[1] & df$group==levels(DF$group)[1],]
    df$x1 <- df$x0 - (1+gap.by)/2
    df$x2 <- df$x0 + (1+gap.by)/2 + nlevels(DF$x) - 1
    df$y2 <- df$y1 + rect.by.height*ifelse(scales %in% c("free","free_y"),df$ry,ry)
    df
  }))
  rownames(datby) <- NULL
  datby[datby$by==levels(DF$by)[1],'x1'] <- -Inf
  datby[datby$by==levels(DF$by)[nlevels(DF$by)],'x2'] <- Inf

  datby2 <- datby[datby$by %in% levels(DF$by)[-nlevels(DF$by)],]

  # Data for labels
  flagtext <- FALSE
  if(!is.null(text)){
    if(is.function(text)){
      funtext <- text
      flagtext <- TRUE
    }else{
      flagtext <- as.logical(text)
      funtext <-  function(x) sprintf("%.3f", x)
    }
  }

  text.just <- as.integer(text.just)
  if(text.scales=="free"){
    DF$ytext <- NA
    for(i in 1:nrow(DF)){
      rr <- c(DF$y1[i], DF$mean[i])
      tmp <- ifelse(text.just==0, rr[1], rr[2])
      DF$ytext[i] <- tmp + ((-1)^text.just)*text.pos*diff(rr)
    }
  }else{
    if(text.just==0){
      tmp <- DF$y1
    }else{
      tmp <- DF$ymax
    }
    DF$ytext <- tmp + ((-1)^text.just)*text.pos*DF$ry
  }
  hjust <- ifelse(text.just==0, 0, 1)

  breaks0 <- unlist(lapply(split(DF,DF$x_by),function(x)x$x0[1]))
  labels0 <- unlist(lapply(strsplit(names(breaks0),":"),function(x)x[2]))

  if(is.null(group.color)){
    group.color <- gg_color_hue(nlevels(DF$group))
    names(group.color) <- levels(DF$group)
  }
  stopifnot(all(levels(DF$group) %in% names(group.color)))

  if(is.null(by.color)){
    if(requireNamespace("RColorBrewer")){
      by.color <- RColorBrewer::brewer.pal(9, name="Blues")[seq(2,9,by=2)][1:nlevels(DF$by)]
    }
    names(by.color) <- levels(DF$by)
  }
  stopifnot(all(levels(DF$by) %in% names(by.color)))

  xmin <- min(DF$x1)
  xmax <- max(DF$x2)
  #xmin <- 1-(1+gap.by)/2
  #xmax <- nlevels(DF$x_by) + gap.by*(nlevels(DF$by)-1) + (1+gap.by)/2

  if(requireNamespace("ggplot2")){
    if(is.null(ylabels)){
      ylabels <- ggplot2:: waiver()
    }

    theme0 <- ggplot2::theme(
                    strip.text.x = ggplot2::element_text(size=7.5, margin=ggplot2::margin(t=1.2,b=1.2)),
                    legend.justification=c("top"),
                    panel.grid.minor.x=ggplot2::element_blank(),
                    panel.grid.major.x=ggplot2::element_blank(),
                    panel.grid.minor.y=ggplot2::element_blank(),
                    legend.key.size=ggplot2::unit(0.85,"line"),
                    legend.margin=ggplot2::margin(0),
                    legend.box.margin=5*ggplot2::margin(1,1,1,-1),
                    legend.title = ggplot2::element_text(size=9),
                    legend.text = ggplot2::element_text(size=7),
                    axis.text = ggplot2::element_text(size=7),
                    axis.title = ggplot2::element_text(size=9)
                  )

    # Transform if secondary axis
    tt <- unlist(lapply(split(DF, DF$group),function(z){
      max(z$mean+z$sd)
    }))
    ylim.prim <- c(min(DF$y1), tt[levels(DF$group)[1]])
    ylim.sec <- c(min(DF$y1), tt[levels(DF$group)[2]])
    b <- diff(ylim.prim)/diff(ylim.sec)
    a <- ylim.prim[1] - b*ylim.sec[1]
    sec_axis0 <- ggplot2:: waiver()
    if(sec.axis){
      if(nlevels(DF$group)>2){
        stop("Data has more than two levels of 'group'")
      }
      tmp <- levels(DF$group)[2]
      DF[DF$group==tmp,'mean'] <- a + b*DF[DF$group==tmp,'mean']
      DF[DF$group==tmp,'sd'] <- a + b*DF[DF$group==tmp,'sd']

      theme0 <- theme0 +
                ggplot2::theme(axis.title.y = ggplot2::element_text(size=10, color=group.color[1], face='bold'),
                      axis.title.y.right = ggplot2::element_text(size=10, color=group.color[2], face='bold')
                )

      sec_axis0 <- ggplot2::sec_axis(~ (. - a)/b, name = ylab2)
    }

    # Making the plot
    pp <- ggplot2::ggplot(DF, ggplot2::aes(y=mean,x=xm)) +
          ggplot2::scale_x_continuous(breaks=breaks0, labels=labels0, limits=c(xmin,xmax),
                                      expand=ggplot2::expansion(mult = c(0.02,0.02)))
    if(!is.null(by)){
      pp <- pp +
          ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=by),
                             data=datby, color=NA, alpha=by.alpha) +
          ggplot2::scale_fill_manual(by.label, values=by.color) +
          ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(size=1.5, color="gray35"))) +
          ggnewscale::new_scale_fill() +
          ggplot2::geom_vline(ggplot2::aes(xintercept=x2), data=datby2, color=rect.by.color)
    }
    if(type == "bar"){
      pp <- pp + ggplot2::geom_rect(ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=mean, fill=group)) +
                 ggplot2::scale_fill_manual(values=group.color)
    }
    if(type == "line"){
      pp <- pp + ggplot2::geom_line(ggplot2::aes(group=by:group, color=group), linewidth=line.width) +
                 ggplot2::geom_point(ggplot2::aes(group=by:group, color=group)) +
                 ggplot2::scale_color_manual(values=group.color)
    }
    pp <- pp +
          ggplot2::geom_errorbar(ggplot2::aes(group=x_by:group, ymin=mean-sd, ymax=mean+sd),
                                 linewidth=errorbar.size, width=errorbar.width) +
          ggplot2::theme_bw() +
          ggplot2::labs(x=xlab, y=ylab, fill=group.label, color=group.label) +
          ggplot2::scale_y_continuous(labels=ylabels,
                                      expand=ggplot2::expansion(mult = expand.y),
                                      sec.axis=sec_axis0) + theme0
     if(flagtext){
       pp <- pp +
            ggplot2::geom_text(ggplot2::aes(y=ytext, label=funtext(mean)),
                size=text.size, vjust=0.5, hjust=hjust, angle=90, color=text.color)
     }

     if(!is.null(facet)){
        pp <- pp + ggplot2::facet_wrap(~facet, scales=scales, labeller=ggplot2::label_parsed)
     }

     if(!is.null(hline)){
        pp <- pp + ggplot2::geom_hline(yintercept=hline, linetype="dashed", color=hline.color)
     }

     pp <- pp +
           ggplot2::guides(fill = ggplot2::guide_legend(override.aes=list(size = 1.5,
                                                                          color = "gray35")))
   }
   pp
}
