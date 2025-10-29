hits <- function(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE, plotting=TRUE){
  
  bagi <- function(gwasm, nmar=10, threshold=1, pick=FALSE, method="cluster", only.mark=FALSE, plotting=TRUE){
    #####
    
    if(is.null(gwasm$W)){
      #cat()
      stop("A GWAS model is needed to create the design matrix for bagging",call.=FALSE)
    }else{ ######## IF MODELS WAS A GWAS MODEL ###########
      
      if(pick){ ############# 'PICK=TRUE' ARGUMENT #################
        cat("Please move the cursor to the GWAS plot, \nclick over the marker dots you want in the design matrix \nand press 'Esc' key when done")
        plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20)
        
        picked <- locator()$x
        marker <- round(picked)
        res <- big.peaks.col(gwasm$W.scores$additive, tre=threshold)
        marker2 <- numeric(length(marker))
        for(t in 1:length(marker)){
          ## in the neighbour of the peak selected
          ff <- which(res$pos %in% c(c(marker[t]-50):c(marker[t]+50)))
          ##
          marker2[t] <- res$pos[ff[which(res$hei[ff] == max(res$hei[ff]))[1]]]
        }
        
        abline(v=marker2, lty=3, col="red")
        legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
        
        XX <- as.matrix(gwasm$W[,marker])
        X1 <- as.data.frame(apply(XX,2, as.factor))
        xvar <- colnames(X1)
        X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
        
      }else{ ############# 'PICK=FALSE' ARGUMENT #################
        big.peaks.col <- function(x, tre){
          r <- rle(x)
          v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
          pos <- v[which(x[v] > tre)] #position of the real peaks
          hei <- x[pos]
          out <- list(pos=pos,hei=hei)
          return(out)
        }
        make.full <- function(X) {
          svd.X <- svd(X)
          r <- max(which(svd.X$d > 1e-08))
          return(as.matrix(svd.X$u[, 1:r]))
        }
        ########################
        # plot(transfft(gwasm$W.scores$additive,.02), type="l")
        # abline(v=res$pos, lty=3,col="red")
        #res <- big.peaks.col(transfft(gwasm$W.scores$additive,.02),0) # smoothed
        #hh <- which(res$hei %in% sort(res$hei, decreasing=TRUE)[1:nmar])
        #res <- list(pos=res$pos[hh], hei=res$hei[hh])
        
        condicion <- gwasm$map
        if(!is.null(condicion)){ # if map is present
          res <- big.peaks.col(as.vector(gwasm$map$p.val), tre=threshold)
        }else{
          res <- big.peaks.col(as.vector(gwasm$W.scores$additive), tre=threshold)
        }
        
        
        if(length(res$pos) >= nmar){ # if enough markers significant
          ## build the number of cluster plus 5 to make sure you have repeatability
          if(length(res$pos) <= nmar+15){
            if(length(res$pos) < nmar){
              res$clus <- kmeans(res$pos,length(res$pos))$cluster 
            }else{
              res$clus <- kmeans(res$pos,nmar-1)$cluster 
            }
            
          }else{
            res$clus <- kmeans(res$pos,nmar+15)$cluster 
          }
          
          heights <- numeric()
          for(i in 1:max(res$clus)){
            vv <- which(res$clus == i)
            heights[i] <- (res$hei[vv][which(res$hei[vv] == max(res$hei[vv]))])[1]
          }
          ## the top clusters (tallest peaks)
          maxo <- which(heights %in% sort(heights, decreasing = TRUE)[1:nmar])
          ## subset the 'res' object
          good <- which(res$clus %in% maxo)
          res2 <- list(pos=res$pos[good], hei=res$hei[good], clus=res$clus[good])
          
          marker <- numeric()
          for(i in unique(res2$clus)){
            vv <- which(res2$clus == i)
            marker[i] <- (res2$pos[vv][which(res2$hei[vv] == max(res2$hei[vv]))])[1]
          }
          marker <- as.numeric(na.omit(marker))
          #marker <- res$pos
          #########
          #########
          
          if(!is.null(condicion)){ #if map was present
            if(plotting){
              col.scheme <- rep((transp(c("cadetblue","red"))),30)
              plot(gwasm$map$p.val, bty="n", col=col.scheme[factor(gwasm$map$Chrom, levels = unique(gwasm$map$Chrom, na.rm=TRUE))], xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), pch=20, cex=2, las=2)
              init.mrks <- apply(data.frame(unique(gwasm$map$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=gwasm$map$Chrom)
              fin.mrks <- apply(data.frame(unique(gwasm$map$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=gwasm$map$Chrom)
              inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
              axis(side=1, at=inter.mrks, labels=paste("Chr",unique(gwasm$map$Chrom),sep=""), cex.axis=.5)
              abline(v=marker, lty=3, col="red")
              legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
            }
            marker <- as.character(gwasm$map$Locus[marker])
          }else{ # if map was not present
            if(plotting){
              plot(gwasm$W.scores$additive, col=transp("cadetblue"), pch=20, cex=1.3)
              abline(v=marker, lty=3, col="red")
              legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
            }
          }
          
          
          
          if(method=="cluster"){
            XX <- as.matrix(gwasm$W[,marker])
          }
          #########################
          if(method == "maximum"){
            n.mark=nmar
            pval <- gwasm$W.scores$additive
            names(pval) <- colnames(gwasm$W)
            marker <- order(pval)[1:n.mark]  
            XX <- as.matrix(gwasm$W[,marker])
          }
          
          X1 <- as.data.frame(apply(XX,2, as.factor))
          
          dada <- data.frame(y=(gwasm$fitted.y[,1]), X1)
          fit <- lm(as.formula(paste("y~ ", paste(c(colnames(dada)[-1]), collapse="+"))), data=dada)
          step <- MASS::stepAIC(fit,direction="both",trace=FALSE)
          # good markers after stepwise
          xvar <-(as.character(attr(summary(step)$terms, "variables")))[-c(1:2)]
          #head(X1)
          #xvar <- colnames(X1)
          #cat(("Markers selected"))
          #cat(paste(xvar))
          X2 <- model.matrix(as.formula(paste("~ ", paste(c("1", xvar), collapse="+"))), data=X1)
          #X2 <- make.full(X2)
          #X2 <- as.data.frame(as.matrix(X2))
          #fit <- lm(gwasm$fitted.y~ -1 + as.matrix(X2))
          #fit <- lm(make.formula(colnames(X1)),data=data.frame(X1,y=gwasm$fitted.y)) #lm using 100 marks
          
          #step <- stepAIC(fit,direction="both",trace=FALSE)  #forward-backward stepwise regression
          
        }else{ # if not enough markers
          #cat()
          stop("Not enough significant markers in your model to create a design matrix \nwith the number of markers specified by you. Please lower the 'threshold' \nargument or the number of markers in the 'nmar' argument",call. = FALSE)
        } #### end of if enough markers
        
      } ############# END OF 'PICK' ARGUMENT #################
    }
    if(only.mark){
      X2 <- colnames(X1)
    }
    return(X2) #X2
  }
  
  #univariate model
  if(names(gwasm)[1] == "var.comp"){ 
    XXX <- bagi(gwasm=gwasm, nmar=nmar, threshold=threshold, pick=pick, method=method, only.mark=only.mark, plotting=plotting)
  }else{ # univariate in parallel
    XXX <- lapply(gwasm, bagi, nmar=nmar, threshold=threshold, pick=pick, method=method, only.mark=only.mark, plotting=plotting)
  }
  
  return(XXX)
}

maxi.qtl <- function(sma, distan=10, no.qtl=5, q95=2.5, LOD.int=FALSE, LODdrop=2, allCI=TRUE){
  param=NULL
  sma <- data.frame(sma)
  rnn <- rownames(sma)
  param <- vector("list", length(unique(sma$chr)))
  param <- lapply(param, function(x){x <- c(q95,no.qtl)})
  ## sma is the dataframe
  ## param is a list with 2 values, lod threshold and #of qtls
  sma <- apply(sma, 2, FUN=as.numeric)
  sma <- as.data.frame(sma)
  rownames(sma) <- rnn
  if(is.null(param)){
    param <- vector("list", length(unique(sma$chr)))
    param <- lapply(param, function(x){x <- c(1,1)})
  }
  
  big.peaks.col <- function(x, tre){
    r <- rle(x)
    v <- which(rep(x = diff(sign(diff(c(-Inf, r$values, -Inf)))) == -2, times = r$lengths))
    pos <- v[which(x[v] > tre)] #position of the real peaks
    hei <- x[pos]
    out <- list(pos=pos,hei=hei)
    return(out)
  }
  
  result <- list(NA)
  for(j in 1:max(sma$chr,na.rm=TRUE)){
    prov <- sma[which(sma$chr == j),]
    para <- param[[j]]
    roxy <- big.peaks.col(prov$lod, tre=para[1]) # first element is the lod threshold, 2nd is the #of qtls
    if(length(roxy$pos) > 0){
      prov2.1 <- prov[roxy$pos,] # to keep
      prov2 <- prov[roxy$pos,] # to discard
      
      w <- numeric()
      res <- data.frame(matrix(NA,nrow=para[2], ncol=3)) # 3 columns because always rqtl gives that
      names(res) <- names(prov2)
      
      GOOD <- TRUE
      k=1
      while(GOOD){
        #for(k in 1:para[2]){# for the peaks required
        mm <- max(prov2$lod, na.rm=TRUE)
        w[k] <- which(prov2$lod == mm)[1]
        provi <- prov2$pos[w[k]]
        
        res[k,c(1:3)] <- prov2[w[k],c(1:3)]
        rownames(res)[k] <- rownames(prov2)[w[k]] #$$$$$$$$$$$$$$$$$$$$
        to.elim <- distan#abs(prov2$pos[w[k]] - 10)
        zz <- setdiff(which(prov2$pos > (provi - to.elim) & prov2$pos < (provi + to.elim) ), w[k])# [-w[k]]
        if(length(zz) > 0){
          prov2 <- prov2[-zz,]
        }else{prov2 <- prov2}
        selected <- which(prov2$lod == mm)[1]
        prov2 <- prov2[-selected,]
        #rownames(prov2.1)[-selected]
        k=k+1
        if(length(which(is.na(prov2))) > 0 | dim(prov2)[1] == 0 | k > no.qtl){GOOD<-FALSE}
        #}
      }
      
      #maxqtls <- sort(prov2$lod, decreasing = TRUE)
      result[[j]] <- res#prov2[which(prov2$lod %in% maxqtls[1:para[2]]),]
    }
    
  }# end for each chromosome
  #############
  result <- lapply(result, function(x){unique(x)})
  good <- which(unlist(lapply(result, is.null)) == FALSE)
  result <- result[good]
  
  good <- which(unlist(lapply(result, function(m){is.null(dim(m))})) == FALSE)
  result <- result[good]
  
  res2 <- do.call("rbind", result)
  
  res2<- res2[which(!is.na(res2$chr)),]
  
  res3 <- res2
  if(LOD.int){
    
    #print(res3)
    
    lod2 <-list() 
    for(i in 1:dim(res3)[1]){
      #apply(res3,1,function(x,sma){
      
      x <- res3[i,]
      #print(x)
      mpp <- sma[which(sma$chr == as.numeric(x[1])),]
      
      babo <- which(rownames(mpp) == rownames(x))#which(mpp$chr == as.numeric(x[1]) & mpp$pos == as.numeric(x[2]))
      toch <- mpp[babo,]
      baba <- babo-1
      babu <- babo+1
      rt=0
      #print(toch)
      while((rt < LODdrop) & (baba > 1)){ # 2-lod interval
        rt <- abs(as.numeric(mpp[baba,] - toch)[3])
        baba <- baba - 1
      }
      ## baba bow tell us what is the lower bound
      mpp[baba,]
      
      rt=0
      while((rt < LODdrop) & (babu < dim(mpp)[1])){
        rt <- abs(as.numeric(mpp[babu,] - toch)[3])
        babu <- babu + 1
      }
      mpp[babu,]
      if(allCI){
        resx <- rbind( mpp[baba:(babo-1),], toch,  mpp[(babo+1):babu,])
      }else{
        resx <- rbind( mpp[baba,], toch,  mpp[babu,])
      }
      
      
      lod2[[i]] <- resx
    }
    
    res2 <- do.call("rbind",lod2)
  }
  return(res2)
  
}

fdr <- function(p, fdr.level=0.05){
  ##### transform to real p-values
  # if maximum value is grater than 1 means that is in -log10 or LOD scale
  # if maximum value is less than one means that the user is using already raw  p.values
  if(max(p, na.rm = TRUE) > 1){ # is in -lod 10 scale
    pval <- 10^-p
  }else{
    pval <- p
  }
  ########## make sure there is a value
  #ro1 <- c(pval)
  #ro2 <- p.adjust(ro1, method="fdr")
  #ro3 <- -log(c(ro2,0.05), 10)
  #ro4 <- 10^-ro3
  #ro5 <-p.adjust(ro4, method="fdr")
  ##### adjust for FDR ---- ADJUSTED IN P.VAL SCALE -----
  pvalA <- p.adjust(pval, method="fdr")
  #plot(pvalA)
  ##### ---- VALS IN LOG.10 SCALE -----
  pvalA.l10 <- -log(pvalA, 10)
  #plot(pvalA.l10)
  ##### ---- FDR IN P.VAL SCALE FOR ADJUSTED ----
  fdr.p <- fdr.level
  ## FDR in original scale
  #1) transform the values to p-values
  # pvalA
  #2) find which value adjusted is equal to 0.05 and go back to the original value
  sortedd <- sort(pvalA, decreasing = TRUE)
  closer <- sortedd[which(sortedd < fdr.level)[1]] # closest value found to the fdr.level indicated by the user
  vv <- which(pvalA == closer)[1]
  
  #vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0){
    fdr.10 <- p[vv]#fdr.or <- min(p[vv])
    #fdr <- 0.05
  }else{
    fdr.10 <- NULL
  }
  ######
  result <- list(p.ad=pvalA, fdr.p=fdr.p, p.log10=pvalA.l10, fdr.10=fdr.10 )
  return(result)
}

fdr2 <- function(p, fdr.level=0.05){
  
  
  ##### transform to real p-values
  # if maximum value is grater than 1 means that is in -log10 or LOD scale
  # if maximum value is less than one means that the user is using already raw  p.values
  if(max(p, na.rm = TRUE) > 1){ # is in -log 10 scale
    pval <- 10^-p
  }else{
    pval <- p
  }
  
  ##for endelmans function
  pen <- -log(pval,10)
  ########## make sure there is a value
  #ro1 <- c(pval)
  #ro2 <- p.adjust(ro1, method="fdr")
  #ro3 <- -log(c(ro2,0.05), 10)
  #ro4 <- 10^-ro3
  #ro5 <-p.adjust(ro4, method="fdr")
  ##### adjust for FDR ---- ADJUSTED IN P.VAL SCALE -----
  pvalA <- p.adjust(pval, method="fdr")
  #plot(pvalA)
  ##### ---- VALS IN LOG.10 SCALE -----
  pvalA.l10 <- -log(pvalA, 10)
  #plot(pvalA.l10)
  ##### ---- FDR IN P.VAL SCALE FOR ADJUSTED ----
  fdr.p <- fdr.level
  ## FDR in original scale
  #1) transform the values to p-values
  # pvalA
  #2) find which value adjusted is equal to 0.05 and go back to the original value
  sortedd <- sort(pvalA, decreasing = TRUE)
  closer <- sortedd[which(sortedd < fdr.level)[1]] # closest value found to the fdr.level indicated by the user
  vv <- which(pvalA == closer)[1]
  
  #vv <- which(pvalA.l10 > fdr.ad)
  if(length(vv)>0 & !is.na(vv)){
    fdr.10 <- p[vv]#fdr.or <- min(p[vv])
    #fdr <- 0.05
  }else{
    
    fdrendel <- function(dd, fdr.level=0.05){
      qvalue <- function(p) {
        smooth.df = 3
        if (min(p) < 0 || max(p) > 1) {
          print("ERROR: p-values not in valid range.")
          return(0)
        }
        lambda = seq(0, 0.9, 0.05)
        m <- length(p)
        pi0 <- rep(0, length(lambda))
        for (i in 1:length(lambda)) {
          pi0[i] <- mean(p >= lambda[i])/(1 - lambda[i])
        }
        spi0 <- smooth.spline(lambda, pi0, df = smooth.df)
        pi0 <- predict(spi0, x = max(lambda))$y
        pi0 <- min(pi0, 1)
        #print(pi0)
        if (pi0 <= 0) {
          #pi0 <- abs(pi0)
          print("ERROR: The estimated pi0 <= 0. Check that you have valid p-values.")
          return(0)
        }
        u <- order(p)
        qvalue.rank <- function(x) {
          idx <- sort.list(x)
          fc <- factor(x)
          nl <- length(levels(fc))
          bin <- as.integer(fc)
          tbl <- tabulate(bin)
          cs <- cumsum(tbl)
          tbl <- rep(cs, tbl)
          tbl[idx] <- tbl
          return(tbl)
        }
        v <- qvalue.rank(p)
        qvalue <- pi0 * m * p/v
        qvalue[u[m]] <- min(qvalue[u[m]], 1)
        for (i in (m - 1):1) {
          qvalue[u[i]] <- min(qvalue[u[i]], qvalue[u[i + 1]],
                              1)
        }
        return(qvalue)
      }
      ## go back to p values
      q.ans <- qvalue(10^-dd)
      temp <- cbind(q.ans, dd)
      temp <- temp[order(temp[, 1]), ]
      
      temp2 <- tapply(temp[, 2], temp[, 1], mean)
      qvals <- as.numeric(rownames(temp2))
      x <- which.min(abs(qvals - fdr.level))
      first <- max(1, x - 2)
      last <- min(x + 2, length(qvals))
      if ((last - first) < 4) {
        last <- first + 3
      }
      #print(qvals[first:last])
      #print(temp2[first:last])
      splin <- smooth.spline(x = qvals[first:last], y = temp2[first:last],
                             df = 3)
      popo <- predict(splin, x = fdr.level)$y
      return(popo)
    }
    
    fdr.10 <- fdrendel( pen,fdr.level = fdr.level)
  }
  ######
  result <- list(p.ad=pvalA, fdr.p=fdr.p, p.log10=pvalA.l10, fdr.10=fdr.10 )
  return(result)
}
