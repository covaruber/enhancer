A.matr <- function (X, min.MAF = NULL) {
  
  X <- as.matrix(X)
  n <- nrow(X)
  frac.missing <- apply(X, 2, function(x) {
    length(which(is.na(x)))/n
  })
  missing <- max(frac.missing) > 0
  freq <- apply(X + 1, 2, function(x) {
    mean(x, na.rm = missing)
  })/2
  MAF <- apply(rbind(freq, 1 - freq), 2, min)
  if (is.null(min.MAF)) {
    min.MAF <- 1/(2 * n)
  }
  max.missing <- 1 - 1/(2 * n)
  markers <- which((MAF >= min.MAF) & (frac.missing <= max.missing))
  m <- length(markers)
  var.A <- 2 * mean(freq[markers] * (1 - freq[markers]))
  one <- matrix(1, n, 1)
  mono <- which(freq * (1 - freq) == 0)
  X[, mono] <- 2 * tcrossprod(one, matrix(freq[mono], length(mono), 
                                          1)) - 1
  freq.mat <- tcrossprod(one, matrix(freq[markers], m, 1))
  W <- X[, markers] + 1 - 2 * freq.mat
  A <- tcrossprod(W)/var.A/m
  return(A)
}

adiag1 <- function (..., pad = as.integer(0), do.dimnames = TRUE){
  args <- list(...)
  if (length(args) == 1) {
    return(args[[1]])
  }
  if (length(args) > 2) {
    jj <- do.call("Recall", c(args[-1], list(pad = pad)))
    return(do.call("Recall", c(list(args[[1]]), list(jj),
                               list(pad = pad))))
  }
  a <- args[[1]]
  b <- args[[2]]
  if (is.null(b)) {
    return(a)
  }
  if (is.null(dim(a)) & is.null(dim(b))) {
    dim(a) <- rep(1, 2)
    dim(b) <- rep(1, 2)
  }
  if (is.null(dim(a)) & length(a) == 1) {
    dim(a) <- rep(1, length(dim(b)))
  }
  if (is.null(dim(b)) & length(b) == 1) {
    dim(b) <- rep(1, length(dim(a)))
  }
  if (length(dim.a <- dim(a)) != length(dim.b <- dim(b))) {
    stop("a and b must have identical number of dimensions")
  }
  s <- array(pad, dim.a + dim.b)
  s <- do.call("[<-", c(list(s), lapply(dim.a, seq_len), list(a)))
  ind <- lapply(seq(dim.b), function(i) seq_len(dim.b[[i]]) +
                  dim.a[[i]])
  out <- do.call("[<-", c(list(s), ind, list(b)))
  n.a <- dimnames(a)
  n.b <- dimnames(b)
  if (do.dimnames & !is.null(n.a) & !is.null(n.b)) {
    dimnames(out) <- mapply(c, n.a, n.b, SIMPLIFY = FALSE)
    names(dimnames(out)) <- names(n.a)
  }
  return(out)
}

add.diallel.vars <- function(df, par1="Par1", par2="Par2",sep.cross="-"){
  # Dummy variables for selfs, crosses, combinations
  df[,"is.cross"] <- ifelse(df[,par1] == df[,par2], 0, 1)
  df[,"is.self"] <- ifelse(df[,par1] == df[,par2], 1, 0)
  df[,"cross.type"] <- ifelse(as.character(df[,par1]) < as.character(df[,par2]), -1,
                              ifelse(as.character(df[,par1]) == as.character(df[,par2]), 0, 1))
  # Dummy variable for the combinations, ignoring the reciprocals
  df[,"cross.id"]<-factor(ifelse(as.character(df[,par1]) <= as.character(df[,par2]),
                                 paste(df[,par1], df[,par2], sep =sep.cross),
                                 paste(df[,par2], df[,par1], sep =sep.cross)) )
  return(df)
}

atcg1234 <- function(data, ploidy=2, format="ATCG", maf=0, multi=TRUE, silent=FALSE, 
                     by.allele=FALSE, imp=TRUE, ref.alleles=NULL){
  
  impute.mode <- function(x) {
    ix <- which(is.na(x))
    if (length(ix) > 0) {
      x[ix] <- as.integer(names(which.max(table(x))))
    }
    return(x)
  }
  ##### START GBS.TO.BISNP DATA ######
  gbs.to.bisnp <- function(x) {
    y <- rep(NA,length(x))
    y[which(x=="A")] <- "AA"
    y[which(x=="T")] <- "TT"
    y[which(x=="C")] <- "CC"
    y[which(x=="G")] <- "GG"
    y[which(x=="R")] <- "AG"
    y[which(x=="Y")] <- "CT"
    y[which(x=="S")] <- "CG"
    y[which(x=="W")] <- "AT"
    y[which(x=="K")] <- "GT"
    y[which(x=="M")] <- "AC"
    y[which(x=="+")] <- "++"
    y[which(x=="0")] <- "NN"
    y[which(x=="-")] <- "--"
    y[which(x=="N")] <- NA
    return(y)
  }
  ##### END GBS.TO.BISNP DATA ######
  imputeSNP <- function(data){
    #######
    data2 <- apply(data,2,function(x){
      areNA <- which(is.na(x))
      if(length(areNA)>0){
        pos.all <- table(data[,1])
        totake <- names(pos.all)[which(pos.all == max(pos.all))]
        x[areNA] <- totake
      }
      return(x)
    })
    #######
    return(data2)
  }
  #### apply with progress bar ######
  apply_pb <- function(X, MARGIN, FUN, ...){
    env <- environment()
    pb_Total <- sum(dim(X)[MARGIN])
    counter <- 0
    pb <- txtProgressBar(min = 0, max = pb_Total,
                         style = 3)
    
    wrapper <- function(...)
    {
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir= env)
      setTxtProgressBar(get("pb", envir= env),
                        curVal +1)
      FUN(...)
    }
    res <- apply(X, MARGIN, wrapper, ...)
    close(pb)
    res
  }
  ###### zero.one function
  zero.one <- function(da){
    # this function takes a matrix of markers in biallelic format and returns a matrix of
    # presense/absense of alleles
    mar.nam <- colnames(da)#unique(gsub("\\.\\d","", names(da))) # find a dot and a number after the dot
    mat.list <- list(NA) # list of matrices for each marker
    wi=0 # counter
    if(!silent){
      count <- 0
      tot <- length(mar.nam)
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
    }
    for(i in 1:length(mar.nam)){ # for each marker
      wi=wi+1
      if(!silent){
        count <- count + 1
      }
      
      v <- which(colnames(da)==mar.nam[i])#grep(mar.nam[i], colnames(da))
      
      if(length(v)==0){
        qqqqq <- grep(mar.nam[i-1],names(da))
        qqqqq2 <- names(da)[qqqqq[length(qqqqq)] + 1]
        
        stop(paste("Marker",qqqqq2,"has a problem"), call.=FALSE)
      }else if(length(v) == 1){ # for markers with a single column
        prov <- matrix(da[,v])
      }else{prov <- da[,v]}
      ##################################
      alls <- unique(unlist(strsplit(prov,"")))
      alls <- alls[which(!is.na(alls))]
      ninds <- dim(prov)[1]
      fff <- apply(data.frame(alls),1,function(h){
        temp <- numeric(length = ninds)
        temp[grep(h,prov)]<-1
        #make sure is full rank
        
        return(temp)
      })#1 # assigning 1's
      #if(FULL){ # if user want to make sure only get the columns that will ensure full rank
      #  fff <- t(unique(t(fff)))
      #}
      colnames(fff) <- paste(mar.nam[i],alls, sep="/")
      
      mat.list[[i]] <- fff
      if(!silent){
        setTxtProgressBar(pb, (count/tot))### keep filling the progress bar
      }
      
    }
    
    fin.mat <- do.call(cbind,mat.list)
    rownames(fin.mat) <- rownames(da)
    #############
    return(fin.mat)
  }
  
  ## remove all markers or columns that are all missing data
  all.na <- apply(data,2,function(x){length(which(is.na(x)))/length(x)})
  bad.na <- which(all.na==1)
  if(length(bad.na) > 0){
    data <- data[,-bad.na]
  }
  
  if(is.null(ref.alleles)){
    #############################
    if(by.allele){ ####&&&&&&&&&&&&&&&&&&&&&& use zero.one function
      ncolsData <- dim(data)[2]
      ncolsData <- max(ncolsData,round(ncolsData/20))
      user.code <- apply(data[,c(1:ncolsData),drop=FALSE], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
      
      AA <- sum(user.code, na.rm = TRUE)/length(user.code)
      if(AA > .9){ # means user is using single letter
        rnd <- rownames(data)
        data <- apply(data,2,gbs.to.bisnp);#W2[1:5,1:5]
        rownames(data) <- rnd
      }
      M <- zero.one(data)
      
    }else{ ###&&&&&&&&&&&&&&&&&&&&&&&&
      n.g <- apply(data,2,function(x){length(table(x))})
      bad <- which(n.g > 3)
      if(length(bad) == dim(data)[2]){
        stop("Error. All your markers are multiallelic. This function requires at least one bi-allelic marker\n")
      }
      
      # tells you which markers have double letter code, i.e. TT instead of T
      # 1: has only one letter
      # 0: has two letters
      ncolsData <- dim(data)[2]
      ncolsData <- max(ncolsData,round(ncolsData/20))
      user.code <- apply(data[,c(1:ncolsData), drop=FALSE], 2, function(x){q <- which(!is.na(x))[1];ss1 <- substr(x[q], start=1,stop=1);ss2 <- substr(x[q], start=2,stop=2);vv1 <-which(c(ss1,ss2)=="");if(length(vv1)>0){y <-1}else{y <- 0}; return(y)})
      AA <- sum(user.code, na.rm = TRUE)/length(user.code)
      if(AA > .9){
        rrn <- rownames(data)
        
        message("Converting GBS or single-letter code to biallelic code\n")
        if(silent){
          data <- apply(data, 2,gbs.to.bisnp)
        }else{
          data <- apply_pb(data, 2,gbs.to.bisnp) 
        }
        rownames(data) <- rrn
        data <- as.data.frame(data)
      }
      #### apply with progress bar ######
      s1 <- rownames(data)
      s2 <- colnames(data)
      data <- as.data.frame(t(data))
      rownames(data) <- s2
      colnames(data) <- s1
      bases <- c("A", "C", "G", "T","l","m","n","p","h","k","-","+","e","f","g","a","b","c","d")
      ## get reference allele function
      get.ref <- function(x, format) {
        if (format == "numeric") {
          ref.alt <- c(0, 1)
        }
        if (format == "AB") {
          ref.alt <- c("A", "B")
        }
        if (format == "ATCG") {
          y <- paste(na.omit(x), collapse = "")
          ans <- apply(array(bases), 1, function(z, y) {
            length(grep(z, y, fixed = T))
          }, y)
          if (sum(ans) > 2) {
            ref.alt <- (bases[which(ans == 1)])[1:2]
            #stop("Error in genotype matrix: More than 2 alleles")
          }
          if (sum(ans) == 2) {
            ref.alt <- bases[which(ans == 1)]
          }
          if (sum(ans) == 1) {
            ref.alt <- c(bases[which(ans == 1)], NA)
          }
        }
        return(ref.alt)
      }
      
      get.multi <- function(x, format) {
        if (format == "numeric") {
          ref.alt <- c(0, 1)
        }
        if (format == "AB") {
          ref.alt <- c("A", "B")
        }
        if (format == "ATCG") {
          y <- paste(na.omit(x), collapse = "")
          ans <- apply(array(bases), 1, function(z, y) {
            length(grep(z, y, fixed = T))
          }, y)
          if (sum(ans) > 2) {
            ref.alt <- TRUE
          }
          if (sum(ans) == 2) {
            ref.alt <- FALSE
          }
          if (sum(ans) == 1) {
            ref.alt <- FALSE
          }
        }
        return(ref.alt)
      }
      
      ####################################
      ## convert to matrix format
      ####################################
      markers <- as.matrix(data)
      ####################################
      # get reference alleles
      ####################################
      message("Obtaining reference alleles\n")
      if(silent){
        tmp <- apply(markers, 1, get.ref, format=format)
      }else{
        tmp <- apply_pb(markers, 1, get.ref, format=format) 
      }
      
      if(multi){ # if markers with multiple alleles should be removed
        message("Checking for markers with more than 2 alleles. If found will be removed.\n")
        if(silent){
          tmpo <- apply(markers, 1, get.multi, format = format)
        }else{
          tmpo <- apply_pb(markers, 1, get.multi, format = format) 
        }
        ###&&&&&&&&&&&& HERE WE MUST INSERT THE NEW FUNCTIONALITY, WHERE WE DETECTED MULTIPLE ALLELES
        multi.allelic <- which(!tmpo) # good markers
        markers <- markers[multi.allelic,,drop=FALSE]
        tmp <- tmp[, multi.allelic,drop=FALSE]
      }
      
      Ref <- tmp[1, ]
      Alt <- tmp[2, ]
      ####################################
      ## bind reference allele and markers and convert to numeric format based on the 
      # reference/alternate allele found
      ####################################
      message("Converting to numeric format\n")
      if(silent){
        M <- apply(cbind(Ref, markers), 1, function(x) {
          y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
          ans <- as.integer(lapply(y, function(z) {
            ifelse(z[1] < 0, ploidy, ploidy - length(z))
          }))
          return(ans)
        })
      }else{
        M <- apply_pb(cbind(Ref, markers), 1, function(x) {
          y <- gregexpr(pattern = x[1], text = x[-1], fixed = T)
          ans <- as.integer(lapply(y, function(z) {
            ifelse(z[1] < 0, ploidy, ploidy - length(z))
          }))
          return(ans)
        })
      }
      gid.geno <- s1 #colnames(geno)
      rownames(M) <- gid.geno
      ####################################
      # identify bad markers
      ####################################
      bad <- length(which(!is.element(na.omit(M), 0:ploidy)))
      if (bad > 0) {
        stop("Invalid marker calls.")
      }
      
    }
    #rownames(M) <- rownames(data)
    ####################################
    rownames(tmp) <- c("Alt","Ref")
  }else{# user provides reference alleles and just want a conversion
    
    common.mark <- intersect(colnames(data), colnames(ref.alleles))
    data <- data[,common.mark, drop=FALSE]
    tmp <- ref.alleles[,common.mark, drop=FALSE]; #rownames(refa) <- c("Alt","Ref")
    message("Converting to numeric format\n")
    M <- apply_pb(data.frame(1:ncol(data)),1,function(k){
      x <- as.character(data[,k])
      x2 <- strsplit(x,"")
      x3 <- unlist(lapply(x2,function(y){length(which(y == tmp[2,k]))}))
      return(x3)
    })
    #M <- M-1
    colnames(M) <- colnames(data)
    
  }
  
  ####################################
  # by column or markers calculate MAF
  ####################################
  message("Calculating minor allele frequency (MAF)\n")
  if(silent){
    MAF <- apply(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }else{
    MAF <- apply_pb(M, 2, function(x) {
      AF <- mean(x, na.rm = T)/ploidy
      MAF <- ifelse(AF > 0.5, 1 - AF, AF)
    })
  }
  ####################################
  # which markers have MAF > 0, JUST GET THOSE
  ####################################
  polymorphic <- which(MAF > maf)
  M <- M[, polymorphic, drop=FALSE]
  ####################################
  # function to impute markers with the mode
  ####################################
  
  # time to impute
  if(imp){
    missing <- which(is.na(M))
    if (length(missing) > 0) {
      message("Imputing missing data with mode \n")
      if(silent){
        M <- apply(M, 2, impute.mode)
      }else{
        M <- apply_pb(M, 2, impute.mode)
      }
    }
  }else{
    message("Imputation not required. Be careful using non-imputed matrices in mixed model solvers\n")
  }
  ## ploidy 2 needs to be adjusted to -1,0,1
  # if(ploidy == 2){
  #   M <- M - 1
  # }
  
  return(list(M=M,ref.alleles=tmp))
}

atcg1234BackTransform <- function(marks, refs){
  marks2 <- matrix(NA, nrow=nrow(marks), ncol = ncol(marks))
  ploidy <- diff(range(marks, na.rm = TRUE))
  # center <- ploidy #/ 2
  for(iMark in 1:ncol(marks)){ # iMark=1
    
    marks2[,iMark] <- apply(as.data.frame(marks[,iMark]),1,function(x){
      if(is.na(x)){
        NA
      }else{
        gsub(pattern=" ",replacement="",
             paste(c( rep(refs["Alt",colnames(marks)[iMark]], abs(x-ploidy) ),
                      rep(refs["Ref",colnames(marks)[iMark]], x) ), collapse = ""
             )
        )
      }
    })
  }
  rownames(marks2) <- rownames(marks)
  colnames(marks2) <- colnames(marks)
  return(marks2)
}

bathy.colors <- function(n, alpha=1){
  return(rgb(seq(0.9,0,len=n), seq(0.9,0,len=n), 1, alpha))
}

bbasis <- function (x, xl, xr, ndx, deg) 
{
  tpower <- function(x, t, p) {
    (x - t)^p * (x > t)
  }
  dx <- (xr - xl)/ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1)/(gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  B
}

build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE, separator=":", n.batch=1000, verbose=TRUE){
  # build hybrid marker matrix
  
  if(!is.null(custom.hyb)){
    pheno <- custom.hyb
    found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
    if(found != 3){
      stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
    }
    return.combos.only=FALSE
  }else{
    a <- rownames(M1)
    b <- rownames(M2)
    pheno <- expand.grid(a,b)
    pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
    pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=separator)
  }
  
  if(!return.combos.only){
    # check that marker matrices are in -1,0,1 format
    checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
    checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))
    
    checkM1[which(checkM1 > 0)] <- 1
    checkM2[which(checkM2 > 0)] <- 1
    
    if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
    }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
      message("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
      message("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      
    }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
      message("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
    }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
      message("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
    }
    
    n.batch <- min(c(n.batch,nrow(pheno)))
    
    if(nrow(pheno)>0){ # if there is hybrids to build
      ## build the marker matrix for batches of n.batch hybrids
      batches <- sort(rep(1:1000,min(c(nrow(pheno),n.batch))))
      pheno$batch <- batches[1:nrow(pheno)]
      data.usedBatches <- split(pheno, pheno$batch)
      
      M1 <- as(M1, "CsparseMatrix")
      M2 <- as(M2, "CsparseMatrix")
      # start the loop
      for(i in 1:length(data.usedBatches)){
        
        ## add markers coming from parents M1
        Z1 <- Matrix::sparse.model.matrix(~Var1-1,data.usedBatches[[i]]);dim(Z1); 
        colnames(Z1) <- gsub("Var1","",colnames(Z1))
        M1r <- M1[colnames(Z1),,drop=FALSE]
        
        ## add markers coming from parents M2
        Z2 <- Matrix::sparse.model.matrix(~Var2-1,data.usedBatches[[i]]);dim(Z2); 
        colnames(Z2) <- gsub("Var2","",colnames(Z2))
        M2r <- M2[colnames(Z2),, drop=FALSE]
        
        
        hyb.names <- data.usedBatches[[i]]$hybrid
        ## marker matrix for hybrids one for each parent
        if(verbose){
          message(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))
          message("Extracting M1 contribution\n")
        }
        
        if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
          Md <- Z1 %*% M1r;  # was already converted to -1,1
        }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
          Md <- 2*Z1 %*% M1r - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
        }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
          Md <- Z1 %*% M1r - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
        }
        
        if(verbose){message("Extracting M2 contribution\n")}
        if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
          Mf <- Z2 %*% M2r;  # was already converted to -1,1
        }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
          Mf <- 2*Z2 %*% M2r - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
        }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
          Mf <- Z2 %*% M2r - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
        }
        
        ## marker matrix coded as additive -1,0,1
        Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
        rownames(Mdf) <- hyb.names
        #hist(Mdf)
        
        ## dominance matrix for hybrids (0,1 coded)
        Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
        rownames(Delta) <- hyb.names
        
        if(i == 1){
          HMM.add <- Mdf
          HMM.dom <- Delta
        }else{
          HMM.add <- rbind(Mdf, HMM.add)
          HMM.dom <- rbind(Delta, HMM.dom)
        }
      }
    }else{
      HMM.add <- HMM.dom <- Matrix::Matrix(NA, nrow=0, ncol=ncol(M1)); colnames(M) <- colnames(M1)
    }
    
    #hist(Delta)
    if(verbose){message("Done!!\n")}
    return(list(HMM.add=HMM.add, HMM.dom=HMM.dom, data.used=pheno))
    
  }else{
    return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
  }
}

imputev <- function (x, method = "median", by=NULL) {
  if (is.numeric(x)) {
    if(is.null(by)){
      by <- rep("A",length(x))
    }
    ms <- aggregate(x~by, FUN=method, na.rm=TRUE)
    rownames(ms) <- ms$by
    y <- ms[by,2]
    x[which(is.na(x))] <- y[which(is.na(x))]
  } else { # if factor
    tt <- table(x)
    x[which(is.na(x))] <- names(tt)[which(tt == max(tt))]
  }
  return(x)
}

I.matr <- function(x){
  labels0 <- sort(unique(x))
  II <- Matrix::Diagonal(n=length(labels0))
  colnames(II) <- rownames(II) <- labels0
  return(II)
}

jet.colors <- function(n, alpha=1) {
  if(n > 0) {
    if(length(alpha) != 1 & length(alpha) != n) {
      message('Warning: using only first alpha value')
      alpha <- alpha[1]
    }
    if(length(alpha) == 1) {
      alpha <- rep(alpha, n)
    }
    ## TODO Include alpha values
    return(colorRampPalette(c('#000066', 'blue', 'cyan', 'yellow',
                              'red', '#660000'))(n))
  } else {
    ## Return an empty character string if they requested nothing.
    character()
  }
}

LD.decay <- function(markers,map,silent=FALSE,unlinked=FALSE,gamma=.95){
  
  #if(is.null(rownames(map))){
  good <- which(!duplicated(map$Locus))
  map <- map[good,]
  rownames(map) <- map$Locus 
  #}
  ## clean markers and map from sd=0 and markers with no position
  markers <- markers[,which(apply(markers,2,sd)>0)]
  map <- map[which(!is.na(map$LG)),] # which(apply(map,1,function(x){length(which(is.na(x)))}) ==0)
  #############
  fullmap <- map
  lgs <- unique(map$LG)
  ############################
  dat.list <- list(NA) # will contain the LD (r) and genetic distance (d) for each LG
  ############################
  if (!silent) {
    count <- 0
    tot <- length(lgs)
    pb <- txtProgressBar(style = 3)
    setTxtProgressBar(pb, 0)
  }
  
  if(!unlinked){
    
    LDM <- list()
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      #ords2 <- ords
      head(ords2);tail(ords2)
      #################
      #intersect markers and map
      #################  
      inn <- intersect(ords2$Locus,colnames(markers))
      
      if(length(inn)>0){
        ords2 <- ords2[inn,]
        
        cor2.mat <- cor(as.matrix(markers[,inn]),use="pairwise.complete.obs")^2
        LDM[[k]] <- cor2.mat
        
        N <- dim(markers[,inn])[1]
        cor2.matp <- 1-pchisq(cor2.mat*N,df=1)
        
        #cor(fofo[,which(apply(fofo,2,sd) > 0)])^2
        # trnsforming to double haploid data
        # this is a mtrix of genetic distances
        mat.dist <- matrix(0,dim(ords2)[1], dim(ords2)[1])
        for(i in 1:dim(ords2)[1]){
          for(j in 1:dim(ords2)[1]){
            mat.dist[i,j] <- ords2$Position[i] - ords2$Position[j]
          }
        }
        # get absolute values
        mat.dist <- abs(mat.dist) 
        # make zero the valuies of upper and lower triangular
        cor2.mat[lower.tri(cor2.mat)] <- 0
        cor2.matp[lower.tri(cor2.mat)] <- 0
        mat.dist[lower.tri(mat.dist)] <- 0
        
        y.dist <- as.vector(mat.dist) #los deshace columna por column
        y.cor <- as.vector(cor2.mat) 
        p.cor <- as.vector(cor2.matp)
        
        dat <- data.frame(d=y.dist,r2=y.cor,p=p.cor)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
        if (!silent) {
          setTxtProgressBar(pb, (count/tot))
        }
      }else{
        dat <- data.frame(d=0,r2=0,p=0)
        head(dat)
        dat2 <- dat[-which(dat$d == 0 & dat$r == 0),]
        dat.list[[k]] <- dat2
      }
      
    }
    # make a big data frame
    if (!silent) {
      setTxtProgressBar(pb, (count/tot))
    }
    
    big <- do.call("rbind",dat.list)
    
    # big contains all the LD and and genetic distances for all groups
    resp <- list(by.LG=dat.list, all.LG=big, LDM=LDM)
    
  }else{ # if WANTS TO KNOW THE THRESHOLD OF UNLINKED
    doo <- intersect(colnames(markers), fullmap$Locus)
    cor2.mat <- cor(as.matrix(markers[,doo]),use="pairwise.complete.obs")^2
    
    dat.list <- numeric()#list()#store the r2 values of unlinked markers with chromosome k
    for(k in lgs){ # for each linkaghe group
      if (!silent) {
        count <- count + 1
      }
      # markers in kth chromosome
      ords2 <- fullmap[which(fullmap[,"LG"] == k),]
      # markers not in kth chromosome
      ords3 <- fullmap[which(fullmap[,"LG"] != k),]
      # markers in kth and present
      sik <- intersect(ords2$Locus,colnames(cor2.mat))
      # markers not in kth and present
      nok <- intersect(ords3$Locus,colnames(cor2.mat))
      #ords2 <- ords
      #head(ords2);tail(ords2)
      #into <- setdiff(nok,ords2$Locus)
      
      step1 <- cor2.mat[sik,nok]
      #boxplot(unlist(step1))
      dat.list[k] <- quantile(sqrt(step1[upper.tri(step1)]),gamma)^2
      
      if (!silent) {
        setTxtProgressBar(pb, (count/tot))
      }
    }
    # make a big data frame
    resp <- list(by.LG=dat.list, all.LG=mean(dat.list))
  }
  return(resp)
}

leg <- function(x,n=1,u=-1,v=1, intercept=TRUE, intercept1=FALSE){
  
  init0 <- as.character(substitute(list(x)))[-1L]
  
  if(system.file(package = "orthopolynom") == ""){
    stop("Please install the orthopolynom package to use this function.",call. = FALSE)
  }
  requireNamespace("orthopolynom",quietly=TRUE)
  (leg4coef <- orthopolynom::legendre.polynomials(n=n, normalized=TRUE))
  leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                  x=orthopolynom::scaleX(x, u=u, v=v))))
  colnames(leg4) <- paste("leg",0:(ncol(leg4)-1),sep="")
  if(!intercept){
    leg4 <- leg4[, 2:ncol(leg4), drop = FALSE]
  }
  if(intercept1){
    leg4 <- leg4*sqrt(2)
    # leg4[,1] <- leg4[,1]*sqrt(2)
  }
  attr(leg4,"variables") <- c(init0)
  return(leg4)
}

logspace <- function (x, p = 2) {
  if(var(x, na.rm=TRUE)>0){
    D = max(x, na.rm=TRUE)
    C = min(x, na.rm=TRUE)
    mysigns <- sign(x)
    y = abs(x)^(1/p)
    y <- y * mysigns
    B = max(y, na.rm=TRUE)
    A = min(y, na.rm=TRUE)
    scale = (D - C)/(B - A)
    offset = -A * (D - C)/(B - A) + C
    return(y * scale + offset)
  }else{
    return(x)
  }
}


map.plot <- function(data, trait=NULL, trait.scale="same", col.chr=NULL, col.trait=NULL, type="hist", cex=0.4, lwd=1, cex.axis=0.4, cex.trait=0.8, jump=5 ){
  sasa <- which(colnames(data) == "Locus")
  if(length(sasa) > 0){
    data$Locus <- as.character(data$Locus) 
  }
  ## transparent function
  ## data needs to have 2 columns; LG and Position
  ## trait needs to indicate the name to plot in the chromosome
  ## the trait can be expressed as "dot", "line" or "hist"
  ## the trait scale can be "same" or "ind", which is same for all or individual
  ## cex is only the cex for the ruler of cM
  ## cex.axis is for the titles
  transp <- function (col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
  }
  #####
  draw.arc <- function (x = 1, y = NULL, radius = 1, angle1 = deg1 * pi/180, 
                        angle2 = deg2 * pi/180, deg1 = 0, deg2 = 45, n = 0.05, col = NA, 
                        lwd = NA, ...) {
    getYmult<-function () {
      if (dev.cur() == 1) {
        warning("No graphics device open.")
        ymult <- 1
      }
      else {
        xyasp <- par("pin")
        xycr <- diff(par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
      }
      return(ymult)
    }
    if (all(is.na(col))) 
      col <- par("col")
    if (all(is.na(lwd))) 
      lwd <- par("lwd")
    xylim <- par("usr")
    ymult <- getYmult()
    devunits <- dev.size("px")
    draw.arc.0 <- function(x, y, radius, angle1, angle2, n, col, 
                           lwd, ...) {
      delta.angle <- (angle2 - angle1)
      if (n != as.integer(n)) 
        n <- as.integer(1 + delta.angle/n)
      delta.angle <- delta.angle/n
      angleS <- angle1 + seq(0, length = n) * delta.angle
      angleE <- c(angleS[-1], angle2)
      if (n > 1) {
        half.lwd.user <- (lwd/2) * (xylim[2] - xylim[1])/devunits[1]
        adj.angle = delta.angle * half.lwd.user/(2 * (radius + 
                                                        half.lwd.user))
        angleS[2:n] = angleS[2:n] - adj.angle
        angleE[1:(n - 1)] = angleE[1:(n - 1)] + adj.angle
      }
      p1x <- x + radius * cos(angleS)
      p1y <- y + radius * sin(angleS) * ymult
      p2x <- x + radius * cos(angleE)
      p2y <- y + radius * sin(angleE) * ymult
      segments(p1x, p1y, p2x, p2y, col = col, lwd = lwd, ...)
    }
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    a1 <- pmin(angle1, angle2)
    a2 <- pmax(angle1, angle2)
    angle1 <- a1
    angle2 <- a2
    args <- data.frame(x, y, radius, angle1, angle2, n, col, 
                       lwd, stringsAsFactors = FALSE)
    for (i in 1:nrow(args)) do.call("draw.arc.0", c(args[i, ], 
                                                    ...))
    invisible(args)
  }
  ## this function takes a dataframe with 2 basic columns:
  ## LG containint the linkage group
  ## Position, the position in cM
  len <- numeric()
  for(i in 1:max(unique(data$LG))){
    len[i] <- max(data[which(data$LG == i),"Position"])
  }
  coree <- which(len == max(len))
  coree1 <- data[which(data$LG == coree),]
  linesss <- coree1$Position / max(coree1$Position)
  # if one trait wants to be plotted or not
  if(!is.null(trait)){
    cols <- 1#(dim(data)[2] - 2)
  }else{cols <- 0}
  extra <- length(len) + cols *length(len)
  fact <- 1/ extra
  fact2 <- fact + (fact*cols) # real separation between chromosomes
  ## colors for traits
  if(!is.null(col.trait)){
    col.trait <- col.trait
  }else{col.trait <- c(1:6,1:6)}
  # colors for chromosomes
  
  
  plot.new()
  for(j in 1:max(unique(data$LG))){
    
    prov <- data[which(data$LG == j),] # extract the jth LG
    ## add zeros to the beggining and end of the LG so the curves look good
    dddd <- prov[1,]
    dddd2 <- prov[dim(prov)[1],]
    dddd[1,which(names(dddd) != "LG")] <- 0
    dddd2[1,which(names(dddd) != "LG" & names(dddd) != "Position")] <- 0
    prov <- rbind(dddd,prov,dddd2)
    ##-----------------------------------------------------------
    ## CHROMOSOME CHUNK
    chr <- prov$Position / max(coree1$Position)
    ### -------------------------------------------------
    ### depending if is a genetic map or physical map we decide how the ruler will be
    #if(max(prov$Position) < 1000){mark <- 10}
    #if(max(prov$Position) > 1000 & max(prov$Position) < 10000){mark <- 100}
    #if(max(prov$Position) > 10000 & max(prov$Position) < 100000){mark <- 1000}
    #if(max(prov$Position) > 100000 & max(prov$Position) < 1000000){mark <- 10000}
    #if(max(prov$Position) < 1000000){mark <- 100000}
    
    ruler <- 1 - (c(seq(0,max(prov$Position), by=jump), round(max(prov$Position),0 )) / max(coree1$Position) ); ruler2 <- c(seq(0,max(prov$Position), by=jump),round(max(prov$Position),0 ))
    if(!is.null(trait)){
      sss <- (fact2*j)-fact# + j*cols 
    }else{sss <- (fact2*j)}# + j*cols}
    
    
    ## heatmap fr density
    dd2 <- density(chr, n=length(chr))$y # regular density
    dd <- sort(density(chr, n=length(chr))$y, decreasing=T)
    ##
    if(!is.null(col.chr)){
      hc <- colorRampPalette(c(col.chr[1], col.chr[2]))( length(dd) )
    }else{hc <- gray.colors(n=length(dd), start = 0, end = 0.6, gamma = 2.2, alpha = NULL)}
    ###### for each chromosome
    for(k in 1:length(dd2)){
      # for each putative point which is the closest position
      ooo <- which(dd == dd2[k])
      lines(y=c(1-chr[k],1-chr[k]), x=c(sss,sss-(fact/3)), lwd=lwd, col=hc[ooo])
    }
    lines(y=c(1,1-max(chr)), x=c(sss,sss), lwd=3)
    lines(y=c(1,1-max(chr)), x=c(sss-(fact/3),sss-(fact/3)), lwd=3)
    text(x=sss-(fact/1.6), y=ruler, labels=ruler2, cex=cex) # cex of the cM ruler
    axis(3,at=(sss-(fact/3)) , labels=paste("LG",j, sep=""), cex.axis=cex.axis, font=2) # cex of the name of LGs
    #axis(2,at=0.275, labels="Density")
    ## --------------------------------------------------------
    draw.arc((sss + sss-(fact/3))/2, 1- max(chr), (sss - (sss-(fact/3)))/2, deg1=180, deg2=360, col="black", lwd=2, lend=1)
    draw.arc((sss + sss-(fact/3))/2, 1, (sss - (sss-(fact/3)))/2, deg1=0, deg2=180, col="black", lwd=2, lend=1)
    ## ----------------
    ## ----------------
    if(!is.null(trait)){
      
      chuy <- which(colnames(prov) %in% trait)
      if(length(chuy)==0){stop("The column indicated as trait is not present in the data provided\nPlease double check your data\n",call. = FALSE)}
      ## ---------------------------
      ## IF TRAIT IS A NUMERIC TRAIT
      ## ---------------------------
      if(is.numeric(prov[,trait])){
        w1 <- which(names(data) == trait) # column of trait to plot dots
        
        if(trait.scale == "same"){
          bobo <- max(data[,trait], na.rm=TRUE)
        }else{bobo <- max(prov[,trait], na.rm=TRUE)}
        
        
        dotss <- fact * ( prov[,trait]/ bobo )
        dotss2 <- sss + (fact/8) + dotss
        ####### for ablines in the trait selected
        sections <- (bobo - min(prov[,trait], na.rm=TRUE))/5
        sections2 <- seq(min(prov[,trait], na.rm=TRUE), bobo, by=sections)
        sections3 <- fact * ( sections2/ bobo )
        sections4 <- sss + (fact/8) + sections3
        ####### if trait is true
        # ablines for different values
        for(d in 1:length(sections4)){
          lines( x=c(sections4[d], sections4[d]), y=c(1,1-max(chr)), col="black", lty=3,lwd=0.5)
          text(x=sections4[d], y=1, labels=round(sections2[d],1), cex=0.4, srt=270)
        }
        
        # plot the trait values
        if(type=="dot"){
          points(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
        }
        if(type == "line"){
          polygon(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.4))
          lines(y=1-chr, x=dotss2, pch=20, cex=cex.trait, col=transp(col.trait[j],0.6))
          # density()
          
        }
        if(type == "hist"){
          for(l in 1:length(dotss2)){
            # y is the position in the chromosome
            # x is how long is the line
            lines( x=c(sss + (fact/8),dotss2[l]), y=c(1-chr[l],1-chr[l]),lwd=cex.trait, col=transp(col.trait[j],0.8))
          }
        }
        axis(3,at=sss + (fact/2), labels=trait, cex.axis=cex.axis) # cex of the scale of the trait
      }
      ## ---------------------------
      ## IF TRAIT IS A FACTOR TRAIT
      ## ---------------------------
      if(is.factor(prov[,trait])){
        riel <- sss + (fact/2)
        lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel+(fact/2), y=yy, labels=prov[ww2,trait], cex=0.3)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex) # cex of the scale of the trait (letters)
      }
      ## ---------------------------
      ## ---------------------------
      
      ## ---------------------------
      ## IF TRAIT IS A character TRAIT
      ## ---------------------------
      if(is.character(prov[,trait])){
        riel <- sss + (fact/1.8)
        #lines( x=c(riel, riel), y=c(1,1-max(chr)), col=transp(col.trait[j],0.8))
        ww2 <- which(!is.na(prov[,trait]))
        ##
        if(length(ww2) > 0){
          yy <- 1-chr[ww2]
          #points(y=yy, x=rep(riel,length(yy)), pch=15, cex=0.9, col=transp(col.trait[j],0.6))
          text(x=riel, y=yy, labels=prov[ww2,trait], cex=0.17)
        }
        ##
        axis(3,at=riel, labels=trait, cex.axis=cex.axis)
      }
      ## ---------------------------
      ## ---------------------------
    }
    ################
    
  }
}


manhattan <- function(map, col=NULL, fdr.level=0.05, show.fdr=TRUE, PVCN=NULL, ylim=NULL,...){
  
  if(!is.null(PVCN)){
    colnames(map)[which(colnames(map)==PVCN)] <- "p.val"
  }
  
  required.names <- c("Chrom","Position","p.val")
  
  che <- which(names(map)%in%required.names)
  if(length(che) < 3){
    stop("Column names; 'Chrom','Position' and 'p.val' need 
         to be present in the data frame provided.",call. = FALSE)
  }
  
  map <- map[with(map, order(Chrom, Position)), ]
  
  
  yylim <- ceiling(max(map$p.val,na.rm=TRUE))
  if(is.null(col)){
    col.scheme <- rep((transp(c("cadetblue","red"))),30)
  }else{
    col.scheme <- rep(col,30)
  }
  ffr <- fdr(map$p.val, fdr.level=fdr.level)$fdr.10
  
  if(!is.null(ylim)){
    yyylim <- ylim
  }else{
    yyylim <- c(0,yylim)
  }
  
  plot(map$p.val, bty="n", 
       col=col.scheme[factor(map$Chrom, levels = unique(map$Chrom, na.rm=TRUE))], 
       xaxt="n", xlab="Chromosome", ylab=expression(paste(-log[10],"(p.value)")), 
       las=2, 
       ylim=yyylim,...)
  init.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x)[1]; return(z)}, y=map$Chrom)
  fin.mrks <- apply(data.frame(unique(map$Chrom)),1,function(x,y){z <- which(y == x);z2 <- z[length(z)]; return(z2)}, y=map$Chrom)
  inter.mrks <- init.mrks + ((fin.mrks - init.mrks)/2)
  axis(side=1, at=inter.mrks, labels=paste("Chr",unique(map$Chrom),sep=""), cex.axis=.5)
  if(show.fdr){
    abline(h=ffr, col="slateblue4", lty=3, lwd=2)
    legend("topright", legend=paste("FDR(",fdr.level,")=",round(ffr,2), sep=""), 
           bty="n", lty=3, lwd=2, col="slateblue4", cex=0.8)
  }
  #abline(v=marker, lty=3, col="red")
  #legend("topleft", bty="n", legend=c("Markers selected"), cex=.6, lty=3, lwd=2, col="red")
  #marker <- as.character(map$Locus[marker])
}

neMarker <- function(M, neExplore=NULL, maxMarker=1000, nSamples=5){
  # maxMarker argument: only used a limited number of markers to avoid this to be too time consuming
  v <- sample(1:ncol(M), min(c(maxMarker, ncol(M))))
  M <- M[,v]
  # calculate the total number of alleles in the population
  nAllelesPop <- apply(M,2,function(x){ifelse(length(table(x)) > 1, 2, 1)})
  nAllelesPopTotal <- sum(nAllelesPop)
  # maxNe argument: define the range to test
  if(is.null(neExplore)){neExplore <- seq(10,100,10)}
  # do the sampling algorithm
  counter <- 1
  allelesCovered <- allelesCoveredSe <- vector(mode="numeric", length = length(neExplore) )
  for(i in neExplore){ # for a possible Ne
    message(paste("Exploring allele coverage (%) at Ne:",i))
    allelesCoveredSample <- vector(mode="numeric", length = nSamples)
    # nSamples argument: take a couple of samples 
    for(j in 1:nSamples){
      ii <- sample(1:nrow(M),i) # sample i individuals
      nAllelesPopI <- apply(M[ii,],2,function(x){ifelse(length(table(x)) > 1, 2, 1)}) # how many alleles we collect in the sample
      allelesCoveredSample[j] <- sum(nAllelesPopI) # sum them up
    }
    allelesCovered[counter] <- mean(allelesCoveredSample)/nAllelesPopTotal # mean across samples
    allelesCoveredSe[counter] <- ( sd(allelesCoveredSample/nAllelesPopTotal) ) # SE across samples 
    counter <- counter+1
  }
  # save results
  result <- data.frame(allelesCovered=allelesCovered, allelesCoveredSe=allelesCoveredSe, Ne=neExplore)
  return(result)
}

overlay<- function (..., rlist = NULL, prefix = NULL, sparse=FALSE){
  init <- list(...) # init <- list(DT$femalef,DT$malef)
  ## keep track of factor variables
  myTypes <- unlist(lapply(init,class))
  init0 <- init
  ##
  init <- lapply(init, as.character)
  namesInit <- as.character(substitute(list(...)))[-1L] # names <- c("femalef","malef")
  dat <- as.data.frame(do.call(cbind, init))
  dat <- as.data.frame(dat)
  ## bring back the levels
  for(j in 1:length(myTypes)){
    if(myTypes[j]=="factor"){
      levels(dat[,j]) <- c(levels(dat[,j]),setdiff(levels(init0[[j]]),levels(dat[,j]) ))
    }
  }
  ##
  if (is.null(dim(dat))) {
    stop("Please provide a data frame to the overlay function, not a vector.\\n",
         call. = FALSE)
  }
  if (is.null(rlist)) {
    rlist <- as.list(rep(1, dim(dat)[2]))
  }
  ss1 <- colnames(dat)
  dat2 <- as.data.frame(dat[, ss1])
  head(dat2)
  colnames(dat2) <- ss1
  femlist <- list()
  S1list <- list()
  for (i in 1:length(ss1)) {
    femlist[[i]] <- ss1[i]
    dat2[, femlist[[i]]] <- as.factor(dat2[, femlist[[i]]])
    if(sparse){
      S1 <- Matrix::sparse.model.matrix(as.formula(paste("~", femlist[[i]],
                                                         "-1")), dat2)
    }else{
      S1 <- model.matrix(as.formula(paste("~", femlist[[i]],
                                          "-1")), dat2)
    }
    colnames(S1) <- gsub(femlist[[i]], "", colnames(S1))
    S1list[[i]] <- S1
  }
  levo <- sort(unique(unlist(lapply(S1list, function(x) {
    colnames(x)
  }))))
  if(sparse){
    S3 <- Matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }else{
    S3 <- matrix(0, nrow = dim(dat2)[1], ncol = length(levo))
  }
  
  rownames(S3) <- rownames(dat2)
  colnames(S3) <- levo
  for (i in 1:length(S1list)) {
    if (i == 1) {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S1list[[i]] *
        rlist[[i]]
    }
    else {
      S3[rownames(S1list[[i]]), colnames(S1list[[i]])] <- S3[rownames(S1list[[i]]),
                                                             colnames(S1list[[i]])] + (S1list[[i]][rownames(S1list[[i]]),
                                                                                                   colnames(S1list[[i]])] * rlist[[i]])
    }
  }
  if (!is.null(prefix)) {
    colnames(S3) <- paste(prefix, colnames(S3), sep = "")
  }
  attr(S3,"variables") <- namesInit
  return(S3)
}

propMissing <- function(x){
  length(which(is.na(x)))/length(x)
}

redmm <- function (x, M = NULL, Lam=NULL, nPC=50, cholD=FALSE, returnLam=FALSE) {
  
  if(system.file(package = "RSpectra") == ""){
    stop("Please install the RSpectra package to use the redmm() function.",call. = FALSE)
  }else{
    requireNamespace("RSpectra",quietly=TRUE)
  }
  
  if(is.null(M)){
    # stop("M cannot be NULL. We need a matrix of features that defines the levels of x")
    smd <- RSpectra::svds(x, k=nPC, which = "LM")
    if(is.null(Lam)){
      Lam0 <- smd$u
      Lam = Lam0[,1:min(c(nPC,ncol(x))), drop=FALSE]
      rownames(Lam) <- rownames(x)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }else{
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
    Zstar <- Lam
  }else{
    
    if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
      notPresentInM <- setdiff(colnames(Z),rownames(M))
      notPresentInZ <- setdiff(rownames(M),colnames(x))
    }else{
      notPresentInM <- setdiff(unique(x),rownames(M))
      notPresentInZ <- setdiff(rownames(M),unique(x))
    }
    if(is.null(Lam)){ # user didn't provide a Lambda matrix
      if(nPC == 0){ # user wants to use the full marker matrix
        Lam <- Lam0 <- M
      }else{ # user wants to use the PCA method
        nPC <- min(c(nPC, ncol(M)))
        if(cholD){
          smd <- try(chol(M) , silent = TRUE)
          if(inherits(smd, "try-error")){smd <- try(chol((M+diag(1e-5,nrow(M),nrow(M))) ) , silent = TRUE)}
          Lam0 = t(smd)
        }else{
          smd <- RSpectra::svds(M, k=nPC, which = "LM")
          Lam0 <- smd$u
        }
        Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
        rownames(Lam) <- rownames(M)
        colnames(Lam) <- paste0("nPC",1:nPC)
      }
    }else{ # user provided it's own Lambda matrix
      Lam0=Lam
      Lam = Lam0[,1:min(c(nPC,ncol(M))), drop=FALSE]
      rownames(Lam) <- rownames(M)
      colnames(Lam) <- paste0("nPC",1:nPC)
    }
  }
  if (inherits(x, "dgCMatrix") | inherits(x, "matrix")) {
    Z <- x
  }else{
    if (!is.character(x) & !is.factor(x)) {
      namess <- as.character(substitute(list(x)))[-1L]
      Z <- Matrix(x, ncol = 1)
      colnames(Z) <- namess
    }else {
      dummy <- x
      levs <- na.omit(unique(dummy))
      if (length(levs) > 1) {
        Z <- Matrix::sparse.model.matrix(~dummy - 1, na.action = na.pass)
        colnames(Z) <- gsub("dummy", "", colnames(Z))
      } else {
        vv <- which(!is.na(dummy))
        Z <- Matrix(0, nrow = length(dummy))
        Z[vv, ] <- 1
        colnames(Z) <- levs
      }
    }
  }
  
  if(is.null(M)){
    Zstar <- Lam
  }else{
    Zstar <- as.matrix(Z %*% Lam[colnames(Z),])
  }
  
  if(returnLam){
    return(list(Z = Zstar, Lam=Lam, Lam0=Lam0)) 
  }else{return(Zstar)}
  
}

replaceValues <- function (Source, Search, Replace) {
  if (length(Search) != length(Replace)) 
    stop("Search and Replace Must Have Equal Number of Items\n")
  Changed <- as.character(Source)
  for (i in 1:length(Search)) {
    Changed <- replace(Changed, Changed == Search[i], Replace[i])
  }
  if (is.numeric(Replace)) 
    Changed <- as.numeric(Changed)
  return(Changed)
}

rrm <- function(x=NULL, H=NULL, nPC=2, returnGamma=FALSE, cholD=TRUE){
  if(is.null(x) ){stop("Please provide the x argument.", call. = FALSE)}
  if(is.null(H) ){stop("Please provide the x argument.", call. = FALSE)}
  # these are called PC models by Meyer 2009, GSE. This is a reduced rank implementation
  # we produce loadings, the Z*L so we can use it to estimate factor scores in mmec()
  
  Y <- apply(H,2, imputev)
  Ys <- scale(Y, scale = TRUE, center = TRUE)
  nans <- which(is.nan(Ys), arr.ind = TRUE)
  if(nrow(nans) > 0){
    Ys[nans]=0
  }
  Sigma <- cov(Ys) # surrogate of unstructured matrix to start with
  Sigma <- as.matrix(Matrix::nearPD(x=Sigma, corr = FALSE, keepDiag = FALSE, base.matrix = FALSE,
                                    do2eigen = TRUE, doSym = FALSE,
                                    doDykstra = TRUE, only.values = FALSE,
                                    ensureSymmetry = !isSymmetric(Sigma),
                                    eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                                    maxit = 100, conv.norm.type = "I", trace = FALSE)$mat)
  # GE <- as.data.frame(t(scale( t(scale(Y, center=T,scale=F)), center=T, scale=F)))  # sum(GE^2)
  if(cholD){
    ## OPTION 2. USING CHOLESKY
    Gamma <- t(chol(Sigma)); # LOADINGS  # same GE=LL' from cholesky  plot(unlist(Gamma%*%t(Gamma)), unlist(GE))
    D=diag(nrow(Gamma))
  }else{
    ## OPTION 1. USING SVD
    U <- svd(Sigma)$u;  # V <- svd(GE)$v
    D <- diag(svd(Sigma)$d)
    Gamma <- U %*% sqrt(D); # LOADINGS
    rownames(Gamma) <- colnames(Gamma) <- rownames(Sigma)
  }
  colnamesGamma <- colnames(Gamma)
  rownamesGamma <- rownames(Gamma)
  Gamma <- Gamma[,1:nPC, drop=FALSE]; 
  colnames(Gamma) <- colnamesGamma[1:nPC]
  rownames(Gamma) <- rownamesGamma
  ##
  rownames(Gamma) <- gsub("v.names_","",rownames(Gamma))#rownames(GE)#levels(dataset$Genotype);  # rownames(Se) <- colnames(GE)#levels(dataset$Environment)
  colnames(Gamma) <- paste("PC", 1:ncol(Gamma), sep =""); # 
  ######### GEreduced = Sg %*% t(Se) 
  # if we want to merge with PCs for environments
  dtx <- data.frame(timevar=x)
  dtx$index <- 1:nrow(dtx)
  Z <- Matrix::sparse.model.matrix(~timevar -1, na.action = na.pass, data=dtx)
  colnames(Z) <- gsub("timevar","",colnames(Z))
  Zstar <- Z%*%Gamma[colnames(Z),] # we multiple original Z by the LOADINGS
  Zstar <- as.matrix(Zstar)
  rownames(Z) <- NULL
  
  if(returnGamma){
    return(list(Gamma=Gamma, H=H, Sigma=Sigma, Zstar=Zstar, D=D))
  }else{
    return(Zstar)
  }
}

simGECorMat <- function (nEnv, nMegaEnv, mu = 0.7, v = 0.2, mu2 = 0, v2 = 0.3) {
  ff <- function(m) {
    m[lower.tri(m)] <- t(m)[lower.tri(m)]
    m
  }
  G = matrix(NA, nEnv, nEnv)
  nEnv2 <- nEnv/nMegaEnv
  
  starts <- seq(1, nEnv, nEnv/nMegaEnv)
  ends <- c((starts - 1)[-1], nEnv)
  for (i in 1:nMegaEnv) {
    corsprov <- rnorm((nEnv2 * (nEnv2 - 1))/2, mu, v)
    counter = 1
    for (j in starts[i]:ends[i]) {
      for (k in j:ends[i]) {
        if (j == k) {
          G[j, k] <- 1
        }
        else {
          G[j, k] <- corsprov[counter]
          counter <- counter + 1
        }
      }
    }
  }
  G <- ff(G)
  tofill <- which(is.na(G), arr.ind = TRUE)
  G[tofill] <- rnorm(nrow(tofill), mu2, v2)
  G[which((G) > 1)] <- 0.98
  G[which((G) < -1)] <- -0.98
  G <- ff(G)
  return(G)
}

simage <- function(data, Var1=NULL, Var2=NULL, ...){
  M <- table(data[,Var1], data[,Var2])
  M2 <- matrix(M, nrow = nrow(M), ncol = ncol(M))
  name1 <- as.character(substitute(list(Var1)))[-1L]
  name2 <- as.character(substitute(list(Var2)))[-1L]
  Matrix::image(as(as(as( t(M2) ,  "dMatrix"), "generalMatrix"), "CsparseMatrix"), #as(t(M2), Class = "dgCMatrix"), 
                xlab=name1, ylab=name2, 
                colorkey=TRUE, ...)
}  

simage2 <- function(X, ...){
  
  Matrix::image(as(as(as( X ,  "dMatrix"), "generalMatrix"), "CsparseMatrix"), #as(t(M2), Class = "dgCMatrix"), 
                # xlab=name1, ylab=name2, 
                colorkey=TRUE, ...)
}  

stan <-function (x, lb=0, ub=1) {
  if(var(x, na.rm=TRUE)>0){
    B=max(x, na.rm=TRUE) # current range
    A=min(x, na.rm=TRUE) # current range 
    D=ub # new range
    C=lb # new range
    
    scale = (D-C)/(B-A)
    offset = -A*(D-C)/(B-A) + C
    return(x*scale + offset)
  }else{
    return(x)
  }
}

tps <- function (columncoordinates, rowcoordinates, nsegments=NULL, 
                 minbound=NULL, maxbound=NULL, degree = c(3, 3), penaltyord = c(2, 2), 
                 nestorder = c(1, 1), asreml = "grp", eigenvalues = "include", 
                 method = "Lee", stub = NULL) 
{
  if (missing(columncoordinates)) 
    stop("columncoordinates argument must be set")
  if (missing(rowcoordinates)) 
    stop("rowcoordinates argument must be set")
  col <- columncoordinates
  nuc <- length(col)
  col.match <- match(columncoordinates, col)
  row <- sort(unique(rowcoordinates))
  nur <- length(row)
  row.match <- match(rowcoordinates, row)
  nv <- length(columncoordinates)
  if (is.null(minbound)) {
    cminval <- min(col)
    rminval <- min(row)
  } else {
    cminval <- min(c(minbound[1], min(col)))
    if (length(minbound) < 2) {
      rminval <- min(c(minbound[1], min(row)))
    }
    else {
      rminval <- min(c(minbound[2], min(row)))
    }
  }
  if (is.null(maxbound)) {
    cmaxval <- max(col)
    rmaxval <- max(row)
  }
  else {
    cmaxval <- max(c(maxbound[1], max(col)))
    if (length(maxbound) < 2) {
      rmaxval <- max(c(maxbound[1], max(row)))
    }
    else {
      rmaxval <- max(c(maxbound[2], max(row)))
    }
  }
  if (is.null(nsegments)) {
    nsegcol <- nuc - 1
    nsegrow <- nur - 1
  }
  else {
    nsegcol <- max(c(nsegments[1], 2))
  }
  if (length(nsegments) < 2) {
    nsegrow <- max(c(nsegments[1], 2))
  }
  else {
    nsegrow <- max(c(nsegments[2], 2))
  }
  nestcol <- floor(nestorder[1])
  if (length(nestorder) < 2) 
    nestrow <- floor(nestorder[1])
  else nestrow <- floor(nestorder[2])
  nsncol <- 0
  if (nestcol > 1) {
    if (nsegcol%%nestcol != 0) 
      warning("Column nesting ignored: number of column segments must be a multiple of nesting order")
    else nsncol <- nsegcol/nestcol
  }
  nsnrow <- 0
  if (nestrow > 1) {
    if (nsegrow%%nestrow != 0) 
      warning("Row nesting ignored: number of row segments must be a multiple of nesting order")
    else nsnrow <- nsegrow/nestrow
  }
  Bc <- bbasis(col, cminval, cmaxval, nsegcol, degree[1])
  nc <- ncol(Bc)
  if (length(degree) < 2) 
    degr <- degree[1]
  else degr <- degree[2]
  Br <- bbasis(row, rminval, rmaxval, nsegrow, degr)
  nr <- ncol(Br)
  if (nsncol > 0) {
    Bcn <- bbasis(col, cminval, cmaxval, nsncol, degree[1])
    ncn <- ncol(Bcn)
  }
  else ncn <- nc
  if (nsnrow > 1) {
    Brn <- bbasis(row, rminval, rmaxval, nsnrow, degr)
    nrn <- ncol(Brn)
  }
  else nrn <- nr
  diff.c <- penaltyord[[1]]
  Dc <- diff(diag(nc), diff = diff.c)
  svd.c <- svd(crossprod(Dc))
  nbc <- nc - diff.c
  U.Zc <- svd.c$u[, c(1:nbc)]
  U.Xc <- svd.c$u[, -c(1:nbc)]
  L.c <- sqrt(svd.c$d[c(1:nbc)])
  diagc <- L.c^2
  BcU <- Bc %*% U.Zc
  BcX <- Bc %*% U.Xc
  BcULi <- BcU %*% diag(1/L.c)
  if ("include" %in% eigenvalues) {
    BcZmat.df <- as.data.frame(BcULi)
    BcZmat <- BcULi
  }
  else {
    BcZmat.df <- as.data.frame(BcU)
    BcZmat <- BcU
  }
  BcZmat.df$TP.col <- col
  mat1c <- matrix(rep(1, nuc), nrow = nuc)
  BcXadj <- BcX - mat1c %*% t(mat1c) %*% BcX/nuc
  Xfc <- (svd(crossprod(BcXadj)))$u[, c(ncol(BcXadj):1)]
  BcX <- BcX %*% Xfc
  if (BcX[1, 1] < 0) 
    BcX[, 1] <- -1 * BcX[, 1]
  if (BcX[1, 2] > 0) 
    BcX[, 2] <- -1 * BcX[, 2]
  if (nsncol > 0) {
    Dcn <- diff(diag(ncn), diff = diff.c)
    svd.cn <- svd(crossprod(Dcn))
    nbcn <- ncn - diff.c
    U.Zcn <- svd.cn$u[, c(1:nbcn)]
    U.Xcn <- svd.cn$u[, -c(1:nbcn)]
    L.cn <- sqrt(svd.cn$d[c(1:nbcn)])
    BcnU <- Bcn %*% U.Zcn
    BcnX <- Bcn %*% U.Xcn
  }
  else {
    nbcn <- nbc
    BcnU <- BcU
    L.cn <- L.c
  }
  if (length(penaltyord) < 2) {
    diff.r <- penaltyord[1]
  }
  else {
    diff.r <- penaltyord[2]
  }
  Dr <- diff(diag(nr), diff = diff.r)
  svd.r <- svd(crossprod(Dr))
  nbr <- nr - diff.r
  U.Zr <- svd.r$u[, c(1:nbr)]
  U.Xr <- svd.r$u[, -c(1:nbr)]
  L.r <- sqrt(svd.r$d[c(1:nbr)])
  diagr <- L.r^2
  BrU <- Br %*% U.Zr
  BrX <- Br %*% U.Xr
  BrULi <- BrU %*% diag(1/L.r)
  if ("include" %in% eigenvalues) {
    BrZmat.df <- as.data.frame(BrULi)
    BrZmat <- BrULi
  }
  else {
    BrZmat.df <- as.data.frame(BrU)
    BrZmat <- BrU
  }
  BrZmat.df$TP.row <- row
  mat1r <- matrix(rep(1, nur), nrow = nur)
  BrXadj <- BrX - mat1r %*% t(mat1r) %*% BrX/nur
  Xfr <- (svd(crossprod(BrXadj)))$u[, c(ncol(BrXadj):1)]
  BrX <- BrX %*% Xfr
  if (BrX[1, 1] < 0) 
    BrX[, 1] <- -1 * BrX[, 1]
  if (BrX[1, 2] > 0) 
    BrX[, 2] <- -1 * BrX[, 2]
  if (nsnrow > 0) {
    Drn <- diff(diag(nrn), diff = diff.r)
    svd.rn <- svd(crossprod(Drn))
    nbrn <- nrn - diff.r
    U.Zrn <- svd.rn$u[, c(1:nbrn)]
    U.Xrn <- svd.rn$u[, -c(1:nbrn)]
    L.rn <- sqrt(svd.rn$d[c(1:nbrn)])
    BrnU <- Brn %*% U.Zrn
    BrnX <- Brn %*% U.Xrn
  }
  else {
    nbrn <- nbr
    BrnU <- BrU
    L.rn <- L.r
  }
  A <- 10^(floor(log10(max(row))) + 1)
  row.index <- rep(row, times = nuc)
  col.index <- rep(col, each = nur)
  index <- A * col.index + row.index
  C.R <- A * columncoordinates + rowcoordinates
  BcrZ1 <- BcnU[col.match, ] %x% matrix(rep(1, nbrn), nrow = 1, 
                                        ncol = nbrn)
  BcrZ2 <- matrix(rep(1, nbcn), nrow = 1, ncol = nbcn) %x% 
    BrnU[row.match, ]
  BcrZ <- BcrZ1 * BcrZ2
  diagrx <- rep(L.cn^2, each = nbrn)
  diagcx <- rep(L.rn^2, times = nbcn)
  if ("Lee" %in% method) {
    diagcr <- diagrx + diagcx
  }
  if ("Wood" %in% method) {
    diagcr <- diagrx * diagcx
  }
  if (!("Lee" %in% method) & !("Wood" %in% method)) {
    stop("Invalid setting of method argument")
  }
  BcrZLi <- BcrZ %*% diag(1/sqrt(diagcr))
  if ("include" %in% eigenvalues) {
    BcrZmat.df <- as.data.frame(BcrZLi)
    BcrZmat <- BcrZLi
  }
  else {
    BcrZmat.df <- as.data.frame(BcrZ)
    BcrZmat <- BcrZ
  }
  BcrZmat.df$TP.CxR <- C.R
  tracelist <- list()
  for (i in 1:diff.c) {
    nm <- paste0("Xc", i, ":Zr")
    tempmat <- (BcX[col.match, i] %x% matrix(rep(1, nbr), 
                                             nrow = 1)) * BrZmat[row.match, ]
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagr), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  for (i in 1:diff.r) {
    nm <- paste0("Zc:Xr", i)
    tempmat <- BcZmat[col.match, ] * (matrix(rep(1, nbc), 
                                             nrow = 1) %x% BrX[row.match, i])
    if ("include" %in% eigenvalues) 
      tempmatsc <- tempmat
    else tempmatsc <- tempmat * (rep(1, nv) %*% matrix((1/diagc), 
                                                       nrow = 1))
    tracelist[nm] <- sum(tempmatsc * tempmat)
  }
  if ("include" %in% eigenvalues) 
    tracelist["Zc:Zr"] <- sum(BcrZmat * BcrZmat)
  else {
    tempmatsc <- BcrZmat * (rep(1, nv) %*% matrix((1/diagcr), 
                                                  nrow = 1))
    tracelist["Zc:Zr"] <- sum(tempmatsc * BcrZmat)
  }
  # outdata <- as.data.frame(data)
  outdata <- data.frame(TP.col=columncoordinates)
  outdata$TP.row <- rowcoordinates
  outdata$TP.CxR <- C.R
  BcrX1 <- BcX[col.match, ] %x% matrix(rep(1, diff.r), nrow = 1)
  BcrX2 <- matrix(rep(1, diff.c), nrow = 1) %x% BrX[row.match, 
  ]
  BcrX <- BcrX1 * BcrX2
  fixed <- list()
  fixed$col <- data.frame(row.names = C.R)
  for (i in 1:diff.c) {
    c.fixed <- paste("TP.C", ".", i, sep = "")
    outdata[c.fixed] <- BcX[col.match, i]
    fixed$col[c.fixed] <- BcX[col.match, i]
  }
  fixed$row <- data.frame(row.names = C.R)
  for (i in 1:diff.r) {
    r.fixed <- paste("TP.R", ".", i, sep = "")
    outdata[r.fixed] <- BrX[row.match, i]
    fixed$row[r.fixed] <- BrX[row.match, i]
  }
  ncolX <- diff.c * diff.r
  fixed$int <- data.frame(row.names = C.R)
  for (i in 1:ncolX) {
    cr.fixed <- paste("TP.CR", ".", i, sep = "")
    outdata[cr.fixed] <- BcrX[, i]
    fixed$int[cr.fixed] <- BcrX[, i]
  }
  if (!missing(stub)) {
    cname <- paste0("BcZ", stub, ".df")
    rname <- paste0("BrZ", stub, ".df")
    crname <- paste0("BcrZ", stub, ".df")
  }
  else {
    cname <- "BcZ.df"
    rname <- "BrZ.df"
    crname <- "BcrZ.df"
  }
  mbftext <- paste0("list(TP.col=list(key=c(\"TP.col\",\"TP.col\"),cov=\"", 
                    cname, "\"),")
  mbftext <- paste0(mbftext, "TP.row=list(key=c(\"TP.row\",\"TP.row\"),cov=\"", 
                    rname, "\"),")
  mbftext <- paste0(mbftext, "TP.CxR=list(key=c(\"TP.CxR\",\"TP.CxR\"),cov=\"", 
                    crname, "\"))")
  mbflist <- eval(parse(text = mbftext))
  if ("grp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    start0 <- start
    scale <- 1
    j <- 1
    for (i in 1:diff.c) {
      nm0 <- paste0(names(fixed$col[i]), "_frow")
      listnames[j] <- nm0
      for (k in 1:nbr) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$col[[i]] * BrZmat[row.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbr, 
                      by = 1)
      start <- start + nbr
      j <- j + 1
    }
    for (i in 1:diff.r) {
      nm0 <- paste0(names(fixed$row[i]), "_fcol")
      listnames[j] <- nm0
      for (k in 1:nbc) {
        nm <- paste0(nm0, "_", k)
        outdata[nm] <- scale * fixed$row[[i]] * BcZmat[col.match, 
                                                       k]
      }
      grp[[j]] <- seq(from = start + 1, to = start + nbc, 
                      by = 1)
      start <- start + nbc
      j <- j + 1
    }
    m <- 0
    nm0 <- "TP_fcol_frow"
    listnames[j] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- scale * BcrZmat[, k]
    }
    grp[[j]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    end <- start + (nbcn * nbrn)
    j <- j + 1
    listnames[j] <- "All"
    grp[[j]] <- seq(from = start0 + 1, to = end, by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("sepgrp" %in% asreml) {
    grp <- list()
    listnames <- list()
    start <- length(outdata)
    nm0 <- "TP_C"
    listnames[1] <- nm0
    for (i in 1:diff.c) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$col[[i]]
    }
    grp[[1]] <- seq(from = start + 1, to = start + diff.c, 
                    by = 1)
    start <- start + diff.c
    nm0 <- "TP_R"
    listnames[2] <- nm0
    for (i in 1:diff.r) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- fixed$row[[i]]
    }
    grp[[2]] <- seq(from = start + 1, to = start + diff.r, 
                    by = 1)
    start <- start + diff.r
    nm0 <- "TP_fcol"
    listnames[3] <- nm0
    for (k in 1:nbc) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcZmat[col.match, k]
    }
    grp[[3]] <- seq(from = start + 1, to = start + nbc, by = 1)
    start <- start + nbc
    nm0 <- "TP_frow"
    listnames[4] <- nm0
    for (k in 1:nbr) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BrZmat[row.match, k]
    }
    grp[[4]] <- seq(from = start + 1, to = start + nbr, by = 1)
    start <- start + nbr
    grp <- structure(grp, names = listnames)
    nm0 <- "TP_fcol_frow"
    listnames[5] <- nm0
    for (k in 1:(nbrn * nbcn)) {
      nm <- paste0(nm0, "_", k)
      outdata[nm] <- BcrZmat[, k]
    }
    grp[[5]] <- seq(from = start + 1, to = start + (nbcn * 
                                                      nbrn), by = 1)
    grp <- structure(grp, names = listnames)
  }
  if ("own" %in% asreml) {
    grp <- list()
    listnames <- list()
    listnames[1] <- "All"
    start <- length(outdata)
    nm0 <- "Xc_Zr"
    Xc_Zr <- (BcX[col.match, ] %x% matrix(rep(1, nbr), nrow = 1)) * 
      (matrix(rep(1, diff.c), nrow = 1) %x% BrZmat[row.match, 
      ])
    nXc_Zr <- ncol(Xc_Zr)
    for (i in 1:nXc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Xc_Zr[, i]
    }
    nm0 <- "Zc_Xr"
    Zc_Xr <- (BcZmat[col.match, ] %x% matrix(rep(1, diff.r), 
                                             nrow = 1)) * (matrix(rep(1, nbc), nrow = 1) %x% BrX[row.match, 
                                             ])
    nZc_Xr <- ncol(Zc_Xr)
    for (i in 1:nZc_Xr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Xr[, i]
    }
    nm0 <- "Zc_Zr"
    Zc_Zr <- BcrZmat
    nZc_Zr <- ncol(Zc_Zr)
    for (i in 1:nZc_Zr) {
      nm <- paste0(nm0, "_", i)
      outdata[nm] <- Zc_Zr[, i]
    }
    grp[[1]] <- seq(from = start + 1, to = start + nXc_Zr + 
                      nZc_Xr + nZc_Zr, by = 1)
    grp <- structure(grp, names = listnames)
  }
  res <- list()
  res$data <- outdata
  res$mbflist <- mbflist
  res[["BcZ.df"]] <- BcZmat.df
  res[["BrZ.df"]] <- BrZmat.df
  res[["BcrZ.df"]] <- BcrZmat.df
  res[["All"]] <- as.matrix(outdata[,grp$All])
  res$dim <- c(diff.c = diff.c, nbc = nbc, nbcn = nbcn, diff.r = diff.r, 
               nbr = nbr, nbrn = nbrn)
  res$trace <- tracelist
  if ("grp" %in% asreml) 
    res$grp <- grp
  if ("sepgrp" %in% asreml) 
    res$grp <- grp
  if ("own" %in% asreml) 
    res$grp <- grp
  if ("mbf" %in% asreml) 
    res$grp <- NULL
  if (!("include" %in% eigenvalues)) 
    res$eigen <- list(diagc = diagc, diagr = diagr, diagcr = diagcr)
  res
}

transp <- function (col, alpha = 0.5) {
  res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, c[2]/255, 
                                                c[3]/255, alpha))
  return(res)
}

wald.test <- function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL,
                       verbose = FALSE) {
  if (is.null(Terms) & is.null(L))
    stop("One of the arguments Terms or L must be used.")
  if (!is.null(Terms) & !is.null(L))
    stop("Only one of the arguments Terms or L must be used.")
  if (is.null(Terms)) {
    w <- nrow(L)
    Terms <- seq(length(b))[colSums(L) > 0]
  }
  else w <- length(Terms)
  if (is.null(H0))
    H0 <- rep(0, w)
  if (w != length(H0))
    stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
  if (is.null(L)) {
    L <- matrix(rep(0, length(b) * w), ncol = length(b))
    for (i in 1:w) L[i, Terms[i]] <- 1
  }
  dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),
                            sep = ""), names(b))
  f <- L %*% b
  V <- Sigma
  mat <- qr.solve(L %*% V %*% t(L))
  stat <- t(f - H0) %*% mat %*% (f - H0)
  p <- 1 - pchisq(stat, df = w)
  if (is.null(df))
    res <- list(chi2 = c(chi2 = stat, df = w, P = p))
  else {
    fstat <- stat/nrow(L)
    df1 <- nrow(L)
    df2 <- df
    res <- list(chi2 = c(chi2 = stat, df = w, P = p), Ftest = c(Fstat = fstat,
                                                                df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
  }
  structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
                 L = L, result = res, verbose = verbose, df = df), class = "wald.test")
}

print.wald.test <- function(x, digits = 2, ...){
  Terms <- x[["Terms"]]; b <- x[["b"]]; H0 <- x[["H0"]]; v <- x[["result"]][["chi2"]]; df <- x[["df"]]
  verbose <- x[["verbose"]]
  namb <- names(b)[Terms]
  message(paste("Wald test:\n", "----------\n", sep = ""))
  if(verbose){
    message("\nCoefficients:\n")
    message(format(b, digits = digits), quote = FALSE)
    message("\nVar-cov matrix of the coefficients:\n")
    message(format(x[["Sigma"]], digits = digits), quote = FALSE)
    message("\nTest-design matrix:\n")
    message(x[["L"]])
    message("\nPositions of tested coefficients in the vector of coefficients:", paste(Terms, collapse = ", "), "\n")
    if(is.null(namb))
      message("\nH0: ", paste(paste(format(b[Terms], digits), format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
    else{
      message("\nH0: ", paste(paste(namb, format(H0, digits = digits), sep = " = "), collapse = "; "), "\n")
    }
    
  }
  message("\nChi-squared test:\n")
  message("X2 = ", format(v["chi2"], digits = digits, nsmall = 1), ", df = ", v["df"],
      ", P(> X2) = ", format(v["P"], digits = digits, nsmall = 1), "\n", sep = "")
  if(!is.null(df)){
    v <- x[["result"]][["Ftest"]]
    message("\nF test:\n")
    message("W = ", format(v["Fstat"], digits = digits, nsmall = 1),
        ", df1 = ", v["df1"],
        ", df2 = ", v["df2"],
        ", P(> W) = ", format(v["P"], digits = digits), "\n", sep = "")
  }
}












