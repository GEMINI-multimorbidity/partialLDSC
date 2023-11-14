#' partial_ldsc - main function to estimate partial genetic correlations
#'
#' Estimates unadjusted and partial genetic correlations (and heritabilities, on observed scale), 
#' and tests for the difference between the two, for all pairs of conditions.
#'
#'
#' @param conditions The path to the files containing munged the GWAS summary statistics 
#'        for the conditions (character)
#' @param confounders The path to the files containing munged the GWAS summary statistics 
#'        for the confounders (character)
#' @param condition.names The names of the conditions, \code{default=NULL} (character)
#' @param confounder.names The names of the confounders, \code{default=NULL} (character)
#' @param ld The path to the folder in which the LD scores used in the analysis are located.
#'        Expects LD scores formated as required by the original LD score regression software.  (character)
#' @param n.blocks  The number of blocks used for block-jackknife, \code{default=200} (numeric)
#' @param log.name  The name of the log file, \code{default=NULL} (character)


#' @details
#' \code{conditions} and \code{confounders} are required arguments.
#' These input files should have been pre-processed using the \code{munge()} function
#' from [\code{GenomicSEM}](https://github.com/GenomicSEM/GenomicSEM).
#' @export

#' @importFrom rlang .data

partial_ldsc <- function(conditions, confounders, 
                         condition.names = NULL, confounder.names = NULL,
                         ld, n.blocks = 200, log.name = NULL){
  
  
  .LOG <- function(..., file, print = TRUE) {
    msg <- paste0(..., "\n")
    if (print) cat(msg)
    cat(msg, file = file, append = TRUE)
  }
  
  begin.time <- Sys.time()
  begin.time.nice <- paste0(lubridate::hour(begin.time), ":", paste0(ifelse(lubridate::minute(begin.time)<10, "0", ""), lubridate::minute(begin.time)))
  #### check the parameters ####
  
  ## 1) conditions / confounder
  
  if (is.character(conditions)){
    if(any(!file.exists(conditions))) stop("conditions : some files do not exist", call. = FALSE)
    # get absolute path
    conditions = normalizePath(conditions)
  } else stop("conditions : wrong format, should be character", call. = FALSE)
  if (is.character(confounders)){
    if(any(!file.exists(confounders))) stop("confounders : some files do not exist", call. = FALSE)
    # get absolute path
    confounders = normalizePath(confounders)
  } else stop("confounders : wrong format, should be character", call. = FALSE)
  
  if(length(conditions)<2) stop("At least two conditions are needed to estimate the partial correlation.")

  n.conditions = length(conditions)
  n.confounders = length(confounders)
  
  traits = c(conditions, confounders)
  
  if(! is.null(condition.names)){
    if(! is.character(condition.names)) stop("condition.names : wrong format, should be character", call. = FALSE)
  } else { # create condition.names name if NULL
    condition.names = paste0("condition", 1:n.conditions)
  }
  if(! is.null(confounder.names)){
    if(! is.character(confounder.names)) stop("confounder.names : wrong format, should be character", call. = FALSE)
  } else { # create condition.names name if NULL
    confounder.names = paste0("confounder", 1:n.confounders)
  }
  trait.names = c(condition.names, confounder.names)
  
  
  ## 2) ld / wld
  if (is.character(ld)){
    if(!dir.exists(ld)) stop("ld : the folder does not exist", call. = FALSE)
    # get absolute path
    ld = normalizePath(ld)
    # check that files are here
    if(any(!file.exists(file.path(ld, paste0(1:22, ".l2.ldscore.gz"))))) stop("ld : expected .ldscore.gz files do not exist in the folder", call. = FALSE)
    if(any(!file.exists(file.path(ld, paste0(1:22, ".l2.M_5_50"))))) stop("ld : expected .l2.M_5_50 files do not exist in the folder", call. = FALSE)
  } else stop("ld : wrong format, should be character", call. = FALSE)
  
  # 3) log.file
  
  if(! is.null(log.name)){
    if(! is.character(log.name)) stop("log.name : wrong format, should be character", call. = FALSE)
  } else { # create log.name name if NULL
    logtraits <- gsub(".*/", "", traits)
    log.name <- paste(logtraits, collapse = "_")
    if(gdata::object.size(log.name) > 200){
      log.name <- substr(log.name,1,100)  
    }
  }
  
  log.file <- file(paste0(log.name, "_ldsc.log"), open = "wt")
  
  #### start analysis ####
  
  .LOG("Multivariate ld-score regression of ", length(traits), " traits ", 
       "(",n.conditions, " conditions: ", paste(condition.names, collapse = ", "), 
       " + ", n.confounders, " confounder(s): ", paste(confounder.names, collapse = ", "), ")", " began at: ", begin.time.nice, file = log.file)

  ## Dimensions
  n.traits <- length(traits)
  # number of elements of the unadjusted h2/gcov matrix
  n.V <- n.traits * (n.traits + 1) / 2
  # number of elements of the partial h2/gcov matrix
  n.Vpartial = n.conditions * (n.conditions + 1) / 2
  
  .LOG("Number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) is ", n.blocks, file=log.file)
  
  ## Storage
  N.vec <- matrix(NA,nrow=1,ncol=n.V)             # sample sizes
  I <- matrix(NA,nrow=n.traits,ncol=n.traits)     # intercepts (not needed, but keep in case)
  
  # leave-one-out values 
  # for h2/gcov
  V.delete <- matrix(NA, nrow = n.blocks, ncol = n.V)   
  # for rg
  V_Stand.delete = matrix(NA, nrow = n.blocks, ncol = n.V)
  # for partial h2/cov
  partial.V.delete  <- matrix(NA, nrow = n.blocks, ncol = n.Vpartial) 
  # for partial rg
  partial.V_Stand.delete = matrix(NA, nrow = n.blocks, ncol = n.Vpartial)
  
  # pseudo-values
  # for h2/gcov
  V.hold <- matrix(NA, nrow = n.blocks, ncol = n.V)  
  # for rg
  V_Stand.hold = matrix(NA, nrow = n.blocks, ncol = n.V)
  # for partial h2/cov
  partial.V.hold <- matrix(NA, nrow = n.blocks, ncol = n.Vpartial)
  # for partial rg
  partial.V_Stand.hold = matrix(NA, nrow = n.blocks, ncol = n.Vpartial)
  
  # names for V / partial.V etc...
  V.names = rep(NA_character_, n.V)
  partial.V.names <- rep(NA_character_, n.Vpartial)
  
  # genome-wide (partial) h2/gcov matrix
  cov <- matrix(NA, nrow = n.traits, ncol = n.traits) 
  partial.cov <- matrix(NA, nrow = n.conditions, ncol = n.conditions) 
  
  # genome-wide (partial) rg matrix - created later
  # S_Stand
  # partial.S_Stand
  
  # variance covariance matrix of cov / partial.cov - created later
  # V
  # partial.V
  # variance covariance matrix of S_Stand / partial.S_Stand - created later
  # V_Stand
  # partial.V_Stand
  
  # covariance between raw/partial rg
  Vcov_Stand <-matrix(NA, nrow = n.conditions, ncol = n.conditions)
  
  
  #### read LD scores / weights / M ####
  .LOG("Reading in LD scores from: ", ld, file = log.file)
  
  chr=22
  # ld scores
  x <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(readr::read_delim(
      file.path(ld, paste0(i, ".l2.ldscore.gz")),
      delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE))
  }))
  x$CM <- NULL
  x$MAF <- NULL
  
  # weights
  w <- x
  w$CM <- NULL
  w$MAF <- NULL
  colnames(w)[ncol(w)] <- "wLD"
  
  # m
  m <- do.call("rbind", lapply(1:chr, function(i) {
    suppressMessages(readr::read_csv(file.path(ld, paste0(i, ".l2.M_5_50")), col_names = FALSE))
  }))
  
  M.tot <- sum(m)
  m <- M.tot
  
  #### read summary stats + merge with LD files ####
  s <- 0
  
  chisq.max = NA # not used as an argument at the moment
  all_y <- lapply(traits, function(chi1) {
    
    ## READ chi2
    y1 <- suppressMessages(stats::na.omit(readr::read_delim(
      chi1, delim = "\t", escape_double = FALSE, trim_ws = TRUE, progress = FALSE)))
    
    .LOG("Read in summary statistics [", s <<- s + 1, "/", n.traits, "] (", trait.names[s], ") from: ", chi1, file=log.file)
    
    # make sure these are munged files, with correct colnames
    if(!all(colnames(y1) %in% c("SNP", "A1", "A2", "Z", "N"))) {
      .LOG(y1, " does not have the expected column names", file=log.file)
      flush(log.file)
      close(log.file)
      stop(paste0(chi1, " : wrong format, make sure the file has been correctly munged"), call. = FALSE)
      
      
    }
    
    ## Merge files
    merged <- merge(y1[, c("SNP", "N", "Z", "A1")], w[, c("SNP", "wLD")], by = "SNP", sort = FALSE)
    merged <- merge(merged, x, by = "SNP", sort = FALSE)
    merged <- merged[with(merged, order(CHR, BP)), ]
    
    .LOG("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files", file=log.file)
    
    ## REMOVE SNPS with excess chi-square:
    
    if(is.na(chisq.max)){
      chisq.max <- max(0.001 * max(merged$N), 80)
    }
    rm <- (merged$Z^2 > chisq.max)
    merged <- merged[!rm, ]
    
    .LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max, "; ", nrow(merged), " remain", file=log.file)
    
    merged
  })
  
  
  #### fit ldsc model ####
  
  # count the total number of runs, both loops
  s <- 1
  
  for(j in 1:n.traits){
    
    chi1 <- traits[j]
    
    y1 <- all_y[[j]]
    y1$chi1 <- y1$Z^2
    
    for(k in j:length(traits)){
      
      ##### HERITABILITY code
      
      if(j == k){
        
        .LOG("     ", "     ", file=log.file, print = FALSE)
        .LOG("Estimating heritability [", s, "/", n.V, "] for: ", trait.names[j], file=log.file)
        
        merged <- y1
        n.snps <- nrow(merged)
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        merged$weights <- merged$initial.w/sum(merged$initial.w)
        
        N.bar <- mean(merged$N)
        
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$chi1*merged$weights)
        
        
        ## Perfrom analysis:
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[,s] <- pseudo.values[,1]
        N.vec[1,s] <- N.bar
        
        V.delete[, s] = delete.values[,1]
        V.names[s] = paste(trait.names[j], trait.names[k], sep = "-")
      
        
        cov[j,j] <- reg.tot
        I[j,j] <- intercept
        
        lambda.gc <- stats::median(merged$chi1) / stats::qchisq(0.5, df = 1)
        mean.Chi <- mean(merged$chi1)
        ratio <- (intercept - 1) / (mean.Chi - 1)
        ratio.se <- intercept.se / (mean.Chi - 1)
        
        .LOG("Heritability Results for trait: ", trait.names[j], file=log.file)
        .LOG("Mean Chi^2 across remaining SNPs: ", round(mean.Chi, 4), file=log.file)
        .LOG("Lambda GC: ", round(lambda.gc, 4), file=log.file)
        .LOG("Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG("Ratio: ", round(ratio, 4), " (", round(ratio.se, 4), ")", file=log.file)
        .LOG("Total Observed Scale h2: ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG("h2 Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
      }
      
      
      ##### GENETIC COVARIANCE code
      
      if(j != k){
        
        .LOG("     ", file=log.file, print = FALSE)
        
        chi2 <- traits[k]
        .LOG("Calculating genetic covariance [", s, "/", n.V, "] for traits: ", trait.names[j], " and ", trait.names[k], file=log.file)
        
        # Reuse the data read in for heritability
        y2 <- all_y[[k]]
        y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")], by = "SNP", sort = FALSE)
        
        y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
        y$ZZ <- y$Z.y * y$Z.x
        y$chi2 <- y$Z.y^2
        merged <- stats::na.omit(y)
        n.snps <- nrow(merged)
        
        .LOG(n.snps, " SNPs remain after merging ", trait.names[j], " and ", trait.names[k], " summary statistics", file=log.file)
        
        ## ADD INTERCEPT:
        merged$intercept <- 1
        merged$x.tot <- merged$L2
        merged$x.tot.intercept <- 1
        
        
        #### MAKE WEIGHTS:
        
        tot.agg <- (M.tot*(mean(merged$chi1)-1))/mean(merged$L2*merged$N.x)
        tot.agg <- max(tot.agg,0)
        tot.agg <- min(tot.agg,1)
        merged$ld <- pmax(merged$L2, 1)
        merged$w.ld <- pmax(merged$wLD, 1)
        merged$c <- tot.agg*merged$N.x/M.tot
        merged$het.w <- 1/(2*(1+(merged$c*merged$ld))^2)
        merged$oc.w <- 1/merged$w.ld
        merged$w <- merged$het.w*merged$oc.w
        merged$initial.w <- sqrt(merged$w)
        
        tot.agg2 <- (M.tot*(mean(merged$chi2)-1))/mean(merged$L2*merged$N.y)
        tot.agg2 <- max(tot.agg2,0)
        tot.agg2 <- min(tot.agg2,1)
        merged$ld2 <- pmax(merged$L2, 1)
        merged$w.ld2 <- pmax(merged$wLD, 1)
        merged$c2 <- tot.agg2*merged$N.y/M.tot
        merged$het.w2 <- 1/(2*(1+(merged$c2*merged$ld))^2)
        merged$oc.w2 <- 1/merged$w.ld2
        merged$w2 <- merged$het.w2*merged$oc.w2
        merged$initial.w2 <- sqrt(merged$w2)
        
        
        merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w + merged$initial.w2 )
        
        N.bar <- sqrt(mean(merged$N.x)*mean(merged$N.y))
        
        ## preweight LD and chi:
        
        weighted.LD <- as.matrix(cbind(merged$L2,merged$intercept)*merged$weights)
        weighted.chi <- as.matrix(merged$ZZ *merged$weights_cov)
        
        ## Perfrom analysis:
        
        
        n.annot <- 1
        
        
        select.from <- floor(seq(from=1,to=n.snps,length.out =(n.blocks+1)))
        select.to <- c(select.from[2:n.blocks]-1,n.snps)
        
        xty.block.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        xtx.block.values <- matrix(data=NA,nrow =((n.annot+1)* n.blocks),ncol =(n.annot+1))
        colnames(xty.block.values)<- colnames(xtx.block.values)<- colnames(weighted.LD)
        replace.from <- seq(from=1,to=nrow(xtx.block.values),by =(n.annot+1))
        replace.to <- seq(from =(n.annot+1),to=nrow(xtx.block.values),by =(n.annot+1))
        for(i in 1:n.blocks){
          xty.block.values[i,] <- t(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.chi[select.from[i]:select.to[i],])
          xtx.block.values[replace.from[i]:replace.to[i],] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],])%*% weighted.LD[select.from[i]:select.to[i],])
        }
        xty <- as.matrix(colSums(xty.block.values))
        xtx <- matrix(data=NA,nrow =(n.annot+1),ncol =(n.annot+1))
        colnames(xtx)<- colnames(weighted.LD)
        for(i in 1:nrow(xtx)){xtx[i,] <- t(colSums(xtx.block.values[seq(from=i,to=nrow(xtx.block.values),by=ncol(weighted.LD)),]))}
        
        reg <- solve(xtx, xty)
        intercept <- reg[2]
        coefs <- reg[1]/N.bar
        reg.tot <- coefs*m
        
        delete.from <- seq(from=1,to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.to <- seq(from=ncol(xtx.block.values),to=nrow(xtx.block.values),by=ncol(xtx.block.values))
        delete.values <- matrix(data=NA,nrow=n.blocks,ncol =(n.annot+1))
        colnames(delete.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){
          xty.delete <- xty-xty.block.values[i,]
          xtx.delete <- xtx-xtx.block.values[delete.from[i]:delete.to[i],]
          delete.values[i,] <- solve(xtx.delete, xty.delete)
        }
        
        tot.delete.values <- delete.values[,1:n.annot]
        pseudo.values <- matrix(data=NA,nrow=n.blocks,ncol=length(reg))
        colnames(pseudo.values)<- colnames(weighted.LD)
        for(i in 1:n.blocks){pseudo.values[i,] <- (n.blocks*reg)-((n.blocks-1)* delete.values[i,])}
        
        jackknife.cov <- cov(pseudo.values)/n.blocks
        jackknife.se <- sqrt(diag(jackknife.cov))
        intercept.se <- jackknife.se[length(jackknife.se)]
        coef.cov <- jackknife.cov[1:n.annot,1:n.annot]/(N.bar^2)
        cat.cov <- coef.cov*(m %*% t(m))
        tot.cov <- sum(cat.cov)
        tot.se <- sqrt(tot.cov)
        
        V.hold[, s] <- pseudo.values[, 1]
        N.vec[1, s] <- N.bar
        
        V.delete[, s] = delete.values[,1]
        V.names[s] = paste(trait.names[j], trait.names[k], sep = "-")
        
        cov[k, j] <- cov[j, k] <- reg.tot
        I[k, j] <- I[j, k] <- intercept
        
        .LOG("Results for genetic covariance between: ", trait.names[j], " and ", trait.names[k], file=log.file)
        .LOG("Mean Z*Z: ", round(mean(merged$ZZ), 4), file=log.file)
        .LOG("Cross trait Intercept: ", round(intercept, 4), " (", round(intercept.se, 4), ")", file=log.file)
        .LOG("Total Observed Scale Genetic Covariance (g_cov): ", round(reg.tot, 4), " (", round(tot.se, 4), ")", file=log.file)
        .LOG("g_cov Z: ", format(reg.tot / tot.se, digits = 3), file=log.file)
        .LOG("g_cov P-value: ", format(2 * stats::pnorm(abs(reg.tot / tot.se), lower.tail = FALSE), digits = 5), file=log.file)
      }
      
      ### Total count
      s <- s + 1
    }
  }
  
  

  #rescale the sampling correlation matrix by the appropriate diagonals
  # V <- v.out * tcrossprod(scaleO)
  V = cov(V.hold) / crossprod(N.vec * (sqrt(n.blocks) / m))
  
  
  # name traits according to trait.names argument
  colnames(cov) <- trait.names
  
  # from covariance to genetic correlation
  if(all(diag(cov) > 0)){ # use cov instead of S + add variance of denominator!!
    
    ##calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(cov))) 
    S_Stand <- cov * ratio
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    
    
    # use JK
    ss = 1
    for(j in 1:n.traits){
      
      for(k in j:n.traits){
        # pseudo values needed to get pseudo values for genetic correlation
        vjj = which(V.names==paste(trait.names[j], trait.names[j], sep = "-"))
        vkk = which(V.names==paste(trait.names[k], trait.names[k], sep = "-"))
        vjk = which(V.names==paste(trait.names[j], trait.names[k], sep = "-"))
        # rg = rho_g / sqrt(h2_1 * h2_2)
        V_Stand.delete[,ss] = V.delete[,vjk]/N.vec[vjk]*m / sqrt(V.delete[,vjj]/N.vec[vjj]*m * V.delete[,vkk]/N.vec[vkk]*m)  
        V_Stand.hold[,ss] = n.blocks*S_Stand[j,k] - ((n.blocks-1)* V_Stand.delete[,ss])
        
        
        ss = ss+1
      }
    }
    V_Stand = cov(V_Stand.hold)/((n.blocks))
    
    
    #enter SEs from diagonal of standardized V
    r<-nrow(cov) # NM # use cov instead of S
    SE_Stand<-matrix(0, r, r)
    SE_Stand[lower.tri(SE_Stand,diag=TRUE)] <-sqrt(diag(V_Stand))
    
    
    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Genetic Correlation Results", file=log.file)
    
    for(j in 1:n.traits){
      if(is.null(trait.names)){
        chi1<-traits[j]
      }else{chi1 <- trait.names[j]}
      for(k in j:length(traits)){
        if(j != k){
          if(is.null(trait.names)){
            chi2<-traits[k]
          }else{chi2 <- trait.names[k]}
          .LOG("Genetic Correlation between ", chi1, " and ", chi2, ": ",
               round(S_Stand[k, j], 4), " (", round(SE_Stand[k, j], 4), ")", file=log.file)
        }
      }
    }
  }else{
    warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
    .LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.", file=log.file, print = FALSE)
    .LOG("Genetic correlation results could not be computed due to negative heritability estimates.", file=log.file)
  }
  
  #### NM modification starts ####
  
  # get partial covariances (and partial h2) (GW)
  for(j in 1:n.conditions){
    for(k in j:n.conditions){
      partial.cov[k, j] <- partial.cov[j, k] <- (cov[j,k] - 
                                                   t(cov[j, (n.conditions+1):n.traits]) %*% 
                                                   solve((cov[(n.conditions+1):n.traits, (n.conditions+1):n.traits])) %*% 
                                                   t(t(cov[k, (n.conditions+1):n.traits])))
      
    }
  }
  
  # get their variances (using JK results, stored in V.hold)
  # for simplicity, use V.names, since V.hold has a different format
  ss = 1
  vCC = matrix(NA, nrow=n.confounders, ncol=n.confounders)
  colnames(vCC) = confounder.names
  rownames(vCC) = confounder.names
  for(row in 1:n.confounders){
    for(col in row:n.confounders){
      vCC[row,col] = vCC[col, row] = which(V.names==paste(rownames(vCC)[row], colnames(vCC)[col], sep = "-"))
    }
  }
  
  for(j in 1:n.conditions){
    
    for(k in j:n.conditions){
      # pseudo values for partial cov
      vjk = which(V.names==paste(trait.names[j], trait.names[k], sep = "-"))
      vjC = matrix(NA, nrow=1, ncol=n.confounders)
      colnames(vjC) = confounder.names
      rownames(vjC) = trait.names[j]
      for(col in 1:n.confounders){
        vjC[1,col] = which(V.names==paste(rownames(vjC)[1], colnames(vjC)[col], sep = "-"))
      }
      vkC = matrix(NA, nrow=n.confounders, ncol=1)
      colnames(vkC) =  trait.names[k]
      rownames(vkC) = confounder.names
      for(row in 1:n.confounders){
        vkC[row,1] = which(V.names==paste(colnames(vkC)[1], rownames(vkC)[row], sep = "-"))
      }
      
      # a bit ugly to use a loop here, but needed to implement something quickly and 
      # didn't have time to look into more elegant options
      for(nb in 1:n.blocks){
        partial.V.delete[nb,ss] =  V.delete[nb,vjk]/N.vec[vjk]*m - 
          (V.delete[nb,vjC]/N.vec[as.numeric(vjC)]*m) %*% 
          solve(matrix((V.delete[nb,vCC]/N.vec[as.numeric(vCC)]*m), nrow = n.confounders)) %*% 
          (V.delete[nb,vkC]/N.vec[as.numeric(vkC)]*m)
        
        partial.V.hold[nb,ss] = n.blocks*partial.cov[j,k] - ((n.blocks-1)* partial.V.delete[nb,ss])
        partial.V.names[ss] = paste(trait.names[j], trait.names[k], sep = "-")
        
      }  
      ss = ss+1
    }
  }
  
  # variance of the partial genetic covariance / h2 estimates
  partial.V <- cov(partial.V.hold) /  ((n.blocks))
  
  # rescale to correlation
  if(all(diag(partial.cov) > 0)){
    
    ##calculate standardized results to print genetic correlations to log and screen
    ratio <- tcrossprod(1 / sqrt(diag(partial.cov)))
    partial.S_Stand <- partial.cov * ratio
    
    #calculate the ratio of the rescaled and original S matrices
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    
    ## MAke sure that if ratio in NaN (devision by zero) we put the zero back in
    # -> not possible because of 'all(diag(S) > 0)'
    # scaleO[is.nan(scaleO)] <- 0
    
    #rescale the sampling correlation matrix by the appropriate diagonals
    # partial.V_Stand <- partial.V * tcrossprod(scaleO)
    
    # use JK
    
    ss = 1
    for(j in 1:n.conditions){
      
      for(k in j:n.conditions){
        # pseudo values needed to get pseudo values for genetic correlation
        vjj = which(partial.V.names==paste(trait.names[j], trait.names[j], sep = "-"))
        vkk = which(partial.V.names==paste(trait.names[k], trait.names[k], sep = "-"))
        vjk = which(partial.V.names==paste(trait.names[j], trait.names[k], sep = "-"))
        # rg = rho_g / sqrt(h2_1 * h2_2)
        partial.V_Stand.delete[,ss] = partial.V.delete[,vjk] / sqrt( partial.V.delete[,vjj] *  partial.V.delete[,vkk])  
        partial.V_Stand.hold[,ss] = n.blocks*partial.S_Stand[j,k] - ((n.blocks-1)* partial.V_Stand.delete[,ss])
        
        
        ss = ss+1
      }
    }
    partial.V_Stand = cov(partial.V_Stand.hold)/((n.blocks))
    
    
    #enter SEs from diagonal of standardized V
    r<-nrow(partial.cov)
    partial.SE_Stand<-matrix(0, r, r)
    partial.SE_Stand[lower.tri(partial.SE_Stand,diag=TRUE)] <-sqrt(diag(partial.V_Stand))
    
    
    .LOG(c("     ", "     "), file=log.file, print = FALSE)
    .LOG("Partial genetic Correlation Results", file=log.file)
    
    for(j in 1:n.conditions){
      if(is.null(condition.names)){
        chi1<-traits[j]
      }else{chi1 <- condition.names[j]}
      for(k in j:length(conditions)){
        if(j != k){
          if(is.null(condition.names)){
            chi2<-traits[k]
          }else{chi2 <- condition.names[k]}
          .LOG("Partial genetic Correlation between ", chi1, " and ", chi2, ": ",
               round(partial.S_Stand[k, j], 4), " (", round(partial.SE_Stand[k, j], 4), ")", file=log.file)
        }
      }
    }
    
    
    
    
  }
  
  
  # covariance between raw/partial genetic covariances / correlations
  for(j in 1:n.conditions){
    
    for(k in j:n.conditions){
      
      # pseudo values for raw rg
      vraw = which(V.names==paste(trait.names[j], trait.names[k], sep = "-"))
      # pseudo values for partial rg
      vpartial = which(partial.V.names==paste(trait.names[j], trait.names[k], sep = "-"))
      
      Vcov_Stand[j,k] <- Vcov_Stand[k,j] <- cov(V_Stand.hold[,vraw], partial.V_Stand.hold[,vpartial])/((n.blocks))
      
      
    }
  }
  row.names(cov) = trait.names
  row.names(S_Stand) = trait.names
  row.names(partial.cov) = condition.names
  colnames(partial.cov) = condition.names
  row.names(partial.S_Stand) = condition.names
  colnames(partial.S_Stand) = condition.names
  row.names(I) = trait.names
  colnames(I) = trait.names
  
  # make a nice table with test statistic for difference & p-values
  res = data.frame(condition.1 = NA_character_,
                   condition.2 = NA_character_,
                   rg = NA_real_,
                   rg.SE = NA_real_,
                   partial_rg = NA_real_,
                   partial_rg.SE = NA_real_,
                   rg_cov = NA_real_)
  
  s = 1
  for(j in 1:n.conditions){
    
    for(k in j:n.conditions){
      if(j != k){
        res[s,1:2] = c(condition.names[j], condition.names[k])
        
        # index for raw corr
        vraw = which(V.names==paste(trait.names[j], trait.names[k], sep = "-"))
        # index for partial corr
        vpartial = which(partial.V.names==paste(trait.names[j], trait.names[k], sep = "-"))
        
        
        res[s,3] = S_Stand[j,k]
        res[s,4] = SE_Stand[k,j]
        
        res[s,5] = partial.S_Stand[j,k]
        res[s,6] = partial.SE_Stand[k,j]
        
        res[s,7] = Vcov_Stand[j,k]
        
        
        s = s+1
      }
    }
  }
  
  res %>%
    dplyr::mutate(diff.T = (.data$rg-.data$partial_rg)/sqrt(.data$rg.SE**2 + .data$partial_rg.SE**2- 2*.data$rg_cov),
           diff.P = 2*stats::pnorm(-abs(.data$diff.T))) -> res
  
  
  end.time <- Sys.time()
  end.time.nice <- paste0(lubridate::hour(end.time), ":", paste0(ifelse(lubridate::minute(end.time)<10, "0", ""), lubridate::minute(end.time)))
  
  total.time <- difftime(time1=end.time,time2=begin.time,units="sec")
  mins <- floor(floor(total.time)/60)
  secs <- floor(total.time-mins*60)
  
  .LOG("     ", file=log.file, print = FALSE)
  .LOG("Analysis finished running at ", end.time.nice, file=log.file)
  .LOG("Runtime: ", mins, " minute(s) and ", secs, " second(s)", file=log.file)
  .LOG("     ", file=log.file, print = FALSE)
  
  flush(log.file)
  close(log.file)
  

  readr::write_tsv(res, paste0(log.name, "_difference.tsv"))
  
  
  # return everything 
  list(res_diff = res,
       S = cov, V = V,
       S_Stand = S_Stand, V_Stand = V_Stand,
       partial.S = partial.cov, partial.V = partial.V,
       partial.S_Stand = partial.S_Stand, partial.V_Stand = partial.V_Stand,
       I=I)
  
  
  
}