#' Predict subtypes, pairwise and by index
#'
#' A function that tests, for each sample, each kTSP model in the family,
#' and each possible pair of subtypes, which subtype is the best choice for the sample
#'
#' @description
#' This is the most basic step in the process of predicting a subtype for a
#' sample. The BCCSclassifier multi-kTSP object is a list of kTSP models, each
#' consisting of a list of kTSP tests for the best subtype between a pair of possible
#' choices, ranging over all kTSP models in the family. This
#' function applies ktspair prediction to each of these inner models and arranges
#' the results in a tibble, respecting the kTSP model index.
#'
#'
#'
#' @param dat An object of class matrix with samples as column names and gene symbols
#'   as row names, and entries the corresponding expression levels. The included tibble
#'   \code{bccs_classifier_genes} contains the gene symbols and ENTREZIDs used in the
#'   BCCSclassifier models, along with Affymetrix (hgu133a) and Illumina
#'   (illuminaHuman v3) probes that can be used to represent the associated genes.
#'   Expression values must \bold{not} be scaled or median-centered.
#'
#' @param model A character string that is one of "erposneg", "erpos", "erneg",
#'   to determine whether BCCSclassifier(ER+/-), BCCSclassifier(ER+), or
#'   BCCSclassifier(ER-), respectively, should be used for subtyping.
#'
#' @return The output is a tibble (data.frame),with columns sample, pair,
#'   pair_prediction, index. pair indicates the pair of subtypes tested in the
#'  kTSP inner model. pair_prediction is the subtype in the pair chosen
#'   by the kTSP model. index (an integer) is the index of the kTSP model in the
#'   Multi-kTSP family. This object is called a PPI object, short for
#'   \emph{pairwise prediction - indexed}
#'
#'
#' @export
make_pairwise_predictions <- function(dat, model) {
  MP1 <- identify_mpred_by_model(model)
  MP <- filter_mpred_by_genes(MP1, dat)
  print(paste0(length(MP), " out of ", length(MP1), " kTSP models used."))
  PPI <- purrr::map(
    1:length(MP),
    ~ predict_for_pair(MP[[.]], B=., dat = dat, model)
  ) %>% dplyr::bind_rows()
  PPI
}

#' Identifies the Multi-kTSP model object corresponding to the \code{model}
#' parameter
#'
#' @param model A character string that is one of "erposneg", "erpos", "erneg",
#'   to determine whether BCCSclassifier(ER+/-), BCCSclassifier(ER+), or
#'   BCCSclassifier(ER-), respectively, should be used for subtyping.
#'
#' @return The multi-kTSP model used for prediction (termed an MPred object for short)
#'   corresponding to the input \code{model}
#'
#' @export
identify_mpred_by_model <- function(model) {
  if(! model %in% c("erposneg", "erpos", "erneg")) {
    stop("model must be one of 'erposneg', 'erpos', or 'erneg'.")
  }
  MP <- NULL
  if(model == "erposneg") {MP <- MPred_erposneg}
  if(model == "erpos") {MP <- MPred_erpos}
  if(model == "erneg") {MP <- MPred_erneg}
  MP
}

# --- Helper functions, not for independent use -----

predict_for_pair <- function(M, B, dat, model) {
  ppred <- vector(mode = "list", length = length(M))
  for (k in 1:length(M))
  {
    p <- suppressWarnings(predict_ktsp(M[[k]], dat, display = F)) # ktspair prediction
    if(is.null(p)) {
      ppred[[k]] <- rep(NA, times = ncol(dat))
    } else {
      ppred[[k]] <- translation_vector(model)[p]
    }
    names(ppred[[k]]) <- NULL
  }
  names(ppred) <- names(M)
  pred_df1 <- dplyr::bind_cols(ppred)
  pred_df1 <- pred_df1 %>% dplyr::mutate(sample = colnames(dat))
  pred_df2 <- pred_df1 %>%
    tidyr::pivot_longer(cols = !sample, names_to = "pair", values_to = "pairwise_prediction")
  pred_df2 <- pred_df2 %>% dplyr::mutate(index = B)
  pred_df2
}


# This translates the subtype labels in kTSP models from the original ones
# to those we used in the publication

translation_vector <- function(model) {
  erposneg_trans <- paste0("BCS", 1:5)
  names(erposneg_trans) <- c("C2", "C3", "C1", "C6", "C5")
  erneg_trans <- paste0("NCS", 1:3)
  names(erneg_trans) <- c("C2", "C6", "C3")
  erpos_trans <- paste0("PCS", 1:4)
  names(erpos_trans) <- paste0("PT", 1:4)
  trans_vector <- list(
    erposneg = erposneg_trans,
    erneg = erneg_trans,
    erpos = erpos_trans
  )
  trans_vector[[model]]
}

# Following is a direct copy of the predict.ktsp function from the
# ktspair package, created by Julien Damond <julien.damond@gmail.com>

predict_ktsp <- function(object, dat = NULL, select = NULL, display = TRUE,...){
  ## Compute prediction of the dataset on which the model is based or of a new dataset if one is inserted.

  ktspobj <- object
  grp <- ktspobj$grp
  k <- ktspobj$k
  grplabels <- character(length(grp))
  grplabels[grp==0] <- ktspobj$labels[1]
  grplabels[grp==1] <- ktspobj$labels[2]
  predict <- c()

  if(!is.null(dat) && ktspobj$med == TRUE){
    dat[c(ktspobj$index),] <- dat[c(ktspobj$index),]-ktspobj$med
  }

  if(is.vector(dat)){
    dat <- as.matrix(dat,length(dat),1)
  }

  if(!is.null(select) && (select <1 || select > k)){stop("The selected pair of genes is not available.")}

  if(is.null(dat) & !is.null(select)){
    z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
    z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
    table_label <- dimnames(table(z1,grplabels))[[2]]

    p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
    p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

    p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
    p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

    max <- which.max(c(abs(p11-p21), abs(p12-p22)))

    if(max==1){
      predict[which(z1 == 0)] <- table_label[which.max(table(z1,grplabels)[1,])]
      predict[which(z1 == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
    }

    if(max==2){
      predict[which(z2 == 0)] <- table_label[which.max(table(z2,grplabels)[1,])]
      predict[which(z2 == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
    }
    return(predict)
  }
  if(is.null(select) && is.null(dat)){

    if(k%%2 == 0){
      k2 <- k-1
      cat("The prediction based on the majority voting procedure cannot be computed for an even value of k. \n")
      cat("The value of k has been reduced to k = ", k2, ". \n")
      k <- k2
    }

    predict <- character(length(grp))
    vote <- numeric(length(grp))
    vote2 <-matrix(nrow=length(grp), ncol=k)
    count <- numeric(length(grp))
    for(i in 1:k){
      z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
      z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

      p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
      p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

      p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
      p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

      max <- which.max(c(abs(p11-p21), abs(p12-p22)))

      if(max==1){
        vote2[which(z1 == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1
        vote2[which(z1 == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1
        count[is.na(z1)==FALSE] <- count[is.na(z1)==FALSE]+1
      }

      if(max==2){
        vote2[which(z2 == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1
        vote2[which(z2 == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1
        count[is.na(z2)==FALSE] <- count[is.na(z2)==FALSE]+1
      }
    }

    table_label <- labels(table(z1,grplabels))[2]
    vote <- apply(vote2, 1, sum)
    predict[(vote/count < 1/2)] <- table_label$grplabels[1]
    predict[(vote/count > 1/2)] <- table_label$grplabels[2]
    nopred <- which((predict=="")==TRUE)
    if(length(nopred)>0){predict2 <- character(length(nopred))
    if(display==TRUE){
      cat("For the observation(s): ", nopred, " a part of the data is missing \n")
      cat("Their prediction was computed on a subset of the k-TSP\n")
    }
    for(i in 1:length(nopred)){
      vote3 <- numeric(length(nopred))
      vote3 <- vote2[i,is.na(vote2[i,])==FALSE]
      pred <- mean(vote3)
      if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
      if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
      else{
        pred <- mean(vote3[-length(vote3)])
        if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
        if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
      }
    }
    predict[nopred] <- predict2
    }

    return(predict)
  }

  if(!is.null(dat) & !is.null(select)){
    ktspnames <- rownames(ktspobj$ktspdat)[c(select,(select + k))]
    predict <- character(dim(dat)[2])
    vote <- numeric(length(grp))
    vote2 <-matrix(nrow=length(grp), ncol=k)
    count <- numeric(length(grp))
    if("matrix" %in% class(dat)){
      if(is.null(rownames(dat))){
        if(display==TRUE){
          cat("No rownames found, using indices \n")}
        ktspnames <- as.numeric(ktspnames)
        if(any(!(ktspnames %in% 1:dim(dat)[1]))){stop("Rownames of new data do not include the TSP names")}
      }
      else{
        if(any(!(ktspnames %in% rownames(dat)))){stop("Rownames of new data do not include the TSP names")}
      }

      z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
      z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
      table_label <- as.vector(dimnames(table(z1,grplabels))[[2]])

      p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
      p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

      p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
      p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

      max <- which.max(c(abs(p11-p21), abs(p12-p22)))

      if(max==1){
        w <- dat[ktspnames[1],] < dat[ktspnames[2],]
        predict[which(w == 0)] <- table_label[which.max(table(z1,grplabels)[1,])]
        predict[which(w == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
      }

      if(max==2){
        w <- dat[ktspnames[1],] > dat[ktspnames[2],]
        predict[which(w == 0)] <- table_label[which.max(table(z2,grplabels)[1,])]
        predict[which(w == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
      }
      nopred <- which((predict=="")==TRUE)

      if(length(nopred)>0){
        if(display==TRUE){
          cat("For the observation(s): ", nopred, " a part of the data is missing \n")
          cat("Their prediction for the selected pair is not possible\n")}
      }
      return(predict)
    }
    if("ExpressionSet" %in% class(dat)){
      if(is.null(featureNames(dat))){
        if(display==TRUE){
          cat("No featureNames info found, using indices \n")
        }
        ktspnames <- as.numeric(ktspnames)
        if(any(!(ktspnames %in% 1:dim(exprs(dat))[1]))){stop("Rownames of new data do not include the TSP names")}
        dat <- exprs(dat)
      }
      else{
        if(any(!(ktspnames %in% featureNames(dat)))){stop("Rownames of new data do not include the TSP names")}
        genenames <- featureNames(dat)
        dat <- exprs(dat)
        rownames(dat) <- genenames
      }
      z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
      z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
      table_label <- dimnames(table(z1,grplabels))[[2]]

      z1 <- ktspobj$ktspdat[select,] < ktspobj$ktspdat[(select + k),]
      z2 <- ktspobj$ktspdat[select,] > ktspobj$ktspdat[(select + k),]
      table_label <- as.vector(dimnames(table(z1,grplabels))[[2]])

      p11 <- mean(ktspobj$ktspdat[select,which(grp==0)]<ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
      p21 <- mean(ktspobj$ktspdat[select,which(grp==1)]<ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

      p12 <- mean(ktspobj$ktspdat[select,which(grp==0)]>ktspobj$ktspdat[(select+k),which(grp==0)], na.rm=TRUE)
      p22 <- mean(ktspobj$ktspdat[select,which(grp==1)]>ktspobj$ktspdat[(select+k),which(grp==1)], na.rm=TRUE)

      max <- which.max(c(abs(p11-p21), abs(p12-p22)))

      if(max==1){
        w <- dat[ktspnames[1],] < dat[ktspnames[2],]
        predict[which(w == 0)] <- table_label[which.max(table(z1,grplabels)[1,])]
        predict[which(w == 1)] <- table_label[which.max(table(z1,grplabels)[2,])]
      }

      if(max==2){
        w <- dat[ktspnames[1],] > dat[ktspnames[2],]
        predict[which(w == 0)] <- table_label[which.max(table(z2,grplabels)[1,])]
        predict[which(w == 1)] <- table_label[which.max(table(z2,grplabels)[2,])]
      }
      nopred <- which((predict=="")==TRUE)
      if(length(nopred)>0){
        if(display==TRUE){
          cat("For the observation(s): ", nopred, " a part of the data is missing \n")
          cat("Their prediction for the selected pair is not possible\n")
        }
      }
      return(predict)
    }
  }

  else{

    if(k%%2 == 0){
      k2 <- k-1
      if(display==TRUE){
        cat("The prediction based on the k-TSP cannot be computed for an even value of k. \n")
        cat("The value of k has been reduced to k = ", k2, ". \n")
      }
      k <- k2
    }
    ktspnames <- rownames(ktspobj$ktspdat)
    if("matrix" %in% class(dat)){
      if(is.null(rownames(dat))){
        if(display==TRUE){
          cat("No rownames found, using indices \n")
        }
        ktspnames <- as.numeric(ktspnames)
        if(any(!(ktspnames %in% 1:dim(dat)[1]))){stop("Rownames of new data do not include the TSP names")}
      }

      else{
        if(any(!(ktspnames %in% rownames(dat)))){stop("Rownames of new data do not include the TSP names")}
      }
      predict <- character(dim(dat)[2])
      vote <- numeric(dim(dat)[2])
      vote2 <-matrix(nrow=dim(dat)[2], ncol=k)
      count <- numeric(dim(dat)[2])
      for(i in 1:k){
        z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
        z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

        p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
        p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

        p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
        p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

        max <- which.max(c(abs(p11-p21), abs(p12-p22)))

        if(max==1){
          w <- dat[ktspnames[i],] < dat[ktspnames[i+k],]
          vote2[which(w == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1
          vote2[which(w == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1
          count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
        }

        if(max==2){
          w <- dat[ktspnames[i],] > dat[ktspnames[i+k],]
          vote2[which(w == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1
          vote2[which(w == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1
          count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
        }
      }
      table_label <- labels(table(z1,grplabels))[2]
      vote <- apply(vote2, 1, sum)
      predict[(vote/count < 1/2)] <- table_label$grplabels[1]
      predict[(vote/count > 1/2)] <- table_label$grplabels[2]
      nopred <- which((predict=="")==TRUE)

      #If predictions were not possible for some observations, the prediction will be based on a restricted list of the pairs in the k-TSP.

      if(length(nopred)>0){
        predict2 <- character(length(nopred))
        if(display==TRUE){
          cat("For the observation(s): ", nopred, " a part of the data is missing \n")
          cat("Their prediction was computed on a subset of the k-TSP\n")
        }
        for(i in 1:length(nopred)){
          vote3 <- numeric(length(nopred))
          vote3 <- vote2[nopred[i],is.na(vote2[nopred[i],])==FALSE]
          if(length(vote3)<1){if(display==TRUE){cat("Prediction not possible for observation", nopred[i], " \n")}}
          else{
            pred <- mean(vote3)
            if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
            if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
            else{
              if(length(vote3)<2){if(display==TRUE){cat("Prediction not possible for observation", nopred[i], "\n")}}
              else{
                pred <- mean(vote3[-length(vote3)])
                if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
                if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
              }
            }
          }
        }
        predict[nopred] <- predict2
      }
      return(predict)
    }

    if("ExpressionSet" %in% class(dat)){
      if(is.null(featureNames(dat))){
        if(display==TRUE){
          cat("No featureNames info found, using indices \n")}
        ktspnames <- as.numeric(ktspnames)
        if(any(!(ktspnames %in% 1:dim(exprs(dat))[1]))){stop("Rownames of new data do not include the TSP names")}
        dat <- exprs(dat)
      }

      else{
        if(any(!(ktspnames %in% featureNames(dat)))){stop("Rownames of new data do not include the TSP names")}
        genenames <- featureNames(dat)
        dat <- exprs(dat)
        rownames(dat) <- genenames
      }
      predict <- character(dim(dat)[2])
      vote <- numeric(dim(dat)[2])
      vote2 <-matrix(nrow=dim(dat)[2], ncol=k)
      count <- numeric(dim(dat)[2])
      for(i in 1:k){
        z1 <- ktspobj$ktspdat[i,] < ktspobj$ktspdat[(i + k),]
        z2 <- ktspobj$ktspdat[i,] > ktspobj$ktspdat[(i + k),]

        p11 <- mean(ktspobj$ktspdat[i,which(grp==0)]<ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
        p21 <- mean(ktspobj$ktspdat[i,which(grp==1)]<ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

        p12 <- mean(ktspobj$ktspdat[i,which(grp==0)]>ktspobj$ktspdat[(i+k),which(grp==0)], na.rm=TRUE)
        p22 <- mean(ktspobj$ktspdat[i,which(grp==1)]>ktspobj$ktspdat[(i+k),which(grp==1)], na.rm=TRUE)

        max <- which.max(c(abs(p11-p21), abs(p12-p22)))

        if(max==1){
          w <- dat[ktspnames[i],] < dat[ktspnames[i+k],]
          vote2[which(w == 0),i] <- which.max(table(z1,grplabels)[1,]) - 1
          vote2[which(w == 1),i] <- which.max(table(z1,grplabels)[2,]) - 1
          count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
        }

        if(max==2){
          w <- dat[ktspnames[i],] > dat[ktspnames[i+k],]
          vote2[which(w == 0),i] <- which.max(table(z2,grplabels)[1,]) - 1
          vote2[which(w == 1),i] <- which.max(table(z2,grplabels)[2,]) - 1
          count[is.na(w)==FALSE] <- count[is.na(w)==FALSE]+1
        }
      }
      table_label <- labels(table(z1,grplabels))[2]
      vote <- apply(vote2, 1, sum)
      predict[(vote/count < 1/2)] <- table_label$grplabels[1]
      predict[(vote/count > 1/2)] <- table_label$grplabels[2]
      nopred <- which((predict=="")==TRUE)
      if(length(nopred)>0){
        predict2 <- character(length(nopred))
        if(display==TRUE){
          cat("For the observation(s): ", nopred, " a part of the data is missing \n")
          cat("Their prediction was computed on a subset of the k-TSP\n")
        }
        for(i in 1:length(nopred)){
          vote3 <- numeric(length(nopred))
          vote3 <- vote2[nopred[i],is.na(vote2[nopred[i],])==FALSE]
          if((length(vote3)<1) & (display == TRUE)){cat("Prediction not possible for observation", nopred[i], " \n")}
          else{
            pred <- mean(vote3)
            if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
            if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
            else{
              if((length(vote3)<2) & (display == TRUE)){cat("Prediction not possible for observation", nopred[i], "\n")}
              else{
                pred <- mean(vote3[-length(vote3)])
                if(pred < 1/2){predict2[i] <- table_label$grplabels[1]}
                if(pred > 1/2){predict2[i] <- table_label$grplabels[2]}
              }
            }
          }
        }
        predict[nopred] <- predict2
      }
      return(predict)
    }
  }
}



