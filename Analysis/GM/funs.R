## =============================================================
## GM (Gaussian Model) pipeline functions
## =============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(mixtools)
library(rmarkdown)
library(readr)
library(rlang)

#' getspliceMax
#'
#' @param dt variant data frame with "SpliceAI_DS_AG", "SpliceAI_DS_AL", "SpliceAI_DS_DG", "SpliceAI_DS_DL" or "SpliceMax"
#' @param splice_cutoff threshold to define splice region for spliceAI values, default is 0.2
#'
#' @return dt with spliceM
#' @export
#'
#' @examples getspliceMax(dt, splice_cutoff = .2)
getspliceMax <- function(dt, splice_cutoff = .2) {
  ## input raw counts in long format with annotation
  ## required column names SpliceAI_DS_AG, SpliceAI_DS_AL, SpliceAI_DS_DG, SpliceAI_DS_DL
  ## return data with SpliceMax and spliceR columns
  if ("SpliceMax" %in% colnames(dt)) {
    cat('NOTE: the input file does have SpliceMax, we will use that directly\n')
    dt <- dt %>% mutate(spliceR = ifelse(SpliceMax > .2, TRUE, FALSE))
  } else {
    if (
      all(
        c(
          "SpliceAI_DS_AG",
          "SpliceAI_DS_AL",
          "SpliceAI_DS_DG",
          "SpliceAI_DS_DL"
        ) %in%
        colnames(dt)
      )
    ) {
      dt$SpliceAI_DS_AG <- as.numeric(dt$SpliceAI_DS_AG)
      dt$SpliceAI_DS_AL <- as.numeric(dt$SpliceAI_DS_AL)
      dt$SpliceAI_DS_DG <- as.numeric(dt$SpliceAI_DS_DG)
      dt$SpliceAI_DS_DL <- as.numeric(dt$SpliceAI_DS_DL)
      
      dt <- dt %>%
        mutate(
          SpliceMax = pmax(
            SpliceAI_DS_AG,
            SpliceAI_DS_AL,
            SpliceAI_DS_DG,
            SpliceAI_DS_DL,
            na.rm = T
          )
        ) %>%
        mutate(spliceR = ifelse(SpliceMax > .2, TRUE, FALSE))
    } else {
      cat(
        'NOTE: the input file does not have these required columns SpliceAI_DS_AG, SpliceAI_DS_AL, SpliceAI_DS_DG, SpliceAI_DS_DL\n'
      )
      next
    }
  }
  return(dt)
}


#' normfreq
#'
#' @param dt raw counts data frame in long format
#'
#' @return freq data frame in wide format with normalized frequency which is either natural log or cubic root or log 2 ratios
#' @export
#'
#' @examples normfreq(dt)
normfreq <- function(dt) {
  ## input raw counts in long format
  ## logbase could be n (natural) or 2
  ## output long format with variant each row, time point freq
  ## dt one exon one replicate
  cat("...Running normalizing freq... \n")
  if (!"use4norm" %in% colnames(dt)) {
    dt$use4norm <- 1
  }
  
  wide <- dt %>%
    arrange(POS, Time) %>%
    tidyr::pivot_wider(
      id_cols = c(
        CHROM,
        POS,
        REF,
        ALT,
        EventType,
        Exon,
        Rep,
        use4norm,
        spliceR
      ),
      names_from = c(Time),
      values_from = EventCount,
      names_glue = "rawCount_T{.name}"
    )
  
  ## create freq columns
  ## create log base ratio
  df <- wide %>%
    group_by(Exon, Rep) %>%
    mutate(
      across(
        .cols = matches("^rawCount_T\\d+"),
        .fns = ~ log(
          (.x + 0.5) / (sum(.x[use4norm == 1], na.rm = TRUE) + 0.5)
        ),
        .names = "log_freq_{str_remove(.col, 'rawCount_')}"
      )
    ) %>%
    ungroup()
  ## create cubic root ratio
  df <- df %>%
    group_by(Exon, Rep) %>%
    mutate(
      across(
        .cols = matches("^rawCount_T\\d+"),
        .fns = ~ (
          (.x + 0.5) / (sum(.x[use4norm == 1], na.rm = TRUE) + 0.5)
        )^(1/3),
        .names = "cr_freq_{str_remove(.col, 'rawCount_')}"
      )
    ) %>%
    ungroup()
  ## create log 2 ratio
  df <- df %>%
    group_by(Exon, Rep) %>%
    mutate(
      across(
        .cols = matches("^rawCount_T\\d+"),
        .fns = ~ log2(
          (.x + 0.5) / (sum(.x[use4norm == 1], na.rm = TRUE) + 0.5)
        ),
        .names = "log2_freq_{str_remove(.col, 'rawCount_')}"
      )
    ) %>%
    ungroup()
  
  return(df)
}


#' getratio
#'
#' @param df freq data frame with log_freq, log2_freq, cr_freq at time 0, 1, or 2
#'
#' @return ratio data frame with additional columns for ratio_T2_T1 (ratio between time 2 vs 1), ratio_T2_T0 (ratio between time 2 vs 0), ratio_T1_T0 (ratio between time 1 vs 0)
#' @export
#'
#' @examples getratio(dt)
getratio <- function(df, transf_pref = "n") {
  ## input : log_freq at time 0, 1, or 2, wide format
  ## output: ratio between times, if no tim0, only time2-1 ratio is calculated
  ##        for each rep
  cat("...Running getting ratio between times... \n")
  df <- df %>% mutate(uPOS = paste(POS, REF, ALT, sep = "_"))
  
  if(transf_pref == "n"){
    t0_col <- ("log_freq_T0")
    t1_col <- ("log_freq_T1")
    t2_col <- ("log_freq_T2")
    
    # T2_T1
    if (all(c(t2_col, t1_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T2_T1") := .data[[t2_col]] - .data[[t1_col]])
    }
    # T2_T0
    if (all(c(t2_col, t0_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T2_T0") := .data[[t2_col]] - .data[[t0_col]])
    }
    # T1_T0
    if (all(c(t1_col, t0_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T1_T0") := .data[[t1_col]] - .data[[t0_col]])
    }
  }else{ ## cubic root
    t0_col <- ("cr_freq_T0")
    t1_col <- ("cr_freq_T1")
    t2_col <- ("cr_freq_T2")
    
    # T2_T1
    if (all(c(t2_col, t1_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T2_T1") := .data[[t2_col]] / .data[[t1_col]])
    }
    # T2_T0
    if (all(c(t2_col, t0_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T2_T0") := .data[[t2_col]] / .data[[t0_col]])
    }
    # T1_T0
    if (all(c(t1_col, t0_col) %in% names(df))) {
      df <- df %>%
        mutate(!!("ratio_T1_T0") := .data[[t1_col]] / .data[[t0_col]])
    }
  }
  
  return(df)
}


#' getsmooth
#'
#' @param dti variant data frame
#' @param test1 subset of dti to fit loess smoother
#' @param sspan loess span if less than 30 data point, 0.4 is used, otherwise use defined value.
#'
#' @return variant data frame with smoothed values
#' @export
#'
#' @examples getsmooth(dti, test1, sspan = 0.2)
getsmooth <- function(dti, test1, sspan = sspan) {
  if (nrow(test1) < 30) {
    sspan <- 0.4
  }
  loess15.1 <- loess(R ~ POS, data = test1, span = sspan)
  
  # Create a result vector
  smoothed_values <- rep(NA, length(dti$POS))
  
  # Find in-range indices
  fit_range <- range(loess15.1$x, na.rm = TRUE)
  in_range_idx <- which(dti$POS >= fit_range[1] & dti$POS <= fit_range[2])
  
  # Predict only for in-range values
  smoothed_result <- predict(
    loess15.1,
    newdata = data.frame(POS = dti$POS[in_range_idx]),
    se = TRUE,
    na.action = na.exclude
  )
  
  # Fill in predicted values for in-range positions
  smoothed_values[in_range_idx] <- smoothed_result$fit
  
  # Fill out-of-range values with the first predicted value
  first_fit <- smoothed_result$fit[1]
  smoothed_values[dti$POS < fit_range[1]] <- first_fit
  
  # Fill out-of-range high end with the last predicted value
  last_fit <- smoothed_result$fit[length(smoothed_result$fit)]
  smoothed_values[dti$POS > fit_range[2]] <- last_fit
  
  return(smoothed_values)
}


#' ladj
#'
#' @param dt variant data frame with ratios
#' @param loessSubset use all variant or subset with cf on time 2v1 ratio
#' @param cf cuff for ratio of time 2 vs 1 to be used to subset the total variant and for lowess fitting.
#' @param sspan loewss span
#' @param ratio which ratio to make the loess adjustment T2_T1 or T2_T0
#'
#' @return variant data frame with ratios with lowess adjusted ratios
#' @export
#'
#' @examples ladj(dt, loessSubset = "all", cf =0.8, sspan = 0.15, ratio = "T2_T1")
ladj <- function(
    dt,
    loessSubset = "all",
    cf = 0.8,
    sspan = 0.15,
    ratio = "T2_T1"
) {
  ## input should have ratio_T2_T1_* ratio_T2_T0_* ratio_T1_T0_*
  ## output will return the ratio_T2_T1_* ratio_T2_T0_* ratio_T1_T0_* adjusted
  ## loessSubset all or ModFindlay2
  ## if ModFindlay2, cf default 0.8 (it is the ratio ratio of T2 and T1)
  ## for each rep
  dt <- dt %>% arrange(POS)
  vtext <- ifelse(
    loessSubset == "all",
    "all variants",
    paste0("Modfindlay2 subset variants with cutoff of", cf)
  )
  cat(
    "...Running lowess adjustment with ",
    vtext,
    paste0(" span of ", sspan),
    " ... \n"
  )
  
  # each exon, each Rep
  exs <- unique(dt$Exon)
  
  out <- NULL
  for (exi in exs) {
    nreps <- max(unique(dt$Rep[dt$Exon == exi]))
    for (repi in 1:nreps) {
      Ri_T2_T1 <- ("ratio_T2_T1")
      Ri_T1_T0 <- ("ratio_T1_T0")
      Ri_T2_T0 <- ("ratio_T2_T0")
      
      T1_T0_loess <- ("T1_T0_loess")
      T2_T1_loess <- ("T2_T1_loess")
      
      T2_T0_ratio_adjusted <- ("ratio_T2_T0_adjusted")
      T2_T1_ratio_adjusted <- ("ratio_T2_T1_adjusted")
      
      dti <- dt %>% filter(Rep == repi & Exon == exi)
      if (loessSubset == "ModFindlay2") {
        ## time2-1 ratio
        test1 <- dti %>% filter(!!sym(Ri_T2_T1) > log2(cf))
      } else {
        ## all v no filter
        test1 <- dti
      }
      
      ############# loess adjusted for Time2-0
      
      if (ratio == "T2_T0") {
        test1 <- test1 %>% mutate(R = !!sym(Ri_T1_T0))
        smoothed15fit <- getsmooth(dti, test1, sspan = sspan)
        dti <- dti %>% mutate(!!sym(T1_T0_loess) := smoothed15fit)
        dti <- dti %>%
          mutate(
            !!sym(T2_T0_ratio_adjusted) := !!sym(Ri_T2_T0) - !!sym(T1_T0_loess)
          )
      }
      
      ################# loess for Time2-1
      
      if (ratio == "T2_T1") {
        test1 <- test1 %>% mutate(R = !!sym(Ri_T2_T1))
        smoothed15fit <- getsmooth(dti, test1, sspan = sspan)
        dti <- dti %>% mutate(!!sym(T2_T1_loess) := smoothed15fit)
        dti <- dti %>%
          mutate(
            !!sym(T2_T1_ratio_adjusted) := !!sym(Ri_T2_T1) - !!sym(T2_T1_loess)
          )
      }
      
      out <- bind_rows(out, dti)
    }
  }
  
  return(out)
}


#' scaling
#'
#' @param x data vector
#' @param pos_ctl.Avg.med pos control average median
#' @param neg_ctl.Avg.med neg control averagemedian
#' @param pos_ctl.med1  pos control median
#' @param neg_ctl.med1  neg control median
#'
#' @return scaled vector
#' @export
#'
#' @examples scaling(x, amp, amn, mp, mn)
scaling <- function(
    x,
    pos_ctl.Avg.med,
    neg_ctl.Avg.med,
    pos_ctl.med1,
    neg_ctl.med1
) {
  StopGain.Avg.med <- pos_ctl.Avg.med
  Syn.Avg.med <- neg_ctl.Avg.med
  StopGain.med1 <- pos_ctl.med1
  Syn.med1 <- neg_ctl.med1
  if (
    any(is_empty(c(Syn.Avg.med, StopGain.Avg.med, Syn.med1, StopGain.med1)))
  ) {
    retrun(x)
  } else {
    a <- (Syn.Avg.med - StopGain.Avg.med) / (Syn.med1 - StopGain.med1)
    b <- Syn.Avg.med - a * Syn.med1
    x2 <- a * x + b
    return(x2)
  }
}


#' norm_win_exon
#'
#' @param dt variant data frame with ratios
#' @param base log base
#' @param pos_group value of positive group in EventType Defined from config file
#' @param neg_group value of negative group in EventType Defined from config file
#' @param ratio value of ratio to test from config file
#' @param ladj if TRUE lowess positional adjustment will be performed, FALSE, no.
#'
#' @return variant data frame with ratio normalized within Exon
#' @export
#'
#' @examples norm_win_exon(dt, base = "n", pos_group, neg_group, ratio, ladj = TRUE)
norm_win_exon <- function(dt, base, pos_group, neg_group, ratio, ladj = TRUE) {
  ### base: could use SNV only
  ### posgroup: EventType use for positive controls, StopGain
  ### neggroup: EventType use for negative controls, Synonymous
  ### ratio: which ratio to normalize eg. ratio_T2_T1, ratio_T2_T0_adjusted, ratio_T2_T1_adjusted
  cat("...Running normalization within Exons...\n")
  outratio <- paste0(ratio, ifelse(ladj, "_adjusted", ""), "_norm_wie")
  ratio <- paste0("ratio_", ratio, ifelse(ladj, "_adjusted", ""))
  
  dt.out <- NULL
  
  for (ei in unique(dt$Exon)) {
    if (base == "SNV") {
      dti <- dt %>%
        filter(Exon == ei) %>%
        filter(nchar(REF) == 1 & nchar(ALT) == 1)
    } else {
      dti <- dt %>% filter(Exon == ei)
    }
    
    meds <- dti %>% ## remove variants in splicing region
      filter(is.na(spliceR) | spliceR == FALSE) %>%
      group_by(EventType, Rep) %>%
      dplyr::summarize(M = median(!!sym(ratio), na.rm = TRUE), .groups = "drop")
    
    meds2 <- meds %>%
      group_by(EventType) %>%
      dplyr::summarize(M2 = mean(M, na.rm = TRUE), .groups = "drop")
    
    dti$score1 <- NA
    for (ri in 1:max(dt$Rep)) {
      if (
        length(meds2$M2[which(meds2$EventType == neg_group)]) > 0 &
        length(meds2$M2[which(meds2$EventType == pos_group)]) > 0 &
        length(meds$M[which(meds$EventType == neg_group & meds$Rep == ri)]) >
        0 &
        length(meds$M[which(meds$EventType == pos_group & meds$Rep == ri)]) >
        0
      ) {
        dti[which(dti$Rep == ri), "score1"] <- scaling(
          x = unlist(dti[which(dti$Rep == ri), ratio]),
          pos_ctl.Avg.med = meds2$M2[which(meds2$EventType == pos_group)],
          neg_ctl.Avg.med = meds2$M2[which(meds2$EventType == neg_group)],
          pos_ctl.med1 = meds$M[which(
            meds$EventType == pos_group & meds$Rep == ri
          )],
          neg_ctl.med1 = meds$M[which(
            meds$EventType == neg_group & meds$Rep == ri
          )]
        )
      } else {
        dti[which(dti$Rep == ri), "score1"] <- dti[which(dti$Rep == ri), ratio]
      }
    }
    dt.out <- rbind(dt.out, dti)
  }
  names(dt.out)[names(dt.out) == "score1"] <- outratio
  return(dt.out)
}


#' norm_x_exon
#'
#' @param dt variant data frame with ratios
#' @param base log base
#' @param pos_group value of positive group in EventType Defined from config file
#' @param neg_group value of negative group in EventType Defined from config file
#' @param ratio value of ratio to test from config file
#' @param ladj if TRUE lowess positional adjustment will be performed, FALSE, no.
#'
#' @return variant data frame with ratio normalized acorss Exons
#' @export
#'
#' @examples norm_x_exon(dt, base, pos_group, neg_group, ratio, ladj = TRUE)
norm_x_exon <- function(dt, base, pos_group, neg_group, ratio, ladj = TRUE) {
  ### base: could use SNV only
  ### posgroup: EventType use for positive controls, Synonymous
  ### neggroup: EventType use for negative controls, StopGain
  ### ratio: which ratio to normalize eg. ratio_T2_T1, ratio_T2_T0_adjusted, ratio_T2_T1_adjusted
  cat("...Running normalization across Exons...\n")
  
  inratio <- paste0(ratio, ifelse(ladj, "_adjusted", ""), "_norm_wie")
  outratio <- paste0(ratio, ifelse(ladj, "_adjusted", ""), "_norm_xe")
  dt.out <- NULL
  
  dt.out <- dt
  if (base == "SNV") {
    dti.temp <- dt %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
  } else {
    dti.temp <- dt
  }
  
  meds <- dti.temp %>% ## remove synonymous variants in splicing region
    filter(!(EventType == neg_group & !is.na(spliceR) & spliceR > 0.2)) %>%
    group_by(Exon, EventType) %>%
    dplyr::summarize(M = median(!!sym(inratio), na.rm = TRUE), .groups = "drop")
  
  meds2 <- meds %>%
    group_by(EventType) %>%
    dplyr::summarize(M2 = mean(M, na.rm = TRUE), .groups = "drop")
  
  dt.out$score2 <- NA
  
  for (ei in unique(dt$Exon)) {
    if (
      length(meds2$M2[which(meds2$EventType == neg_group)]) > 0 &
      length(meds2$M2[which(meds2$EventType == pos_group)]) > 0 &
      length(meds$M[which(meds$EventType == neg_group & meds$Exon == ei)]) >
      0 &
      length(meds$M[which(meds$EventType == pos_group & meds$Exon == ei)]) > 0
    ) {
      dt.out[which(dt.out$Exon == ei), "score2"] <- scaling(
        x = unlist(dt.out[which(dt.out$Exon == ei), inratio]),
        pos_ctl.Avg.med = meds2$M2[which(meds2$EventType == pos_group)],
        neg_ctl.Avg.med = meds2$M2[which(meds2$EventType == neg_group)],
        pos_ctl.med1 = meds$M[which(
          meds$EventType == pos_group & meds$Exon == ei
        )],
        neg_ctl.med1 = meds$M[which(
          meds$EventType == neg_group & meds$Exon == ei
        )]
      )
    } else {
      dt.out[which(dt.out$Exon == ei), "score2"] <- dt.out[
        which(dt.out$Exon == ei),
        inratio
      ]
    }
  }
  names(dt.out)[names(dt.out) == "score2"] <- outratio
  return(dt.out)
}


#' calculate_auc_balance
#'
#' @param predictions vector with prediction numeric values
#' @param labels value of case or control
#' @param plot if TRUE roc plot will be made
#' @param plot_title Title for the ROC plot
#'
#' @return a list of object : auc, balance_threshold, sensitivity, specificity, 2x2 cross tab and a plot
#' @export
#'
#' @examples calculate_auc_balance(predictions, labels, plot = TRUE, plot_title = "ROC Curve")
calculate_auc_balance <- function(
    predictions,
    labels,
    plot = TRUE,
    plot_title = "ROC Curve"
) {
  # Input validation
  if (!all(labels %in% c("case", "control"))) {
    stop("Labels must be 'case' or 'control'")
  }
  
  # Convert to binary format
  binary_labels <- ifelse(labels == "case", 1, 0)
  
  # Calculate ROC curve
  roc_obj <- roc(response = binary_labels, predictor = predictions)
  
  # Get AUC value
  auc_value <- auc(roc_obj)
  
  # Extract threshold metrics
  coords_df <- coords(
    roc_obj,
    x = "all",
    ret = c("threshold", "sensitivity", "specificity")
  )
  
  # Find optimal balance point
  coords_df$diff <- abs(coords_df$sensitivity - coords_df$specificity)
  balance_row <- coords_df[which.min(coords_df$diff), ]
  
  # Generate ROC plot if requested
  if (plot) {
    plot_output <- ggroc(roc_obj, color = "#2c7bb6", size = 1) +
      geom_segment(
        aes(x = 1, xend = 0, y = 0, yend = 1),
        color = "darkgray",
        linetype = "dashed"
      ) +
      geom_point(
        aes(x = balance_row$specificity, y = balance_row$sensitivity),
        color = "#d7191c",
        size = 3
      ) +
      annotate(
        "text",
        x = balance_row$specificity - 0.1,
        y = balance_row$sensitivity - 0.2,
        label = paste0(
          "Balance Point\n",
          "Sens: ",
          round(balance_row$sensitivity, 2),
          "\n",
          "Spec: ",
          round(balance_row$specificity, 2)
        ),
        color = "#d7191c"
      ) +
      labs(
        title = plot_title,
        subtitle = paste("AUC =", round(auc_value, 3)),
        x = "1 - Specificity",
        y = "Sensitivity"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)
      )
  } else {
    plot_output <- NULL
  }
  
  pred <- ifelse(predictions >= balance_row$threshold, "case", "control")
  tb <- table(labels, pred)
  
  # Return results
  return(list(
    auc = round(auc_value, 4),
    balance_threshold = balance_row$threshold,
    cross_tb = tb,
    sensitivity = balance_row$sensitivity,
    specificity = balance_row$specificity,
    plot = plot_output
  ))
}

#' calFS
#'
#' @param dt variant data frame with cross excon normalized scores
#' @param scorename which normalized scores to use
#'
#' @return variant data frame with functional score calculated
#' @export
#'
#' @examples calFS(dt, "T2_T1_norm_xe")
calFS <- function(dt, scorename) {
  ## get functional score, ie average, from score name column
  ## input: long format by Rep
  ## output: wide format per exon, per variant.
  cat("...Running FS calculation...\n")
  if (!scorename %in% colnames(dt)) {
    stop("ratio is not present in the data provided")
  }
  dt <- dt %>% mutate(uPOS = paste(POS, REF, ALT, sep = "_"))
  id_cols_vec <- setdiff(names(dt), c("Rep", scorename))
  value_cols <- setdiff(
    names(dt),
    c(
      "uPOS",
      "Rep",
      "CHROM",
      "POS",
      "REF",
      "ALT",
      "EventType",
      "Exon",
      "use4norm",
      "spliceR"
    )
  )
  
  wide <- dt %>%
    arrange(uPOS) %>%
    tidyr::pivot_wider(
      names_from = c(Rep),
      values_from = all_of(value_cols),
      names_glue = "{.value}_R{Rep}"
    )
  wide <- wide %>%
    rowwise() %>%
    mutate(functional_score = mean(c_across(matches(scorename)), na.rm = TRUE))
  return(wide)
}


#' posteriorProbabilityFromGaussianModels
#'
#' @param effect effect to be used
#' @param standardError Standard error to be used
#' @param meanAndSD object with mu and sd
#' @param prior prior to be used
#'
#' @return a list of posteriorP and  likelihoodRatio
#' @export
#'
#' @examples posteriorProbabilityFromGaussianModels(effect, standardError,meanAndSD, prior)
posteriorProbabilityFromGaussianModels = function(
    effect,
    standardError,
    meanAndSD,
    prior
) {
  # calculate the posterior probabilities
  
  # effect: the estimated effect
  # standardError: the estimated standard error of the estimated effect
  # meanAndSD: a list with mu and sd for each component
  # prior: the prior vector for each variant
  
  mu = meanAndSD$mu
  sd = meanAndSD$sd
  # calculate the likelihood
  logLikelihood1 = dnorm(
    effect,
    mu[1],
    sqrt(sd[1]^2 + standardError^2),
    log = T
  )
  logLikelihood2 = dnorm(
    effect,
    mu[2],
    sqrt(sd[2]^2 + standardError^2),
    log = T
  )
  
  logLikelihoodRatio = logLikelihood1 - logLikelihood2
  likelihoodRatio = exp(logLikelihoodRatio)
  
  posteriorP = (likelihoodRatio * prior) / ((likelihoodRatio - 1) * prior + 1)
  
  return(list(posteriorP = posteriorP, likelihoodRatio = likelihoodRatio))
}


#' getfs
#'
#' @param p1.th theoretical posterior probability
#' @param fs.th theoretical functional score
#' @param pc probability cutoff
#'
#' @returns corresponding functional score given the cutoff
#' @export
#'
#' @examples getfs(p1.th, fs.th, pc = 0.98)
getfs <- function(p1.th, fs.th, pc) {
  ## get the functional score corresponding to posterior probability
  minp <- min(unlist(p1.th['posteriorP']), na.rm = T) ## minimal of probability
  maxf <- min(fs.th[which(unlist(p1.th['posteriorP']) == minp)], na.rm = T) ## corresponding minimal functional score
  fs.th.short <- fs.th[which(fs.th <= maxf)]
  p1.th.short <- unlist(p1.th['posteriorP'])[which(fs.th <= maxf)]
  if (minp > pc) {
    return(NA) ##return(maxf) changed this after discussion with Wenan 4/8/2025
  } else {
    cutoff <- fs.th[tail(which(p1.th.short > pc), 1)]
    if (length(cutoff) == 0) {
      return(NA)
    } else {
      return(cutoff)
    }
  }
}


#' mm
#'
#' @param data
#' @param typecolumn Column contains values to indicate groups
#' @param neg_group what value to use to identify negative group in typecolumn
#' @param pos_group what value to use to identify positive group in typecolumn
#' @param lambda prior, ratio of pos_group to overall, if "data" use data to calculate ratio
#' @param MM_OPTION value of "normal"or "robust", pre-estimate of two groups mean and sd using normal or robust methods.
#' @param base what to use for the training data, "allv" or "SNV" all variants or snv only used for traning.
#' @param se.Qcutoff subset training variants to those below qth quantile of se
#'
#' @returns a list of objects containing:
#' \describe{
#'   \item{error.occured}{TRUE of FALSE, if any error occurs during the modeling.}
#'   \item{f1}{NA, or A numeric number, functional score at 1% probablity of being pathogentic}
#'   \item{f99}{NA, or A numeric number, functional score at 99% probablity of being pathogentic}
#'   \item{f2}{NA, or A numeric number, functional score at 2% probablity of being pathogentic}
#'   \item{f98}{NA, or A numeric number, functional score at 98% probablity of being pathogentic}
#'   \item{f5}{NA, or A numeric number, functional score at 5% probablity of being pathogentic}
#'   \item{f95}{NA, or A numeric number, functional score at 95% probablity of being pathogentic}
#'   \item{f10}{NA, or A numeric number, functional score at 10% probablity of being pathogentic}
#'   \item{f90}{NA, or A numeric number, functional score at 90% probablity of being pathogentic}
#'   \item{f20}{NA, or A numeric number, functional score at 20% probablity of being pathogentic}
#'   \item{f80}{NA, or A numeric number, functional score at 80% probablity of being pathogentic}
#'   \item{data}{data used for analysis}
#'   \item{plot1}{plot of functional score (x) vs posterior probablity of non-functional/pathogentic}
#'   \item{plot2}{plot of density of functional scores of positive and negative groups and cutoffs of fitted values}
#' }
#'
#' @export
#'
#' @examples mm(data=data, typecolumn = "EventType", neg_group = "StopGain", pos_group = "Synonymous", lambda = "data", MM = "normal", base = "allv")
mm <- function(
    data,
    typecolumn,
    neg_group,
    pos_group,
    lambda = "data",
    MM = "normal",
    base = "allv",
    se.Qcutoff = NA
) {
  ## plot mixture model
  ## Get thresholds and return posterior prob,
  ## require data have functional_score
  ## must have functional_score_se
  cat("...Running Mixture Model on Funtional Scores...\n")
  if (!"functional_score_se" %in% colnames(data)) {
    data$functional_score_se <- 0
  }
  
  data <- data %>%
    mutate(uPOS = paste(POS, REF, ALT, sep = "_")) %>%
    mutate(posteriorP = NA, likelihoodRatio = NA)
  f1A <- f2A <- f3A <- f4A <- f5A <- f6A <- f7A <- f8A <- f9A <- f10A <- NA
  
  data0 <- data %>% filter(!is.na(functional_score))
  
  temp <- data0 %>%
    filter(!!sym(typecolumn) %in% c(neg_group, pos_group)) %>%
    filter(is.na(spliceR) | spliceR == FALSE)
  
  if (base == "SNV") {
    temp <- temp %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
  }
  
  ## check if functional_score_se is in the file.
  
  if (!is.na(se.Qcutoff)) {
    temp <- temp %>%
      filter(
        functional_score_se <=
          quantile(
            temp$functional_score_se,
            probs = c(se.Qcutoff),
            na.rm = TRUE
          )
      )
  }
  
  myind <- match(temp$uPOS, data0$uPOS) ## index of variants used for mm
  
  m1 <- mean(
    temp$functional_score[which(pull(temp, typecolumn) == neg_group)],
    na.rm = TRUE
  )
  tmp <- temp$functional_score[which(pull(temp, typecolumn) == neg_group)]
  tmp2 <- tmp
  
  v1 <- sd(tmp2, na.rm = TRUE)
  
  m2 <- mean(temp$functional_score[which(pull(temp, typecolumn) == pos_group)])
  v2 <- sd(temp$functional_score[which(pull(temp, typecolumn) == pos_group)])
  l1 <- sum(pull(temp, typecolumn) == pos_group) /
    sum(pull(temp, typecolumn) == pos_group)
  
  if (MM_OPTION == "Robust") {
    # robust calculation
    m2 = mean(
      temp$functional_score[which(pull(temp, typecolumn) == pos_group)],
      trim = 0.05
    )
    m1 = mean(
      temp$functional_score[which(pull(temp, typecolumn) == neg_group)],
      trim = 0.05
    )
    v2 = IQR(temp$functional_score[which(
      pull(temp, typecolumn) == pos_group
    )]) /
      1.349
    v1 = IQR(temp$functional_score[which(
      pull(temp, typecolumn) == neg_group
    )]) /
      1.349
  }
  error.occured <- FALSE
  tryCatch(
    {
      x <- normalmixEM(
        temp$functional_score,
        lambda = l1,
        k = 2,
        mean.constr = c(m2, m1),
        sd.constr = c(v2, v1),
        mu = c(m2, m1),
        sigma = c(v2, v1)
      )
      m1 = list(mu = x$mu, sd = x$sigma)
      
      prior = ifelse(
        lambda == "0.1",
        0.1,
        ifelse(lambda == "0.2", 0.2, x$lambda[1])
      )
      print(paste0("MM prior used is: ", prior))
      
      p1 = posteriorProbabilityFromGaussianModels(
        data0$functional_score,
        rep(0, nrow(data0)),
        m1,
        prior = prior
      )
      fs.th <- seq(
        from = min(data0$functional_score),
        to = max(data0$functional_score),
        length.out = 10000
      )
      p1.th <- posteriorProbabilityFromGaussianModels(
        fs.th,
        rep(0, 10000),
        m1,
        prior = prior
      )
      data0$posteriorP <- unlist(p1["posteriorP"])
      data0$likelihoodRatio <- unlist(p1["likelihoodRatio"])
      data0$mmused <- 0
      data0$mmused[myind] <- 1
      suppressMessages({
        temp1 <- data0[, c(
          "uPOS",
          "posteriorP",
          "likelihoodRatio",
          "mmused",
          typecolumn
        )] %>%
          unique()
        temp2 <- left_join(
          data[,
               !names(data) %in% c("posteriorP", "likelihoodRatio", typecolumn)
          ],
          temp1,
          by = "uPOS"
        )
      })
      
      data <- temp2 %>% mutate(pos_mean = x$mu[1], pos_sd = x$sigma[1], neg_mean = x$mu[2], neg_sd = x$sigma[2])
      
      print(paste0(
        "MM Difference from data posterior < 1e-08: ",
        max(
          abs(unlist(p1['posteriorP'])[myind] - x$posterior[, 1]),
          na.rm = T
        ) <
          1e-08
      ))
      
      f1A <- getfs(p1.th, fs.th, 0.01)
      f2A <- getfs(p1.th, fs.th, 0.99)
      
      f3A <- getfs(p1.th, fs.th, 0.02)
      f4A <- getfs(p1.th, fs.th, 0.98)
      
      f5A <- getfs(p1.th, fs.th, 0.05)
      f6A <- getfs(p1.th, fs.th, 0.95)
      
      f7A <- getfs(p1.th, fs.th, 0.10)
      f8A <- getfs(p1.th, fs.th, 0.90)
      
      f9A <- getfs(p1.th, fs.th, 0.20)
      f10A <- getfs(p1.th, fs.th, 0.80)
      
      plot_df <- data.frame(
        x = seq(
          from = min(data0$functional_score),
          to = max(data0$functional_score),
          length.out = 10000
        ),
        y = unlist(p1.th['posteriorP'])
      )
      
      # Create points data frame
      points_df <- data.frame(
        x = temp$functional_score,
        y = unlist(p1['posteriorP'])[myind],
        type = ifelse(
          pull(temp, typecolumn) == neg_group,
          "Negative",
          "Positive"
        )
      )
      # Build plot
      gg1 <- ggplot() +
        # Main probability curve
        geom_line(data = plot_df, aes(x = x, y = y), linewidth = 0.8) +
        
        # Points with conditional coloring
        geom_point(
          data = points_df,
          aes(x = x, y = y, color = type),
          size = 1.5,
          shape = 1
        ) +
        scale_color_manual(
          values = c("Positive" = "red", "Negative" = "blue"),
          labels = c(neg_group, pos_group)
        )
      
      # Vertical lines (grouped by linetype)
      if (!is.na(f1A) & !is.na(f2A)) {
        gg1 <- gg1 +
          geom_vline(
            xintercept = c(f1A, f2A),
            linetype = "solid",
            linewidth = 0.8
          )
      }
      if (!is.na(f3A) & !is.na(f4A)) {
        gg1 <- gg1 +
          geom_vline(
            xintercept = c(f3A, f4A),
            linetype = "dashed",
            linewidth = 0.8
          )
      }
      if (!is.na(f5A) & !is.na(f6A)) {
        gg1 <- gg1 +
          geom_vline(
            xintercept = c(f5A, f6A),
            linetype = "dotted",
            linewidth = 0.8
          )
      }
      if (!is.na(f7A) & !is.na(f8A)) {
        gg1 <- gg1 +
          geom_vline(
            xintercept = c(f7A, f8A),
            linetype = "dotdash",
            linewidth = 0.8
          )
      }
      if (!is.na(f9A) & !is.na(f10A)) {
        gg1 <- gg1 +
          geom_vline(
            xintercept = c(f9A, f10A),
            linetype = "longdash",
            linewidth = 0.8
          )
      }
      # Horizontal lines (grouped by linetype)
      gg1 <- gg1 +
        geom_hline(
          yintercept = c(0.01, 0.99),
          linetype = "solid",
          linewidth = 0.8
        ) +
        geom_hline(
          yintercept = c(0.02, 0.98),
          linetype = "dashed",
          linewidth = 0.8
        ) +
        geom_hline(
          yintercept = c(0.05, 0.95),
          linetype = "dotted",
          linewidth = 0.8
        ) +
        geom_hline(
          yintercept = c(0.10, 0.90),
          linetype = "dotdash",
          linewidth = 0.8
        ) +
        geom_hline(
          yintercept = c(0.20, 0.80),
          linetype = "longdash",
          linewidth = 0.8
        ) +
        
        # Labels and theme
        labs(
          x = "Functional score",
          y = "Posterior probability of non-functional",
          color = "Group"
        ) +
        theme_minimal() +
        theme(
          legend.position = "top",
          panel.grid = element_blank() # Removes both major and minor grid lines
        )
      
      mix_result <- x
      plot_data <- data.frame(x = mix_result$x)
      components <- data.frame(
        mu = mix_result$mu,
        sigma = mix_result$sigma,
        lambda = mix_result$lambda
      )
      
      # Create ggplot
      gg2 <- ggplot(plot_data, aes(x = x)) +
        geom_histogram(
          aes(y = ..density..),
          bins = 30,
          fill = "gray",
          alpha = 0.7
        ) +
        stat_function(
          fun = function(x) {
            components$lambda[1] *
              dnorm(x, components$mu[1], components$sigma[1])
          },
          color = "red"
        ) +
        stat_function(
          fun = function(x) {
            components$lambda[2] *
              dnorm(x, components$mu[2], components$sigma[2])
          },
          color = "blue"
        ) +
        geom_density(linetype = "dashed", linewidth = 0.8)
      
      # Vertical lines with different linetypes
      if (!is.na(f1A) & !is.na(f2A)) {
        gg2 <- gg2 +
          geom_vline(
            xintercept = c(f1A, f2A),
            linetype = "solid",
            linewidth = 0.8
          )
      }
      if (!is.na(f3A) & !is.na(f4A)) {
        gg2 <- gg2 +
          geom_vline(
            xintercept = c(f3A, f4A),
            linetype = "dashed",
            linewidth = 0.8
          )
      }
      if (!is.na(f5A) & !is.na(f6A)) {
        gg2 <- gg2 +
          geom_vline(
            xintercept = c(f5A, f6A),
            linetype = "dotted",
            linewidth = 0.8
          )
      }
      if (!is.na(f7A) & !is.na(f8A)) {
        gg2 <- gg2 +
          geom_vline(
            xintercept = c(f7A, f8A),
            linetype = "dotdash",
            linewidth = 0.8
          )
      }
      if (!is.na(f9A) & !is.na(f10A)) {
        gg2 <- gg2 +
          geom_vline(
            xintercept = c(f9A, f10A),
            linetype = "longdash",
            linewidth = 0.8
          )
      }
      gg2 <- gg2 +
        labs(
          title = "Two components prob(nonfunctional) with cutoffs",
          x = "Functional Scores"
        ) +
        theme_minimal()
    },
    error = function(e) {
      error.occured <<- TRUE
      print(e)
    }
  )
  if (error.occured) {
    gg1 <- NULL
    gg2 <- NULL
  }
  return(list(
    error.occured = error.occured,
    f1 = as.numeric(f1A),
    f99 = as.numeric(f2A),
    f2 = as.numeric(f3A),
    f98 = as.numeric(f4A),
    f5 = as.numeric(f5A),
    f95 = as.numeric(f6A),
    f10 = as.numeric(f7A),
    f90 = as.numeric(f8A),
    f20 = as.numeric(f9A),
    f80 = as.numeric(f10A),
    data = data,
    plot1 = gg1,
    plot2 = gg2
  ))
}


#' render_report
#'
#' @param Input_file result file from mixture model
#' @param cfgfile configuration file
#'
#' @returns NULL
#' @export
#'
#' @examples render_report(Input_file, cfgfile)
render_report <- function(Input_file, cfgfile) {
  #define main directory
  source(cfgfile)
  main_dir <- OUTDIR
  
  #create directory if it doesn't exist
  if (!dir.exists(main_dir)) {
    dir.create(main_dir)
  }
  
  #define sub directory
  sub_dir <- "reports"
  
  #define full directory
  full_dir <- file.path(main_dir, sub_dir)
  
  #create directory if it doesn't exist
  if (!dir.exists(full_dir)) {
    dir.create(full_dir)
  }
  
  rmarkdown::render(
    "./GM_report.Rmd",
    params = list(
      inputfile = Input_file,
      cfgfile = cfgfile
    ),
    output_file = paste0(full_dir, "/", "MM.html"),
    quiet = TRUE
  )
}


#' MASSAGE_GM
#'
#' GM-only version of the MASSAGE() wrapper. Runs the mixture-model
#' pipeline: read data -> getspliceMax -> normfreq -> getratio ->
#' (optional) ladj -> norm_win_exon -> norm_x_exon -> calFS -> mm ->
#' (optional) render_report. The VarCall (JAGS) and Deseq branches
#' from the original MASSAGE() have been removed.
#'
#' @param cfg path to a config file that is source()'d and must define:
#'   DATAFILE, OUTDIR, LOWESS_ADJ, LOWESS_SUBSET, ModFindlay_CF, LOWESS_SPAN,
#'   RATIO_COMP, NORM_BASE, POS_GROUP, NEG_GROUP, MM_LAMBDA, MM_REPORT
#'
#' @return NULL (writes MM.RData and MM_result.csv to OUTDIR; a report
#'   in OUTDIR/reports/MM.html if MM_REPORT is TRUE)
#' @export
#'
#' @examples MASSAGE_GM("config.R")
MASSAGE_GM <- function(cfg) {
  ## workflow to analyse the data based on config file, MM analysis only.
  source(cfg)
  cat("Analysis is started... \n\n")
  cat("...Reading inputs... \n")
  
  ## readin data
  testdata <- read_delim(
    DATAFILE,
    delim = "\t",
    escape_double = FALSE,
    trim_ws = TRUE,
    show_col_types = FALSE
  )
  
  testdata0 <- testdata %>%
    dplyr::select(
      CHROM,
      POS,
      REF,
      ALT,
      EventCount,
      EventType,
      SpliceMax,
      Rep,
      Time,
      Exon
    )
  testdata0 <- getspliceMax(testdata0)
  
  ## terminate if wrong setting
  if (!LOWESS_ADJ & RATIO_COMP == "T2_T0") {
    return(
      "...Exiting for wrong setting.... T2_T0 need LOWESS_ADJ to be TRUE..."
    )
  }
  
  ## normalize counts
  norm_counts <- normfreq(testdata0)
  
  ## get ratio
  raw_ratios <- getratio(norm_counts)
  
  ## check config so that we know what to do with the data
  if (LOWESS_ADJ) {
    ## T2_T1 with lowess adjustment
    adj_ratio <- ladj(
      raw_ratios,
      loessSubset = LOWESS_SUBSET,
      cf = ModFindlay_CF,
      sspan = LOWESS_SPAN,
      ratio = RATIO_COMP
    )
  } else {
    adj_ratio <- raw_ratios
  }
  
  ## norm wi exon
  tempdt <- norm_win_exon(
    adj_ratio,
    base = NORM_BASE,
    pos_group = POS_GROUP,
    neg_group = NEG_GROUP,
    ratio = RATIO_COMP,
    ladj = LOWESS_ADJ
  )
  
  ## norm xexon
  norm_ratios <- norm_x_exon(
    tempdt,
    base = NORM_BASE,
    pos_group = POS_GROUP,
    neg_group = NEG_GROUP,
    ratio = RATIO_COMP,
    ladj = LOWESS_ADJ
  )
  
  ## calculate FS
  FS_data <- calFS(
    norm_ratios,
    paste0(RATIO_COMP, ifelse(LOWESS_ADJ, "_adjusted", ""), "_norm_xe")
  )
  
  ## Analysis: MM
  mymm <- mm(
    data = FS_data,
    typecolumn = "EventType",
    neg_group = NEG_GROUP,
    pos_group = POS_GROUP,
    lambda = MM_LAMBDA,
    MM = "normal",
    base = "allv"
  )
  
  ## output
  save(norm_ratios, FS_data, mymm, file = paste0(OUTDIR, "/MM.RData"))
  
  ## generate report
  if (MM_REPORT) {
    cat("...Running report generation... \n")
    render_report(Input_file = paste0(OUTDIR, "/MM.RData"), cfgfile = cfg)
  }
  
  write.csv(mymm$data, file = paste0(OUTDIR, "/MM_result.csv"))
  cat(
    "\nCheck your result in this file MM_result.csv, need more details regarding model inputs check MM.RData, more visualization check report.\n"
  )
  
  cat("\nAnalyisis is DONE... Thank you.\n")
}
