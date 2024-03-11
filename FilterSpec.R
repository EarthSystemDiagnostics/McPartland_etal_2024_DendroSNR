#' Filter a Power Spectrum Object
#'
#' @param spec A spec object
#' @param keep_low_f Keep filtered (smoothed) low frequencies or replace with unfiltered
#' @inheritParams stats::spec.pgram
#' @inheritParams PaleoSpec::ApplyFilter
#' @return A spec object (list)
#' @export
#'
#' @examples
#' ## Comparison of the four methods - for power spectra, methods 0, 2 or 3 make the most sense
#' library(PaleoSpec)
#' 
#' a <- 100
#' b <- 1
#' N <- 1e03
#' set.seed(20230625)
#' ts1 <- SimPLS(N, beta = b, alpha = a)
#' sp1 <- SpecMTM(ts(ts1), bin.width = 1)
#' LPlot(sp1)
#' abline(log10(a), -b, col = "green")

#' fl <- seq(3, 9, by = 2)
#' sp1_f3_0 <- FilterSpec(sp1, spans = fl, method = 0)
#' sp1_f3_1 <- FilterSpec(sp1, spans = fl, method = 1)
#' sp1_f3_2 <- FilterSpec(sp1, spans = fl, method = 2)
#' sp1_f3_3 <- FilterSpec(sp1, spans = fl, method = 3)
#' sp1_f3_4 <- FilterSpec(sp1, spans = fl, method = 4)
#'  
#' LPlot(sp1)
#' LLines(sp1_f3_0, col = "blue")
#' LLines(sp1_f3_1, col = "green", lty = 2)
#' LLines(sp1_f3_2, col = "red", lty = 3)
#' LLines(sp1_f3_3, col = "orange", lty = 4)
#' LLines(sp1_f3_4, col = "gold", lty = 5)
#' 
#' ## Comparison of keeping the filtered values in the reflected end portions or not 
#' sp1_f3_0T <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = TRUE)
#' sp1_f3_0F <- FilterSpec(sp1, spans = fl, method = 0, keep_low_f = FALSE)

#' LPlot(sp1_f3_0F)
#' LLines(sp1_f3_0T, col = "red")

#' sp1_f3_2T <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = TRUE)
#' sp1_f3_2F <- FilterSpec(sp1, spans = fl, method = 2, keep_low_f = FALSE)

#' LPlot(sp1_f3_2F)
#' LLines(sp1_f3_2T, col = "red")
FilterSpec <- function(spec, spans, method = 3, keep_low_f = TRUE) {
  if (length(spec$dof) == 1) {
    spec$dof <- rep(spec$dof, length(spec$freq))
  }
  
  dof0 <- spec$dof
  
  kernel <- stats::kernel("modified.daniell", spans %/% 2)
  filter <- kernel[-kernel$m:kernel$m]
  
  spec_filt <- ApplyFilter(spec$spec, filter = filter, method = method)
  
  if (keep_low_f == FALSE) {
    # replace filtered spec with original in area where freqs have been reflected
    i <- 1:ceiling(length(filter) / 2)
    spec_filt[i] <- spec$spec[i]
    
    iend <- length(spec$freq) - (i-1)
    
    spec_filt[iend] <- spec$spec[iend]
    
  }
  
  spec$spec <- as.numeric(spec_filt)
  
  df.kern <- stats::df.kernel(kernel)
  
  spec$dof <- df.kern * spec$dof / 2
  
 if (keep_low_f == FALSE) {
   
   i <- 1:ceiling(length(filter) / 2)
   spec$dof[i] <- dof0[i]
   
   iend <- length(spec$freq) - (i-1)
   spec$dof[iend] <- dof0[iend]
   
 }
  
  # Adjust DOF in reflected filter region
  if (keep_low_f == TRUE){
    
    fl <- length(filter)
    i <- 1:ceiling(fl / 2)
    iend <- length(spec$freq) - (i-1)
    
    
    if (method %in% c(2,3)){
     scl <- 2 * (fl - (i - 1)) / fl
     spec$dof[i] <- spec$dof[i] / scl
     spec$dof[iend] <- spec$dof[iend] / scl
     
    }
    
    if (method == 0){
      # remove NA portion 
      spec$freq <- spec$freq[is.na(spec$spec) == FALSE]
      spec$dof <- spec$dof[is.na(spec$spec) == FALSE]
      spec$shape <- spec$shape[is.na(spec$spec) == FALSE]
      spec$spec <- spec$spec[is.na(spec$spec) == FALSE]
    }
    
    }
  
   spec$shape <- spec$dof / 2
 
   spec <- PaleoSpec::AddConfInterval(spec)
   
  
  return(spec)
}



#' Smooth a Spectrum with Evenly Spaced Bins in Logspace
#'
#' @param spec A spec object
#' @inheritParams PaleoSpec::LogSmooth
#' @inheritParams PaleoSpec::ApplyFilter#'
#' @return A spec object (list)
#' @export
#' @examples
#' library(PaleoSpec)
#' 
#' # simulate a timeseries with powerlaw power spectrum
#' a <- 100
#' b <- 1
#' N <- 1e03
#' 
#' set.seed(20230625)
#' ts1 <- SimPLS(N, beta = b, alpha = a)
#' sp1 <- SpecMTM(ts(ts1), bin.width = 1)
#' LPlot(sp1)
#' abline(log10(a), -b, col = "green")
#' # 
#' sp1_f3_0 <- FilterSpecLog(sp1, method = 0)
#' sp1_f3_2 <- FilterSpecLog(sp1, method = 2)
#'  
#' LPlot(sp1)
#' LLines(sp1_f3_0, col = "blue")
#' LLines(sp1_f3_2, col = "green", lty = 3)
#' 
#' sp1_df0.05 <- FilterSpecLog(sp1)
#' sp1_df0.1 <- FilterSpecLog(sp1, df.log = 0.1)
#'  
#' LPlot(sp1)
#' LLines(sp1_df0.05, col = "blue")
#' LLines(sp1_df0.1, col = "red")
#' 
#' ## A combination of FilterSpec and FilterSpecLog
#' 
#' sp1_FSL <- FilterSpecLog(sp1)
#' sp1_FSL_FS <- FilterSpec(FilterSpecLog(sp1), spans = c(3, 5))
#' sp1_FS_FSL <- FilterSpecLog(FilterSpec(sp1, spans = c(3, 5)))
#' LPlot(sp1)
#' LLines(sp1_FSL, col = "blue")
#' LLines(sp1_FSL_FS, col = "red")
#' LLines(sp1_FS_FSL, col = "green")
FilterSpecLog <- function(spec,
                          df.log = 0.05,
                          spans = NULL,
                          method = 3, f.res = 10){
  
  GetFW <- function(spec, df.log) {
    ((exp(df.log) - 1) * max(spec$freq)) / min(spec$freq)
  }
  
  if (length(spec$dof) == 1){
    spec$dof <- rep(spec$dof, length(spec$freq))
  }
  
  if (is.null(spans)){
    spans <- GetFW(spec, df.log = df.log)
  }
  
  # interpolate spectrum onto equal in log space freq axis
  delta_f <- min(spec$freq) 
  logfreq <- log(spec$freq)
  
  freq_logspace <- (seq(min(logfreq), max(logfreq)+delta_f, length.out = f.res*length(spec$freq)))
  spec_loginterp <- stats::approx(logfreq, spec$spec, xout = freq_logspace, rule = 2)$y
  
  spans_adj <- spans * f.res
  
  # DOF of filter
  kernel <- stats::kernel("daniell", spans_adj %/% 2)
  filter <- kernel[-kernel$m:kernel$m]
  df.kern <- stats::df.kernel(kernel)
  
  # DOF of a boxcar filter the same width
  kernal.flat <- stats::kernel("daniell", length(filter) %/% 2)
  df.kern.flat <- stats::df.kernel(kernal.flat)
  
  # modify for non boxcar filters
  df.mod <- df.kern / df.kern.flat
  
  # smooth/filter in log space
  spec_filt <- ApplyFilter(spec_loginterp, filter = filter, method = method)
  
  # re-interpolate back to original freq axis
  spec3 <- stats::approx(freq_logspace, spec_filt, xout = logfreq)$y
  
  # overwrite spec with filtered spec
  spec$spec <- spec3
 
  # keep old DOF
  dof0 <- spec$dof
  
  # Gets the difference in delta_f for the log and standard freq axis
  NpF <- function(freq, fw, df){
    
    posdiff <- (exp(log(freq) + df) - freq)
    negdiff <- (freq - exp(log(freq) - df))
    
    fdiff <- rowMeans(cbind(negdiff, posdiff))
    
    2 * fw * (fdiff/df) * 1/(2*max(freq))
  }
  df.logkern <- NpF(spec$freq, length(filter), df = diff(freq_logspace[1:2]))
  
  spec$dof <- spec$dof + df.mod * df.logkern * spec$dof/2
  spec$shape <- spec$dof/2
  spec$spans <- paste(spans, collapse = ",")
  
  if (method == 0){
    # remove NA portion 
    spec$freq <- spec$freq[is.na(spec$spec) == FALSE]
    spec$dof <- spec$dof[is.na(spec$spec) == FALSE]
    spec$shape <- spec$shape[is.na(spec$spec) == FALSE]
    spec$spec <- spec$spec[is.na(spec$spec) == FALSE]
  }
  
  spec <- PaleoSpec::AddConfInterval(spec)
  
  
  
  return(spec)
}
