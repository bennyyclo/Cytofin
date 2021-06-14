#' Normalize CyTOF data panel
#'
#' This function normalize CyTOF data panel using external anchors.
#'
#' @param control_metadata_path metadata table of anchor CyTOF files (.fcs) 
#' (must be in the current directory).
#' @param anchor_statistics .rds object containing anchor reference statistics 
#' (must be in the current directory).
#' @param sample_metadata_path metadata table of homogenized CyTOF files 
#' (.fcs) (must be in the current directory).
#' @param panel_path standardized antigen panel table file (.xlsx/.csv)
#' (must be in the current directory).
#' @param input_path folder directory containing input homogenized CyTOF 
#' data file.
#' @param val_path folder directory containing validation homogenized CyTOF 
#' data file (optional).
#' @param output_path folder directory containing output normalized CyTOF data file.
#' @param mode transformation function used for normalization 
#' (meanshift, meanshift_bulk, variance, z_score, or beadlike).
#'
#' @return normalized files with specified panel
#' 
#' 
#' @export
cytofin_normalize <-
  function(
    control_metadata_path, 
    anchor_statistics, 
    sample_metadata_path, 
    panel_path, 
    input_path, 
    val_path = "none",
    output_path, 
    mode = c("meanshift", "meanshift_bulk", "variance", "z_score", "beadlike")
  ) {
    
    # create output directory
    dir.create(output_path)
    
    #read metadata tables
    md_control <- cytofin_read_metadata(control_metadata_path)
    md <- cytofin_read_metadata(sample_metadata_path)
    
    # if anchor_statistics is a file path
    if (is.character(anchor_statistics)) { 
      anchor_statistics_list <- readr::read_rds(anchor_statistics)
    } else if (is.list(anchor_statistics)) { 
      anchor_statistics_list <- anchor_statistics
    } else { 
      stop("anchor_statistics must be either a character vector (file path) or a list")
    }
    
    # extract needed values from anchor_statistics_list
    var_uni <- anchor_statistics_list$var_uni
    mean_uni <- anchor_statistics_list$mean_uni
    var_uni_mean <- anchor_statistics_list$var_uni_mean
    mean_uni_mean <- anchor_statistics_list$mean_uni_mean
    
    # read in standardized panel
    ref_panel <- cytofin_read_panel_info(panel_path = panel_path)
    
    # TO DO: are all of these values used somewhere???
    lineage_markers <- as.character(ref_panel$desc[ref_panel$Lineage == 1])
    functional_markers <- as.character(ref_panel$desc[ref_panel$Functional == 1])
    all_markers <- c(lineage_markers, functional_markers)
    
    # create transformation functions
    norm_1 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] - mean_ctr[m] + mean_uni[m]
      return(z)
    } #meanshift
    
    norm_2 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m] - mean(mean_ctr[m]) + mean(mean_uni[m])
      return(z)
    } #meanshift bulk
    
    norm_3 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- (z[m] - mean_ctr[m] + mean_uni[m])*sqrt(var_uni[m])/sqrt(var_ctr[m])
      return(z)
    } #variance
    
    norm_4 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- (z[m] - mean_ctr[m])*sqrt(var_uni[m])/sqrt(var_ctr[m]) + mean_uni[m]
      return(z)
    } #z-score
    
    norm_5 <- function(x) {
      y <- mean_uni[all_markers]
      z <- x
      m <- match(names(y), names(x))
      z[m] <- z[m]*lm(mean_uni[m]~mean_ctr[m])$coefficient[[2]]
      return(z)
    } #beadlike
    
    
    # Define a function for arsinh transformation (based on definition from Nikolay).
    # TO DO: Doesn't this just get rewritten later before it's ever used???
    asinhNik <- function(value) {
      value <- value - 1
      # TO DO: Why use a for-loop instead of vectorize???
      for(i in 1:length(value)) {
        if((value[i] < 0) | is.na(value[i])) {
          value[i] <- rnorm(1, mean = 0, sd = 0.01)
        }
      }  # can also use max(value, 0)
      value <- value / 5
      value <- asinh(value) 
      return(value)
    }
    
    # Reverses arsinh transformation with a cofactor of 5.
    rev_asinh <- function(x) {
      e.x <- exp(x)
      rev.x <- (e.x^2 - 1) / (2 * e.x)
      x <- rev.x * 5
      return(x)
    }
    
    # find the user-specified normalization function
    if (mode == "meanshift") {
      norm <- norm_1 
    } else if (mode == "meanshift_bulk") {
      norm <- norm_2
    } else if (mode == "variance") {
      norm <- norm_3
    } else if (mode == "z_score") {
      norm <- norm_4
    } else if (mode == "beadlike") {
      norm <- norm_5
    }
    
    # for each file being batch-normalized...
    for (i in 1:length(md$filename)) {
      #calculate adjustment parameters from control plate
      
      # find the control file corresponding to the same plate and cohort 
      # as the file being batch normalized
      filename_ctr <- 
        md_control$filename[
          which(
            (md_control$plate_number == md$plate_number[i]) & 
              (md_control$cohort == md$cohort[i])
          )
        ]
      fcs <- 
        flowCore::read.flowSet(
          file.path(input_path, paste0("homogenized_", filename_ctr)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      asinhNik <- flowCore::arcsinhTransform(a = 0, b = 0.2)
      colname <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_ctr <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
      mean_ctr <- apply(expr_ctr, 2, mean)
      mean_ctr_mean <- mean(mean_ctr)
      var_ctr <- apply(expr_ctr, 2, var)
      var_ctr_mean <- mean(var_ctr)
      markername <- flowCore::pData(flowCore::parameters(fcs[[1]]))$desc
      
      #normalize the control plate
      expr_ctr_norm <- t(flowCore::fsApply(fcs_asinh, function(x) apply(x, 1, norm), use.exprs=TRUE))
      mean_ctr_norm <- apply(expr_ctr_norm, 2, mean)
      var_ctr_norm <- apply(expr_ctr_norm, 2, var)
      mean_ctr_norm_mean <- mean(mean_ctr_norm)
      var_ctr_norm_mean <- mean(var_ctr_norm)
      
      # normalize the target plate
      ## before
      filename <- md$filename[i]
      fcs <- 
        flowCore::read.flowSet(
          file.path(input_path, paste0("homogenized_", filename)), 
          transformation = FALSE, 
          truncate_max_range = FALSE
        )
      
      # doesn't this overwrite the function definition above?
      asinhNik <- flowCore::arcsinhTransform(a = 0, b = 0.2)
      colname <- flowCore::colnames(fcs)
      tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
      fcs_asinh <- flowCore::transform(fcs, tlist)
      expr_b4norm <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
      mean_b4norm <- apply(expr_b4norm, 2, mean)
      var_b4norm <- apply(expr_b4norm, 2, var)
      mean_b4norm_mean <- mean(mean_b4norm)
      var_b4norm_mean <- mean(var_b4norm)
      
      ## after
      expr_norm <- t(flowCore::fsApply(fcs_asinh, function(x) apply(x, 1, norm), use.exprs = TRUE))
      mean_norm <- apply(expr_norm, 2, mean)
      var_norm <- apply(expr_norm, 2, var)
      mean_norm_mean <- mean(mean_norm)
      var_norm_mean <- mean(var_norm)
      fcs_norm <- flowCore::flowFrame(expr_norm)
      
      
      # normalization completed, reverse transformation
      tlist2 <- flowCore::transformList(from = colname, tfun = rev_asinh)
      fcs_asinh_rev <- flowCore::transform(fcs_norm, tlist2)
      flowCore::pData(flowCore::parameters(fcs_asinh_rev))$desc <- flowCore::pData(flowCore::parameters(fcs_asinh[[1]]))$desc
      fcs_name <- paste0(output_path,"normalized_",filename)
      flowCore::write.FCS(fcs_asinh_rev, fcs_name)
      
      # compare to validation
      # TO DO: What does "validation" mean in this context? should this be more modular?
      if (val_path != "none") {
        filename_val <- md$validation[i]
        fcs <- 
          flowCore::read.flowSet(
            paste0(val_path,filename_val), 
            transformation = FALSE, 
            truncate_max_range = FALSE
          )
        asinhNik <- flowCore::arcsinhTransform(a = 0, b = 0.2)
        colname <- flowCore::colnames(fcs)
        tlist <- flowCore::transformList(from = colname, tfun = asinhNik)
        fcs_asinh <- flowCore::transform(fcs, tlist)
        expr_val <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
        mean_val <- apply(expr_val, 2, mean)
        var_val <- apply(expr_val, 2, var)
        mean_val_mean <- mean(mean_val)
        var_val_mean <- mean(var_val)
      }
      
      
      
      # visualize antigen panel comparing universal vs plate
      # TO DO: these visualizations are a little confusing...
      # should we abstract them into their own new function?
      par(mfrow = c(2, 4))
      
      len <- length(mean_uni[all_markers])
      
      #expression (mean)
      #plot 1
      plot(
        mean_uni[all_markers], 
        col = 'red', 
        xlab = 'antigen', 
        ylab = 'universal expression (mean)', 
        xlim = c(0, len), 
        ylim = c(-5,10), 
        main = 'overall', 
        cex.main = 1
      )
      
      legend(1, 10, legend = c("universal"), col = c("red"), lty=1:2, cex = 0.8)
      
      #plot 2
      plot(mean_ctr[all_markers], col='cyan', xlab='antigen', ylab='control expression (mean)', xlim=c(0,len), ylim=c(-5,10), main=filename_ctr, cex.main=0.8)
      par(new=TRUE)
      plot(mean_ctr_norm[all_markers], col='blue', xlab='antigen', ylab='control expression (mean)', xlim=c(0, len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("blue", "cyan"), lty=1:2, cex=0.8)
      
      #plot 3
      plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='plate expression (mean)', xlim=c(0,len), ylim=c(-5,10), main=filename, cex.main=0.8)
      par(new=TRUE)
      plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='plate expression (mean)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("darkgreen", "green"), lty=1:2, cex=0.8)
      
      #plot 4
      if (val_path != "none"){
        plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_val[all_markers], col='purple', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized", "validation"),col=c("green", "darkgreen", "purple"), lty=1:2, cex=0.8)
      } else {
        plot(mean_b4norm[all_markers], col='green', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(mean_norm[all_markers], col='darkgreen', xlab='antigen', ylab='overlay expression (mean)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized"), col=c("green", "darkgreen"), lty=1:2, cex=0.8)
      }
      
      
      #expression (std)
      #plot 5
      plot(sqrt(var_uni[all_markers]), col='red', xlab='antigen', ylab='universal expression (std)', xlim=c(0,len), ylim=c(-5,10), main='overall', cex.main=1)
      legend(1,10,legend=c("universal"),col=c("red"), lty=1:2, cex=0.8)
      
      #plot 6
      plot(sqrt(var_ctr[all_markers]), col='cyan', xlab='antigen', ylab='control expression (std)', xlim=c(0,len), ylim=c(-5,10), main=filename_ctr, cex.main=0.8)
      par(new=TRUE)
      plot(sqrt(var_ctr_norm[all_markers]), col='blue', xlab='antigen', ylab='control expression (std)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("blue", "cyan"), lty=1:2, cex=0.8)
      
      #plot 7
      plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='plate expression (std)', xlim=c(0,len), ylim=c(-5,10), main=filename, cex.main=0.8)
      par(new=TRUE)
      plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='plate expression (std)', xlim=c(0,len), ylim=c(-5,10), cex.main=0.8)
      legend(1,10,legend=c("normalized", "original"),col=c("darkgreen", "green"), lty=1:2, cex=0.8)
      
      #plot 8
      if (val_path != "none"){
        plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(var_val[all_markers], col='purple', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized", "validation"), col=c("green", "darkgreen", "purple"), lty=1:2, cex=0.8)
      } else {
        plot(sqrt(var_b4norm[all_markers]), col='green', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        plot(sqrt(var_norm[all_markers]), col='darkgreen', xlab='antigen', ylab='overlay expression (std)', xlim=c(0,len), ylim=c(-5,10))
        par(new=TRUE)
        legend(1,10,legend=c("original", "normalized"), col=c("green", "darkgreen"), lty=1:2, cex=0.8)
      }
      
      
    }
  }
