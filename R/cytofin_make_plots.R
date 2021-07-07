#' Make diagnostic plots to evaluate CytofIn batch normalization
#' 
#' When given the output data structure from `cytofin_normalize` or `cytofin_normalize_nrs`, 
#' this function plots mean and variance plots for all normalized .fcs files and their 
#' associated anchors.
#'
#' @param normalization_result An output data.frame produced by the `cytofin_normalize` or
#' `cytofin_normalize_nrs` function. 
#' 
#' The following columns should be present: `filename`,
#' `cohort`, `plate_number`, `patient_id`, `condition`, `is_anchor`, `validation`,
#'  `universal_var`, `anchor_mean`, `anchor_var`, `mean_b4norm`, `var_b4norm`,
#'  `mean_norm`, `var_norm`, `mean_ctr_norm`, `var_ctr_norm`.
#'  
#'  @param which_rows A numeric vector indicating which rows of `normalization_result`
#'  (i.e. which .fcs files in the combined dataset) should be used for plotting. Defaults 
#'  to 1:nrow(normalization_result), which will make all possible plots.
#'
#' @param val_path The folder directory containing validation (i.e. bead-normalized)
#' .fcs files corresponding to the input .fcs files in the metadata table. (Optional).
#'
#' @return 8 diagnostic plots are made for each input .fcs file that was batch 
#' normalized (i.e. each .fcs file represented as a row in `normalization_result`). 
#' From left-to-right (and top-to-bottom), these plots represent the following: 
#'    1) The entry in the universal mean vector corresponding to each antigen in the 
#'    consensus antigen panel. X-axis: antigen index in the universal mean vector. 
#'    Y-axis: Arcsinh-transformed entry in the universal mean vector corresponding 
#'    to each antigen. 
#'    2) The mean (across all cells) antigen expression vector for the anchor 
#'    associated with each input .fcs file both before and after normalization. 
#'    X-axis: antigen index (as in plot 1). Y-axis: Mean antigen expression in the 
#'    anchor .fcs file.
#'    3) The mean (across all cells) antigen expression vector for each input 
#'    .fcs file both before and after normalization. 
#'    X-axis: antigen index (as in plot 1). Y-axis: Mean antigen expression in the 
#'    input .fcs file.
#'    4) The mean (across all cells) antigen expression vector for each "validation"
#'    (i.e. bead-normalized) .fcs file both before and after bead-normalization. 
#'    This plot can be used to compare CytofIn batch normalization with gold-
#'    standard approaches. If `val_path` is "none", this plot will be identical to 
#'    plot 3 (see above).  
#'    X-axis: antigen index (as in plot 1). Y-axis: Mean antigen expression in the 
#'    validation .fcs file.
#'    5) The entry in the universal standard deviation vector corresponding to each antigen in the 
#'    consensus antigen panel. X-axis: antigen index in the universal standard deviation vector. 
#'    Y-axis: Arcsinh-transformed entry in the universal standard deviation vector corresponding 
#'    to each antigen.
#'    6) The standard deviation (across all cells) antigen expression vector for the anchor 
#'    associated with each input .fcs file both before and after normalization. 
#'    X-axis: antigen index (as in plot 1). Y-axis: the standard deviation of all 
#'    antigen expression values in the anchor .fcs file.
#'    7) The standard deviation (across all cells) antigen expression vector for each input 
#'    .fcs file both before and after normalization. 
#'    X-axis: antigen index (as in plot 1). Y-axis: the standard deviation of all 
#'    antigen expression values in the input .fcs file.
#'    8) The standard deviation (across all cells) antigen expression vector for each "validation"
#'    (i.e. bead-normalized) .fcs file both before and after bead-normalization. 
#'    This plot can be used to compare CytofIn batch normalization with gold-
#'    standard approaches. If `val_path` is "none", this plot will be identical to 
#'    plot 3 (see above).  
#'    X-axis: antigen index (as in plot 1). Y-axis: Standard deviation of all 
#'    antigen expression values in the validation .fcs file.
#'
#' @export
#'
cytofin_make_plots <- 
  function(
    normalization_result, 
    which_rows = 1:nrow(normalization_result),
    val_path = "none"
  ) {
    
    # extract needed values from the normalization_result attributes
    all_markers <- attr(normalization_result, which = "all_markers")
    shift_factor <- attr(normalization_result, which = "shift_factor")
    scale_factor <- attr(normalization_result, which = "scale_factor")
    
    # filter out rows that we aren't interested in plotting
    normalization_result <- normalization_result[which_rows,]
    
    # for all rows in the normalization result
    for (i in 1:nrow(normalization_result)) {
      
      # extract needed values for the current file
      filename <- normalization_result$filename[[i]]
      plate_number <- normalization_result$plate_number[[i]]
      cohort <- normalization_result$cohort[[i]]
      universal_mean <- normalization_result$universal_mean[[i]]
      universal_var <- normalization_result$universal_var[[i]]
      anchor_mean <- normalization_result$anchor_mean[[i]]
      anchor_var <- normalization_result$anchor_var[[i]]
      mean_b4norm <- normalization_result$mean_b4norm[[i]]
      var_b4norm <- normalization_result$var_b4norm[[i]]
      mean_norm <- normalization_result$mean_norm[[i]]
      var_norm <- normalization_result$var_norm[[i]]
      anchor_mean_norm <- normalization_result$anchor_mean_norm[[i]]
      anchor_var_norm <- normalization_result$anchor_var_norm[[i]]
      
      # find name of the anchor that corresponds to each file
      md_control <- dplyr::filter(normalization_result, is_anchor == 1)
      
      filename_anchor <-
        md_control$filename[
          which(
            (md_control$plate_number == plate_number) &
              (md_control$cohort == cohort)
          )
        ]
      
      # read in validation .fcs file
      if (val_path != "none") {
        filename_val <- normalization_result$validation[i]
        fcs <-
          flowCore::read.flowSet(
            file.path(val_path, filename_val),
            transformation = FALSE,
            truncate_max_range = FALSE
          )
        
        # arcsinh-transform validation .fcs file
        asinh_transform <-
          flowCore::arcsinhTransform(a = shift_factor, b = scale_factor)
        col_names <- flowCore::colnames(fcs)
        tlist <- flowCore::transformList(from = col_names, tfun = asinh_transform)
        fcs_asinh <- flowCore::transform(fcs, tlist)
        expr_val <- flowCore::fsApply(fcs_asinh, flowCore::exprs)
        
        # find the mean and variance vector from the validation file
        mean_val <- apply(expr_val, 2, mean)
        var_val <- apply(expr_val, 2, var)
        
        # find the bulk mean and variance from the validation file
        mean_val_mean <- mean(mean_val)
        var_val_mean <- mean(var_val)
      }
      
      # make visualizations
      par(mfrow = c(2, 4))
      len <- length(universal_mean[all_markers])
      
      # expression (mean)
      # plot 1
      plot(
        universal_mean[all_markers],
        col = "red",
        xlab = "antigen",
        ylab = "universal expression (mean)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        main = "overall",
        cex.main = 1
      )
      
      legend(1, 10, legend = c("universal"), col = c("red"), lty = 1:2, cex = 0.8)
      
      # plot 2
      plot(
        anchor_mean[all_markers],
        col = "cyan",
        xlab = "antigen",
        ylab = "control expression (mean)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        main = filename_anchor,
        cex.main = 0.8
      )
      
      par(new = TRUE)
      plot(
        anchor_mean_norm[all_markers],
        col = "blue",
        xlab = "antigen",
        ylab = "control expression (mean)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        cex.main = 0.8
      )
      legend(1, 10, legend = c("normalized", "original"), col = c("blue", "cyan"), lty = 1:2, cex = 0.8)
      
      # plot 3
      plot(
        mean_b4norm[all_markers],
        col = "green",
        xlab = "antigen",
        ylab = "sample expression (mean)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        main = filename,
        cex.main = 0.8
      )
      par(new = TRUE)
      plot(
        mean_norm[all_markers],
        col = "darkgreen",
        xlab = "antigen",
        ylab = "sample expression (mean)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        cex.main = 0.8
      )
      legend(1, 10, legend = c("normalized", "original"), col = c("darkgreen", "green"), lty = 1:2, cex = 0.8)
      
      # plot 4
      if (val_path != "none") {
        plot(
          mean_b4norm[all_markers],
          col = "green",
          xlab = "antigen",
          ylab = "overlay expression (mean)",
          xlim = c(0, len),
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        plot(
          mean_norm[all_markers], 
          col = "darkgreen", 
          xlab = "antigen", 
          ylab = "overlay expression (mean)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        plot(
          mean_val[all_markers], 
          col = "purple", 
          xlab = "antigen", 
          ylab = "overlay expression (mean)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        par(new = TRUE)
        legend(1, 10, legend = c("original", "normalized", "validation"), col = c("green", "darkgreen", "purple"), lty = 1:2, cex = 0.8)
        
      } else {
        
        plot(
          mean_b4norm[all_markers], 
          col = "green", 
          xlab = "antigen", 
          ylab = "overlay expression (mean)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        par(new = TRUE)
        plot(mean_norm[all_markers], col = "darkgreen", xlab = "antigen", ylab = "overlay expression (mean)", xlim = c(0, len), ylim = c(-5, 10))
        par(new = TRUE)
        legend(1, 10, legend = c("normalized", "original"), col = c("darkgreen", "green"), lty = 1:2, cex = 0.8)
      }
      
      # expression (std)
      # plot 5
      plot(
        sqrt(universal_var[all_markers]),
        col = "red",
        xlab = "antigen",
        ylab = "universal expression (std)",
        xlim = c(0, len),
        ylim = c(-5, 10),
        main = "overall",
        cex.main = 1
      )
      legend(1, 10, legend = c("universal"), col = c("red"), lty = 1:2, cex = 0.8)
      
      # plot 6
      plot(sqrt(anchor_var[all_markers]), col = "cyan", xlab = "antigen", ylab = "control expression (std)", xlim = c(0, len), ylim = c(-5, 10), main = filename_anchor, cex.main = 0.8)
      par(new = TRUE)
      plot(sqrt(anchor_var_norm[all_markers]), col = "blue", xlab = "antigen", ylab = "control expression (std)", xlim = c(0, len), ylim = c(-5, 10), cex.main = 0.8)
      legend(1, 10, legend = c("normalized", "original"), col = c("blue", "cyan"), lty = 1:2, cex = 0.8)
      
      # plot 7
      plot(
        sqrt(var_b4norm[all_markers]), 
        col = "green", 
        xlab = "antigen", 
        ylab = "sample expression (std)", 
        xlim = c(0, len), 
        ylim = c(-5, 10), 
        main = filename, 
        cex.main = 0.8
      )
      
      par(new = TRUE)
      plot(
        sqrt(var_norm[all_markers]), 
        col = "darkgreen", 
        xlab = "antigen", 
        ylab = "sample expression (std)", 
        xlim = c(0, len), 
        ylim = c(-5, 10), 
        cex.main = 0.8
      )
      
      legend(1, 10, legend = c("normalized", "original"), col = c("darkgreen", "green"), lty = 1:2, cex = 0.8)
      
      # plot 8
      if (val_path != "none") {
        plot(
          sqrt(var_b4norm[all_markers]), 
          col = "green", 
          xlab = "antigen", 
          ylab = "overlay expression (std)",
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        plot(
          sqrt(var_norm[all_markers]), 
          col = "darkgreen", 
          xlab = "antigen", 
          ylab = "overlay expression (std)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        plot(
          var_val[all_markers], 
          col = "purple", 
          xlab = "antigen", 
          ylab = "overlay expression (std)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        legend(1, 10, legend = c("original", "normalized", "validation"), col = c("green", "darkgreen", "purple"), lty = 1:2, cex = 0.8)
        
      } else {
        plot(
          sqrt(var_b4norm[all_markers]), 
          col = "green", 
          xlab = "antigen", 
          ylab = "overlay expression (std)", 
          xlim = c(0, len), 
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        plot(
          sqrt(var_norm[all_markers]), 
          col = "darkgreen", 
          xlab = "antigen", 
          ylab = "overlay expression (std)", 
          xlim = c(0, len),
          ylim = c(-5, 10)
        )
        
        par(new = TRUE)
        legend(
          1, 
          10, 
          legend = c("original", "normalized"), 
          col = c("green", "darkgreen"), 
          lty = 1:2, cex = 0.8
        )
      }
    }
  }
