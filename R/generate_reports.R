#' Generate differential expression analysis reports
#' @export
#'
generate_DEA_reports <- function(lfqdata, GRP2, prot_annot, ZIPDIR) {
  # Compute all possible 2 Grps to avoid specifying reference.


  levels  <- lfqdata$factors() |>
    dplyr::select(Group_ = Group_,  control = starts_with("control", ignore.case = TRUE)) |>
    dplyr::distinct()
  logger::log_info("levels : ", paste(levels, collapse = " "))
  if(! length(levels$Group_) > 1){
    logger::log_error("not enough group levels_ to make comparisons.")
  }

  # Contrasts are already defined
  if(!is.null(GRP2$pop$Contrasts)) {
    logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = " "))
    lfqwork <- lfqdata$get_copy()
    lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)

    grp2 <- prolfquapp::make_DEA_report(lfqwork,
                                       prot_annot,
                                       GRP2)

    fname <- paste0("DE_Groups_vs_Controls")
    qcname <- paste0("QC_Groups_vs_Controls")
    outpath <- file.path( ZIPDIR, fname)

    logger::log_info("writing into : ", outpath, " <<<<")
    prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")

    bb <- grp2$RES$transformedlfqData
    if(sum(!grepl("^control",bb$config$table$factorKeys(), ignore.case = TRUE))  > 1) {
      prolfquapp::writeLinesPaired(bb, outpath)
    } else{
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }
    return()
  }

  ## Generate contrasts from dataset
  if(!is.null(levels$control)) {
    Contrasts <- character()
    Names <- character()
    for (i in 1:nrow(levels)) {
      for (j in 1:nrow(levels)) {
        if (i != j) {
          if(levels$control[j] == "C"){
            cat(levels$Group_[i], levels$Group_[j], "\n")
            Contrasts <- c(Contrasts, paste0("Group_",levels$Group_[i], " - ", "Group_",levels$Group_[j]))
            Names <- c(Names, paste0(levels$Group_[i], "_vs_", levels$Group_[j]))
          }
        }
      }
    }

    names(Contrasts) <- Names
    GRP2$pop$Contrasts <- Contrasts

    logger::log_info("CONTRAST : ", paste( GRP2$pop$Contrasts, collapse = " "))
    lfqwork <- lfqdata$get_copy()
    lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ %in% levels$Group_)

    grp2 <- prolfquapp::make_DEA_report(lfqwork,
                                       prot_annot,
                                       GRP2)

    fname <- paste0("DE_Groups_vs_Controls")
    qcname <- paste0("QC_Groups_vs_Controls")
    outpath <- file.path( ZIPDIR, fname)

    logger::log_info("writing into : ", outpath, " <<<<")
    prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname)
    prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")

    bb <- grp2$RES$transformedlfqData
    if(sum(!grepl("^control",bb$config$table$factorKeys(), ignore.case = TRUE))  > 1) {
      prolfquapp::writeLinesPaired(bb, outpath)
    } else{
      pl <- bb$get_Plotter()
      pl$write_boxplots(outpath)
    }

  } else {
    # create all possible 2 grp comparisons
    for (i in seq_along(levels$Group_)) {
      for (j in seq_along(levels$Group_)) {
        if (i != j) {
          cat(levels$Group_[i], levels$Group_[j], "\n")
          GRP2$pop$Contrasts <- paste0("Group_",levels$Group_[i], " - ", "Group_",levels$Group_[j])
          names(GRP2$pop$Contrasts) <- paste0(levels$Group_[i], "_vs_", levels$Group_[j])
          logger::log_info("CONTRAST : ", GRP2$pop$Contrasts, collapse = " ")
          lfqwork <- lfqdata$get_copy()
          lfqwork$data <- lfqdata$data |> dplyr::filter(.data$Group_ == levels$Group_[i] | .data$Group_ == levels$Group_[j])
          grp2 <- prolfquapp::make_DEA_report(lfqwork,
                                             prot_annot,
                                             GRP2)

          fname <- paste0("Group_" , levels$Group_[i], "_vs_", levels$Group_[j])
          qcname <- paste0("QC_" , levels$Group_[i], "_vs_", levels$Group_[j])
          outpath <- file.path( ZIPDIR, fname)
          logger::log_info("writing into : ", outpath, " <<<<")

          prolfquapp::write_DEA(grp2, outpath = outpath, xlsxname = fname)
          prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = fname)
          prolfquapp::render_DEA(grp2, outpath = outpath, htmlname = qcname, markdown = "_DiffExpQC.Rmd")
          bb <- grp2$RES$transformedlfqData
          if(sum(!grepl("^control",bb$config$table$factorKeys(), ignore.case = TRUE)) > 1) {
            prolfquapp::writeLinesPaired(bb, outpath)
          } else{
            pl <- bb$get_Plotter()
            pl$write_boxplots(outpath)
          }
        }
      }
    }
  }
}

