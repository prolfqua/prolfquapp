
.find_cleavage_sites <- function(sequence, pattern = "(K|R)(?!P|$)") {
  # Convert sequence to uppercase
  sequence <- toupper(sequence)
  # Find positions of tryptic sites
  positions <- stringr::str_locate_all(sequence, pattern)[[1]][,"end"] # gregexpr(pattern, sequence, perl = TRUE)[[1]]
  # stri_locate_all_regex(sequence, pattern)
  # cleavage_sites <- unlist(positions)
  return(positions)
}

.compute_peptide_lengths <- function(sequence, cleavage_sites) {
  # Compute the lengths of the peptides
  start_positions <- c(0, cleavage_sites)
  end_positions <- c(cleavage_sites, nchar(sequence))
  peptide_lengths <- end_positions - start_positions
  return(peptide_lengths)
}


#' Compute number of tryptic peptides
#' @export
#' @examples
#' # example code
#'
#' sequence <- "MKGLPRAKSHGSTGWGKRKRNKPK"
#' nr_tryptic_peptides(sequence, min_length=5)
#'
nr_tryptic_peptides <- function(sequence, min_length = 6, max_length = 30){
  peptide_lengths <- .compute_peptide_lengths(sequence, .find_cleavage_sites(sequence))
  res <- sum(peptide_lengths >= min_length & peptide_lengths < max_length)
  return(res)
}


#' extract gene names from uniprot 1sp fasta.headers
#' @param fasta.headers vector with
#' @export
extract_GN <- function(fasta.headers){
  res <- vector(mode = "character", length = length(fasta.headers))
  for (i in seq_along(fasta.headers)) {
    header <- fasta.headers[i]
    res[i] <- if (grepl(".+ GN=(.+) PE=.+", header)) {
      gsub(".+ GN=(.+) PE=.*","\\1",header)} else {""}
  }
  return(res)
}



.getSequences <- function(x){

  fasta_str <- ">sp|A0A385XJL2|YGDT_ECOLI Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1
MLSTESWDNCEKPPLLFPFTALTCDETPVFSGSVLNLVAHSVDKYGIG
>sp|A5A615|YNCL_ECOLI Uncharacterized protein YncL OS=Escherichia coli (strain K12) OX=83333 GN=yncL PE=1 SV=1
MNVSSRTVVLINFFAAVGLFTLISMRFGWFI
>sp|P03018|UVRD_ECOLI DNA helicase II OS=Escherichia coli (strain K12) OX=83333 GN=uvrD PE=1 SV=1
MDVSYLLDSLNDKQREAVAAPRSNLLVLAGAGSGKTRVLVHRIAWLMSVENCSPYSIMAV
>sp|P04982|RBSD_ECOLI D-ribose pyranase OS=Escherichia coli (strain K12) OX=83333 GN=rbsD PE=1 SV=3
MKKGTVLNSDISSVISRLGHTDTLVVCDAGLPIPKSTTRIDMALTQGVPSFMQVLGVVTN
>sp|P04994|EX7L_ECOLI Exodeoxyribonuclease 7 large subunit OS=Escherichia coli (strain K12) OX=83333 GN=xseA PE=1 SV=2
MLPSQSPAIFTVSRLNQTVRLLLEHEMGQVWISGEISNFTQPASGHWYFTLKDDTAQVRC
>zz|Y-FGCZCont00001|  zz_FGCZCont0000_P61626_LYSC_HUMAN blastpHomologue_5.0e-107
MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRA
>zz|Y-FGCZCont00002|  zz_FGCZCont0001_P02534_K1M1_SHEEP blastpHomologue_0.0
SFNFCLPNLSFRSSCSSRPCVPSSCCGTTLPGACNIPANVGSCNWFCEGSFDGNEKETMQ
>REV_sp|Q13515|BFSP2_HUMAN Phakinin OS=Homo sapiens OX=9606 GN=BFSP2 PE=1 SV=1
GSEERDLLAHYSAVDKQLQCKRALLHAREQQQQEAEARIERLEAELRGVVAGLNQLEMDH
>REV_sp|Q14183|DOC2A_HUMAN Double C2-like domain-containing protein alpha OS=Homo sapiens OX=9606 GN=DOC2A PE=1 SV=5
ASSLAGAAPPLESTLTHWRELAADPQQLCDSWHKRAEGRAGPGLSVGGIFDNSKGIDYDW
>REV_tr|A0A075B6W8|A0A075B6W8_HUMAN T cell receptor alpha joining 17 (Fragment) OS=Homo sapiens OX=9606 GN=TRAJ17 PE=4 SV=1
PKVLVRTGGGFTLKNGAAKIX
>REV_sp|A0A385XJL2|YGDT_ECOLI Protein YgdT OS=Escherichia coli (strain K12) OX=83333 GN=ygdT PE=4 SV=1
GIGYKDVSHAVLNLVSGSFVPTEDCTLATFPFLLPPKECNDWSETSLM
"
  return(fasta_str)
}

#' get_annot_from_fasta
#'
#' @export
#' @examples
#'
#' fasta_conn <- textConnection(prolfquapp:::.getSequences())
#' testthat::expect_error(prolfquapp::get_annot_from_fasta(fasta_conn, pattern_decoys = "" ))
#' close(fasta_conn)
#' fasta_conn <- textConnection(prolfquapp:::.getSequences())
#' prolfquapp::get_annot_from_fasta(fasta_conn, pattern_decoys = "^REV_|^rev" )
#' close(fasta_conn)
#'
get_annot_from_fasta <- function(
    fasta.files,
    pattern_decoys = "^REV_|^rev_",
    isUniprot = TRUE,
    min_length = 7,
    max_length = 30,
    include_seq = FALSE) {
  fasta <- list()

  if ("connection" %in% class(fasta.files) ) {
    fasta <- seqinr::read.fasta(file = fasta.files, as.string = TRUE, seqtype = "AA")
  } else {
    for (fasta.file in fasta.files) {
      logger::log_info("get_annot : ", fasta.file)
      x <- seqinr::read.fasta(file = fasta.file, as.string = TRUE, seqtype = "AA")
      fasta <- c(fasta, x)
    }
  }

  logger::log_info("get_annot : finished reading")
  fasta_annot <- data.frame(annot = vapply(fasta, seqinr::getAnnot, ""), sequence = as.character(fasta))
  fasta_annot <- fasta_annot |> tidyr::separate(.data$annot,
                                                c("fasta.id","fasta.header"),
                                                sep = "\\s", extra = "merge")
  fasta_annot <- fasta_annot |> dplyr::mutate(fasta.id = stringr::str_replace_all(.data$fasta.id, "^>", "") )

  logger::log_info("get_annot : extract headers")

  logger::log_info("get_annot : all seq : ", nrow(fasta_annot))
  if (!is.null(pattern_decoys) && pattern_decoys != "") {
    logger::log_info("removing decoy sequences usin patter : ", pattern_decoys)
    pcdecoy <- mean(grepl(pattern_decoys, fasta_annot$fasta.id))
    if (pcdecoy < 0.1) {
      logger::log_warn("Only ", pcdecoy, " found using pattern : ", pattern_decoys)
      logger::log_warn("Please specify empty string if no decoy's in fasta.")
      warning("no decoys found")
    }
    fasta_annot <- fasta_annot |> dplyr::filter( !grepl(pattern_decoys, .data$fasta.id))
    logger::log_info("get_annot nr seq after decoy removal: ", nrow(fasta_annot))
  }

  logger::log_info("get_annot : isUniprot : ", isUniprot)
  if (isUniprot) {
    fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = gsub(".+\\|(.+)\\|.*","\\1", .data$fasta.id) )
  } else {
    fasta_annot <- fasta_annot |> dplyr::mutate(proteinname = .data$fasta.id )
  }
  if (any(duplicated(fasta_annot$proteinname))) {
    logger::log_error("there are duplicated protein ID's , mean duplicate :", mean(duplicated(fasta_annot$proteinname)))
    logger::log_error("the pattern_decoys is : [", pattern_decoys, "]")
    stop("wrong decoys pattern.", pattern_decoys, "\n")
  }

  if (mean(grepl(".+ GN=(.+) PE=.+",fasta_annot$fasta.header)) > 0.5) {
    fasta_annot <- fasta_annot |> dplyr::mutate(gene_name = extract_GN(.data$fasta.header))
    logger::log_info("get_annot : extracted gene names")
  }

  # remove duplicated id's
  # fasta_annot <- fasta_annot[!duplicated(fasta_annot$proteinname),]

  fasta_annot$protein_length <- vapply(fasta_annot$sequence, nchar, 0)
  logger::log_info("get_annot : protein length")

  fasta_annot$nr_tryptic_peptides <- vapply(fasta_annot$sequence, nr_tryptic_peptides, 0, min_length = min_length, max_length = max_length)
  logger::log_info("get_annot : nr of tryptic peptides per protein computed.")

  if ( !include_seq ) {
    fasta_annot$sequence <- NULL
  }
  return(fasta_annot)
}

