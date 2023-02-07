# Loading and management of a XENA database.

# Loading -----

#' Load a XENA data base.
#'
#' @description Loads a XENA database from the disc and converts it to a
#' tidyverse - compatible format.
#' @details The XENA object enhances functionality of tibbles to hold
#' mutation data. You may manipulate the XENA object with a plethora of
#' tidyverse tools in a very similar way as an ordinary tibble.
#' @param path path to the XENA database file; accepts .tar or .gz archives.
#' @return a XENA object, with a data frame interface.
#' @export

  load_xena <- function(path) {

    xena_db <- readr::read_tsv(path)

    xena_db <- dplyr::mutate(xena_db,
                             sample_id = Sample_ID,
                             mutation_id = paste(chrom,
                                                 start,
                                                 end,
                                                 ref,
                                                 alt,
                                                 sep = '_'))

    xena_db <- dplyr::select(xena_db, - Sample_ID)

    xena::xena(xena_db)

  }

# Mutation counting ------

#' Count mutations per sample.
#'
#' @description Counts all mutations per sample provided in the XENA class
#' database.
#' @param xena_db a XENA database object.
#' @param ... extra arguments, currently none.
#' @return a data frame listing the samples and the number of mutations per
#' sample.
#' @export

  count_mutations <- function(xena_db, ...) {

    ## entry control

    if(!xena::is_xena(xena_db)) {

      stop('Please provide a valid XENA database object.', call. = FALSE)

    }

    ## counting

    xena_counts <- plyr::dlply(xena_db,
                               'sample_id',
                               dplyr::filter,
                               !duplicated(mutation_id))

    purrr::map2_dfr(xena_counts,
                    names(xena_counts),
                    ~tibble::tibble(sample_id = .y,
                                    mut_count = nrow(.x)))

  }

# Mutation table/matrix ------

#' List mutated genes.
#'
#' @description Creates a matrix or table with the mutated genes for each
#' sample (0: wt, 1: a mutation listed in the database). Specificity towards
#' specific mutations (e.g. missense ones) and their quantity (Variant Allele
#' Frequency) can be acheived by simple filtering of the input xena_bd object
#' (inherits the transformation interface from data.frame).
#' @param xena_db a XENA database object.
#' @param as_matrix logical, should the output be in a matrix form?
#' @param sparse logical or NULL, specifying if the result should be sparse or
#' not. By default, it is made sparse when more than half of the entries are 0.
#' @param .parallel logical, should the operation be run in parallel? Runs via
#' the foreach backend, the user need to register the process prior to calling
#' the function.
#' If FALSE, a data.frame is returned (slower).
#' @param ... extra arguments, currently none.
#' @return a (sparse) matrix or a data frame.
#' @export

  tab_mutations <- function(xena_db,
                            as_matrix = TRUE,
                            sparse = NULL,
                            .parallel = FALSE, ...) {

    ## entry control

    if(!xena::is_xena(xena_db)) {

      stop('Please provide a valid XENA database object.', call. = FALSE)

    }

    stopifnot(is.logical(as_matrix))
    stopifnot(is.logical(.parallel))

    ## benchmarking

    start_time <- Sys.time()
    message(paste('Mutation tabulation for n =', nrow(xena_db), 'records.'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## re-shaping the table

    mut_tbl <- plyr::dlply(xena_db,
                           'gene',
                           .parallel = .parallel,
                           dplyr::filter,
                           !duplicated(sample_id))

    ## the simple reduction with merging won't work because of the computation
    ## time. Instead, an iterative approach using a matrix

    av_samples <- unique(xena_db$sample_id)
    av_genes <- names(mut_tbl)

    mut_mtx <- Matrix::Matrix(0,
                              nrow = length(av_samples),
                              ncol = length(av_genes),
                              dimnames = list(av_samples, av_genes),
                              sparse = sparse)

    for(i in names(mut_tbl)) {

      mut_record <- mut_tbl[[i]]

      mut_mtx[mut_record$sample_id, i] <- 1

    }

    if(as_matrix) return(mut_mtx)

    mut_mtx <- as.matrix(mut_mtx)

    mut_mtx <- tibble::rownames_to_column(data.frame(mut_mtx), 'sample_id')

    tibble::as_tibble(mut_mtx)

  }

# Mutation frequency -----

#' Gen frequency of mutated genes.
#'
#' @description Calculates frequencies of the mutated genes in the database.
#' @param xena_db a XENA database object.
#' @param .parallel logical, should the computation be run in parallel?
#' @param ... extra arguments, currently none.
#' @return a data frame with the frequency and percent of mutations.
#' @export

  freq_mutations <- function(xena_db, .parallel = FALSE, ...) {

    ## entry control

    stopifnot(xena::is_xena(xena_db))

    stopifnot(is.logical(.parallel))

    ## mutation matrix

    mut_mat <- xena::tab_mutations(xena_db,
                                   as_matrix = TRUE,
                                   sparse = TRUE,
                                   .parallel = FALSE)

    av_genes <- colnames(mut_mat)

    ## frequencies

    start_time <- Sys.time()
    message('Mutation frequencies')
    on.exit(paste('Elapsed:', Sys.time() - start_time))

    if(.parallel) {

      future::plan('multisession')

      freqs <- furrr::future_map_dbl(av_genes,
                                     ~mean(mut_mat[, .x], na.rm = TRUE))

      future::plan('sequential')

    } else {

      freqs <- purrr::map_dbl(av_genes, ~mean(mut_mat[, .x], na.rm = TRUE))

    }

    tibble::tibble(gene = av_genes,
                   frequency = freqs,
                   percent = freqs * 100)

  }

# END -----
