# S3 OOP interface.

# Constructors ------

#' Create a xena database object.
#'
#' @description Creates a XENA object on the top of a data frame.
#' @details The XENA object enhances functionality of tibbles to hold
#' mutation data. You may manipulate the XENA object with a plethora of
#' tidyverse tools in a very similar way as an ordinary tibble.
#' @param x a data frame.
#' @param ... extra arguments, currently none.
#' @return a xena database which inherits most of its methods from the
#' data.frame class.
#' @export

  xena <- function(x, ...) {

    if(!is.data.frame(x)) {

      stop('Please provide a data frame.', call. = FALSE)

    }

    structure(x,
              class = c('xena',
                        'tbl_df',
                        'tbl',
                        'data.frame'))

  }

# Class checker -----

#' Cehcks for the XENA class.
#'
#' @description Checks of the object is an instance of the XENA class.
#' @param x an object.
#' @return a logical value.
#' @export

  is_xena <- function(x) {

    any(class(x) == 'xena')

  }

# Searching in the XENA object -----

#' Query a local XENA database.
#'
#' @description Searches for samples carrying the specified mutation.
#' Lists all the samples with a mutation of the given gene, chromosome, with
#' the effect provided or mutation ID. Searches directly via the key value or
#' with a regular expression. The option 'detailed' returns a complete subset
#' of the xena_mut object. Otherwise, only the sample ids and vaf value are
#' returned.
#' @param x a XENA object.
#' @param key name of the database field to search in. Possible values are:
#' 'gene' (gene with the mutation), 'chrom' (chromosome), 'effect' (effect of
#' the mutation), 'mutation_id' (mutation ID).
#' @param value a key value or a vector for the direct search.
#' @param regex a regular expression to scan the key, ignored if value is
#' provided.
#' @param unique_samples logical, should duplicated samples be excluded from
#' the output?
#' @param join logical, should the output be merged into a single data frame by
#' the sample ID?
#' @param detailed logical, if TRUE a XENA object being a subset of the input
#' is returned.
#' @return a data frame or a XENA object, if detailed = TRUE.
#' @importFrom gseaTools search
#' @export search.xena
#' @export

  search.xena <- function(x,
                          key = c('gene', 'chrom', 'effect', 'mutation_id'),
                          value = NULL,
                          regex = NULL,
                          unique_samples = FALSE,
                          join = FALSE,
                          detailed = FALSE) {

    ## entry control

    stopifnot(xena::is_xena(x))

    key <- match.arg(key[1], c('gene', 'chrom', 'effect', 'mutation_id'))

    if(!is.null(value)) regex <- NULL

    stopifnot(is.logical(unique_samples))
    stopifnot(is.logical(join))
    stopifnot(is.logical(detailed))

    ## benchmarking

    start_time <- Sys.time()
    message(paste('Checking mutations for', length(value), key, 'values'))
    on.exit(message(paste('Elapsed:', Sys.time() - start_time)))

    ## quering

    if(!is.null(regex)) {

      x <- dplyr::filter(x, stringi::stri_detect(.data[[key]], regex = regex))

    }

    if(!is.null(value)) {

      x <- dplyr::filter(x, .data[[key]] %in% value)

    }

    ## searching the data base

    if(!detailed) {

      query <- plyr::dlply(x, key, function(x) x[c('sample_id', 'dna_vaf')])

      query <- purrr::map2(query,
                           names(query),
                           ~rlang::set_names(.x, c('sample_id', .y)))

    } else {

      query <- plyr::dlply(x, key)

    }

    if(unique_samples) {

      query <- purrr::map(query,
                          ~dplyr::filter(.x, !duplicated(.x)))


    }

    if(!join){

      purrr::map(query, tibble::as_tibble)

    } else {

      query <- purrr::reduce(query, dplyr::full_join, by = 'sample_id')

      query <- tibble::as_tibble(query)

      purrr::map_dfc(query, ~ifelse(is.na(.x), 0, .x))

    }

  }

# XENA object: counting the feature of interest -------

#' Count unique occurrences of a feature in XENA object.
#'
#' @description Counts occurrences of unique categories of a variable in
#' a XENA object.
#' @details Functions in a similar way as the genuine `count()` function from
#' the dplyr package but handles also 'fuzzy' categories such as mutation
#' effects with multiple categories per sample seperated by semicolon.
#' A sample is counted once even if it e.g. has two missense mutations.
#' @return a data frame with counts, frequencies and percentages for each
#' category of a variable.
#' @param x a XENA class object.
#' @param ft a variable to count.
#' @param sort logical, if TRUE, will show the largest groups at the top.
#' @param name The name of the new column in the output.
#' If omitted, it will default to n.
#' @importFrom dplyr count
#' @export count.xena
#' @export

  count.xena <- function(x,
                         ft,
                         sort = FALSE,
                         name = NULL) {

    ## entry control -----

    stopifnot(is_xena(x))

    ft_catch <- rlang::enexpr(ft)

    ft_name <- rlang::as_string(ft_catch)

    if(!ft_name %in% names(x)) {

      stop('The ft variable is absent from the XENA object.',
           call. = FALSE)

    }

    ## counting -------

    n_total <- length(unique(x$sample_id))

    cats <- stringi::stri_split_fixed(x[[ft_name]],
                                      pattern = ';')

    cats <- unique(unlist(cats))

    cat_lst <-
      purrr::map(rlang::set_names(cats, cats),
                 ~dplyr::filter(x,
                                stringi::stri_detect(.data[[ft_name]],
                                                     fixed = .x)))

    cat_lst <- purrr::map(cat_lst,
                          ~dplyr::filter(.x, !duplicated(sample_id)))

    cat_counts <- purrr::map_dbl(cat_lst, nrow)

    if(is.null(name)) name <- 'n'

    cat_counts <- tibble::tibble(!!ft_name := names(cat_counts),
                                 !!name := cat_counts)

    dplyr::mutate(cat_counts,
                  n_total = n_total,
                  frequency = .data[[name]]/n_total,
                  percentage = .data[[name]]/n_total * 100)

  }

# END --------
