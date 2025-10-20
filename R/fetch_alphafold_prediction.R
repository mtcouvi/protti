#' Fetch AlphaFold prediction
#'
#' Fetches atom level data for AlphaFold predictions either for selected proteins or whole
#' organisms.
#'
#' @param uniprot_ids optional, a character vector of UniProt identifiers for which predictions
#' should be fetched. This argument is mutually exclusive to the \code{organism_name} argument.
#' @param organism_name optional, a character value providing the name of an organism for which
#' all available AlphaFold predictions should be retreived. The name should be the capitalised
#' scientific species name (e.g. "Homo sapiens"). **Note:** Some organisms contain a lot of
#' predictions which might take a considerable amount of time and memory to fetch. Therefore, you
#' should be sure that your system can handle fetching predictions for these organisms. This
#' argument is mutually exclusive to the \code{uniprot_ids} argument.
#' @param version a character value that specifies the alphafold version that should be used. This
#' is regularly updated by the database. We always try to make the current version the default version.
#' Available version can be found here: https://ftp.ebi.ac.uk/pub/databases/alphafold/
#' @param timeout a numeric value specifying the time in seconds until the download of an organism
#' archive times out. The default is 3600 seconds.
#' @param max_tries a numeric value that specifies the number of times the function tries to download
#' the data in case an error occurs. The default is 5. This only applies if `uniprot_ids` were provided.
#' @param return_data_frame a logical value that specifies if true, a data frame instead of a list
#' is returned. It is recommended to only use this if information for few proteins is retrieved.
#' Default is FALSE.
#' @param show_progress a logical value that specifies if true, a progress bar will be shown.
#' Default is TRUE.
#'
#' @return A list that contains atom level data for AlphaFold predictions. If return_data_frame is
#' TRUE, a data frame with this information is returned instead. The data frame contains the
#' following columns:
#'
#' * label_id: Uniquely identifies every atom in the prediction following the standardised
#' convention for mmCIF files.
#' * type_symbol: The code used to identify the atom species representing this atom type.
#' This code is the element symbol.
#' * label_atom_id: Uniquely identifies every atom for the given residue following the
#' standardised convention for mmCIF files.
#' * label_comp_id: A chemical identifier for the residue. This is the three- letter code
#' for the amino acid.
#' * label_asym_id: Chain identifier following the standardised convention for mmCIF files.
#' Since every prediction only contains one protein this is always "A".
#' * label_seq_id: Uniquely and sequentially identifies residues for each protein. The
#' numbering corresponds to the UniProt amino acid positions.
#' * x: The x coordinate of the atom.
#' * y: The y coordinate of the atom.
#' * z: The z coordinate of the atom.
#' * prediction_score: Contains the prediction score for each residue.
#' * auth_seq_id: Same as \code{label_seq_id}. But of type character.
#' * auth_comp_id: Same as \code{label_comp_id}.
#' * auth_asym_id: Same as \code{label_asym_id}.
#' * uniprot_id: The UniProt identifier of the predicted protein.
#' * score_quality: Score annotations.
#'
#' @import dplyr
#' @import progress
#' @import purrr
#' @import tidyr
#' @importFrom tibble tibble
#' @importFrom utils download.file untar
#' @importFrom readr read_tsv
#' @importFrom stringr str_replace_all str_detect
#' @importFrom curl has_internet
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @export
#'
#' @examples
#' \donttest{
#' alphafold <- fetch_alphafold_prediction(
#'   uniprot_ids = c("F4HVG8", "O15552"),
#'   return_data_frame = TRUE
#' )
#'
#' head(alphafold, n = 10)
#' }


fetch_alphafold_prediction <- function(uniprot_ids = NULL,
                                       organism_name = NULL,  # not implemented; kept for compatibility
                                       version = "v4",        # ignored (API provides version); kept for compatibility
                                       timeout = 3600,
                                       max_tries = 5,
                                       return_data_frame = FALSE,
                                       show_progress = TRUE) {

  if (!curl::has_internet()) {
    message("No internet connection.")
    return(invisible(NULL))
  }
  if (!missing(uniprot_ids) & !missing(organism_name)) {
    stop(strwrap("Please only provide either a list of UniProt identifiers or one organism name!",
                 prefix = "\n", initial = ""))
  }
  if (is.null(uniprot_ids) || length(uniprot_ids) == 0) {
    message("No UniProt IDs supplied.")
    return(invisible(NULL))
  }

  # sanitize IDs
  uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]

  # helper: GET JSON with retries
  .get_json <- function(url, tries = max_tries, timeout_sec = timeout) {
    for (i in seq_len(tries)) {
      resp <- tryCatch(
        httr::GET(url, httr::timeout(timeout_sec)),
        error = function(e) e
      )
      if (inherits(resp, "error")) next
      if (httr::http_error(resp)) {
        # backoff a bit on server-side errors
        if (httr::status_code(resp) >= 500) Sys.sleep(min(5 * i, 30))
        next
      }
      return(jsonlite::fromJSON(httr::content(resp, as = "text", encoding = "UTF-8"), simplifyVector = TRUE))
    }
    return(paste0("Client error: Failed to fetch JSON from ", url))
  }

  # helper: read a remote text file into tibble(X1 = lines)
  .read_lines_tbl <- function(file_url, tries = max_tries, timeout_sec = timeout) {
    for (i in seq_len(tries)) {
      con <- NULL
      out <- tryCatch({
        con <- curl::curl(file_url, open = "rb", handle = curl::new_handle(timeout = timeout_sec))
        lines <- readr::read_lines(con, progress = FALSE)
        tibble::tibble(X1 = lines)
      }, error = function(e) e, finally = if (!is.null(con)) close(con))
      if (inherits(out, "error")) {
        Sys.sleep(min(2 * i, 10))
        next
      }
      return(out)
    }
    return(paste0("Client error: Failed to download file: ", file_url))
  }

  # map IDs to preferred file URLs via the API
  api_urls <- setNames(
    paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot_ids),
    uniprot_ids
  )

  if (isTRUE(show_progress)) {
    pb <- progress::progress_bar$new(
      total = length(api_urls),
      format = "  Resolving AFDB file URLs [:bar] :current/:total (:percent) :eta"
    )
  }

  file_map <- lapply(names(api_urls), function(id) {
    url <- api_urls[[id]]
    meta <- .get_json(url)
    if (isTRUE(show_progress)) pb$tick()

    if (is.character(meta)) {
      # error string
      return(list(id = id, error = meta))
    }

    # API returns a list of prediction objects; pick F1 by default
    # (To use all fragments, drop the filtering line below.)
    preds <- meta
    if (is.data.frame(preds)) preds <- as.list(as.data.frame(t(preds), stringsAsFactors = FALSE))
    # If it's a list of predictions (most common case)
    if (is.list(preds) && length(preds) > 0) {
      # find F1 if present, else first
      pick <- NULL
      for (p in preds) {
        # common fields include 'id' like "AF-<ACC>-F1-model_v4"
        if (!is.null(p$id) && grepl("-F1-", p$id)) { pick <- p; break }
      }
      if (is.null(pick)) pick <- preds[[1]]

      # Prefer mmCIF; then bCIF; then PDB
      file_url <- pick$cifUrl %||% pick$bcifUrl %||% pick$pdbUrl
      if (is.null(file_url)) {
        return(list(id = id, error = "No downloadable model URL (cifUrl/bcifUrl/pdbUrl) in API response"))
      } else {
        return(list(id = id, file_url = file_url))
      }
    } else {
      return(list(id = id, error = "Unexpected API response format"))
    }
  })

  # split successes and errors
  errors <- Filter(function(x) !is.null(x$error), file_map)
  ok     <- Filter(function(x) !is.null(x$file_url), file_map)

  if (length(errors)) {
    error_table <- tibble::tibble(
      id    = vapply(errors, `[[`, "", "id"),
      error = vapply(errors, `[[`, "", "error")
    ) |> dplyr::distinct()
    message("Some IDs could not be resolved to file URLs:")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  if (!length(ok)) {
    message("No valid file URLs could be resolved from the API.")
    return(invisible(NULL))
  }

  # Now fetch & parse structures (as before)
  if (isTRUE(show_progress)) {
    pb2 <- progress::progress_bar$new(
      total = length(ok),
      format = "  Fetching AlphaFold predictions [:bar] :current/:total (:percent) :eta"
    )
  }

  query_result <- setNames(lapply(ok, function(rec) {
    q <- .read_lines_tbl(rec$file_url)
    if (isTRUE(show_progress)) pb2$tick()

    if (is.character(q)) {
      # pass through error string
      return(q)
    }

    # your original parsing, unchanged
    q %>%
      dplyr::filter(stringr::str_detect(X1, pattern = "^ATOM\\s+\\d|^HETATM\\s+\\d")) %>%
      dplyr::mutate(X2 = stringr::str_replace_all(X1, pattern = "\\s+", replacement = " ")) %>%
      tidyr::separate(
        X2,
        sep = " ",
        into = c(
          "x1","label_id","type_symbol","label_atom_id","x2","label_comp_id","label_asym_id",
          "entity_id","label_seq_id","x3","x","y","z","site_occupancy","prediction_score",
          "formal_charge","auth_seq_id","auth_comp_id","auth_asym_id","x4","pdb_model_number",
          "uniprot_id","x5","x6","x7"
        ),
        fill = "right", remove = TRUE
      ) %>%
      dplyr::select(-c(
        "X1","x1","x2","x3","x4","x5","x6","x7",
        "formal_charge","site_occupancy","entity_id","pdb_model_number"
      )) %>%
      dplyr::mutate(
        label_id       = suppressWarnings(as.numeric(.data$label_id)),
        label_seq_id   = suppressWarnings(as.numeric(.data$label_seq_id)),
        x              = suppressWarnings(as.numeric(.data$x)),
        y              = suppressWarnings(as.numeric(.data$y)),
        z              = suppressWarnings(as.numeric(.data$z)),
        prediction_score = suppressWarnings(as.numeric(.data$prediction_score)),
        auth_seq_id    = .data$auth_seq_id
      ) %>%
      dplyr::mutate(score_quality = dplyr::case_when(
        .data$prediction_score > 90 ~ "very_good",
        .data$prediction_score > 70 ~ "confident",
        .data$prediction_score > 50 ~ "low",
        .data$prediction_score <= 50 ~ "very_low",
        TRUE ~ NA_character_
      ))
  }), vapply(ok, `[[`, "", "id"))

  # report any per-file download/parse errors
  error_list <- purrr::keep(query_result, ~ is.character(.x))
  if (length(error_list) != 0) {
    error_table <- tibble::tibble(
      id = names(error_list),
      error = unlist(error_list)
    ) %>% dplyr::distinct()
    message("The following IDs have not been retrieved correctly.")
    message(paste0(utils::capture.output(error_table), collapse = "\n"))
  }

  # keep only successful data
  query_result <- purrr::keep(query_result, ~ !is.character(.x))
  if (length(query_result) == 0) {
    message("No valid information could be retrieved!")
    return(invisible(NULL))
  }

  if (isFALSE(return_data_frame)) {
    return(query_result)
  } else {
    return(purrr::list_rbind(query_result))
  }
}

