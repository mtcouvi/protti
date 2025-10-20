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
                                       organism_name = NULL,   # kept for compatibility; not used
                                       version = "v4",         # kept for compatibility; API provides current
                                       timeout = 3600,
                                       max_tries = 5,
                                       return_data_frame = FALSE,
                                       show_progress = TRUE,
                                       include_all_fragments = FALSE,
                                       prefer = c("cif","bcif","pdb")) {
  # ---- quick guards ----
  if (!curl::has_internet()) {
    message("No internet connection."); return(invisible(NULL))
  }
  if (!missing(uniprot_ids) & !missing(organism_name)) {
    stop(strwrap("Please only provide either a list of UniProt identifiers or one organism name!",
                 prefix = "\n", initial = ""))
  }
  if (is.null(uniprot_ids) || length(uniprot_ids) == 0) {
    message("No UniProt IDs supplied."); return(invisible(NULL))
  }
  uniprot_ids <- uniprot_ids[!is.na(uniprot_ids)]
  prefer <- match.arg(prefer)

  # ---- helpers ----
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  .get_json <- function(url, tries = max_tries, timeout_sec = timeout) {
    for (i in seq_len(tries)) {
      resp <- tryCatch(httr::GET(url, httr::timeout(timeout_sec)), error = identity)
      if (inherits(resp, "error")) { Sys.sleep(min(2 * i, 10)); next }
      if (httr::http_error(resp))  { if (httr::status_code(resp) >= 500) Sys.sleep(min(5 * i, 30)); next }
      txt <- httr::content(resp, as = "text", encoding = "UTF-8")
      return(tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) paste0("Client error: ", e$message)))
    }
    paste0("Client error: Failed to fetch JSON from ", url)
  }

  .read_lines_tbl <- function(file_url, tries = max_tries, timeout_sec = timeout) {
    for (i in seq_len(tries)) {
      con <- NULL
      out <- tryCatch({
        con <- curl::curl(file_url, open = "rb", handle = curl::new_handle(timeout = timeout_sec))
        tibble::tibble(X1 = readr::read_lines(con, progress = FALSE))
      }, error = identity, finally = if (!is.null(con)) close(con))
      if (!inherits(out, "error")) return(out)
      Sys.sleep(min(2 * i, 10))
    }
    paste0("Client error: Failed to download file: ", file_url)
  }

  .parse_mmcif_atom_site <- function(lines) {
    # locate loop_ with _atom_site.* tags
    is_loop <- trimws(lines) == "loop_"
    loop_idxs <- which(is_loop)
    start_tags <- end_tags <- start_rows <- end_rows <- NULL
    for (li in loop_idxs) {
      i <- li + 1
      tags <- character()
      while (i <= length(lines) && grepl("^_", lines[i])) {
        if (grepl("^_atom_site\\.", lines[i])) tags <- c(tags, lines[i])
        i <- i + 1
      }
      if (length(tags)) {
        r <- i
        while (r <= length(lines) &&
               nzchar(trimws(lines[r])) &&
               !grepl("^(loop_|data_|#)", trimws(lines[r])) &&
               !grepl("^_", lines[r])) {
          r <- r + 1
        }
        start_tags <- li + 1; end_tags <- i - 1
        start_rows <- i;      end_rows <- r - 1
        break
      }
    }
    if (is.null(start_rows)) stop("Could not locate _atom_site loop in mmCIF.")

    tag_lines <- lines[start_tags:end_tags]
    colnames <- sub("^_atom_site\\.", "", trimws(tag_lines))
    dat_txt <- paste0(lines[start_rows:end_rows], collapse = "\n")

    df <- readr::read_table2(dat_txt, col_names = colnames, na = c(".", "?"),
                             progress = FALSE, show_col_types = FALSE)

    want <- intersect(c(
      "id","type_symbol","label_atom_id","label_comp_id",
      "label_asym_id","label_entity_id","label_seq_id",
      "Cartn_x","Cartn_y","Cartn_z","B_iso_or_equiv",
      "auth_seq_id","auth_comp_id","auth_asym_id"
    ), names(df))
    out <- df[, want, drop = FALSE]

    num_cols <- intersect(c("id","label_seq_id","Cartn_x","Cartn_y","Cartn_z","B_iso_or_equiv","auth_seq_id"), names(out))
    out[num_cols] <- lapply(out[num_cols], function(x) suppressWarnings(as.numeric(x)))

    dplyr::transmute(
      out,
      label_id         = .data$id,
      type_symbol      = .data$type_symbol,
      label_atom_id    = .data$label_atom_id,
      label_comp_id    = .data$label_comp_id,
      label_asym_id    = .data$label_asym_id,
      label_seq_id     = .data$label_seq_id,
      x                = .data$Cartn_x,
      y                = .data$Cartn_y,
      z                = .data$Cartn_z,
      prediction_score = .data$B_iso_or_equiv,     # AF pLDDT lives here
      auth_seq_id      = .data$auth_seq_id,
      auth_comp_id     = .data$auth_comp_id,
      auth_asym_id     = .data$auth_asym_id
    ) |>
      dplyr::mutate(
        score_quality = dplyr::case_when(
          .data$prediction_score > 90 ~ "very_good",
          .data$prediction_score > 70 ~ "confident",
          .data$prediction_score > 50 ~ "low",
          is.na(.data$prediction_score) ~ NA_character_,
          TRUE ~ "very_low"
        )
      )
  }

  .pick_url <- function(pred) {
    # prefer user choice, fallback through others
    if (prefer == "cif")  return(pred$cifUrl %||% pred$bcifUrl %||% pred$pdbUrl)
    if (prefer == "bcif") return(pred$bcifUrl %||% pred$cifUrl %||% pred$pdbUrl)
    if (prefer == "pdb")  return(pred$pdbUrl %||% pred$cifUrl %||% pred$bcifUrl)
  }

  # ---- resolve file URLs from API ----
  api_urls <- setNames(paste0("https://alphafold.ebi.ac.uk/api/prediction/", uniprot_ids), uniprot_ids)

  if (isTRUE(show_progress)) {
    pb <- progress::progress_bar$new(total = length(api_urls),
      format = "  Resolving AFDB file URLs [:bar] :current/:total (:percent) :eta")
  }

  # for each UniProt: collect 1 (F1) or many (all fragments) file URLs
  resolved <- lapply(names(api_urls), function(id) {
    meta <- .get_json(api_urls[[id]])
    if (isTRUE(show_progress)) pb$tick()

    if (is.character(meta)) return(list(id = id, errors = meta, files = NULL))

    preds <- meta
    # normalize to list-of-lists
    if (is.data.frame(preds)) preds <- split(preds, seq_len(nrow(preds)))
    if (!length(preds)) return(list(id = id, errors = "No predictions", files = NULL))

    # choose subset
    if (!include_all_fragments) {
      # prefer F1 (id like "AF-<ACC>-F1-model_vX")
      pick <- NULL
      for (p in preds) if (!is.null(p$id) && grepl("-F1-", p$id)) { pick <- p; break }
      if (is.null(pick)) pick <- preds[[1]]
      preds <- list(pick)
    }

    # build files
    files <- lapply(preds, function(p) {
      file_url <- .pick_url(p)
      if (is.null(file_url)) return(NULL)
      frag <- NA_character_
      if (!is.null(p$id)) {
        m <- regexpr("-F\\d+-", p$id)
        if (m > 0) frag <- sub("^-|-$", "", substring(p$id, m, m + attr(m, "match.length") - 1))
      }
      list(file_url = file_url,
           fragment = frag,
           model_id = p$id %||% NA_character_,
           model_version = p$modelVersion %||% NA_character_,
           uniprot = id)
    })
    files <- Filter(Negate(is.null), files)
    if (!length(files)) return(list(id = id, errors = "No downloadable model URL in API response", files = NULL))
    list(id = id, errors = NULL, files = files)
  })

  # collect errors & the full download plan
  api_errors <- Filter(function(x) !is.null(x$errors), resolved)
  if (length(api_errors)) {
    err_tbl <- tibble::tibble(
      id    = vapply(api_errors, `[[`, "", "id"),
      error = vapply(api_errors, `[[`, "", "errors")
    ) |> dplyr::distinct()
    message("Some IDs could not be resolved to file URLs:")
    message(paste0(utils::capture.output(err_tbl), collapse = "\n"))
  }

  plan <- unlist(lapply(resolved, function(x) x$files), recursive = FALSE)
  if (!length(plan)) {
    message("No valid file URLs could be resolved from the API.")
    return(invisible(NULL))
  }

  # ---- fetch & parse structures ----
  if (isTRUE(show_progress)) {
    pb2 <- progress::progress_bar$new(total = length(plan),
      format = "  Fetching AlphaFold predictions [:bar] :current/:total (:percent) :eta")
  }

  fetched <- lapply(plan, function(rec) {
    tbl <- .read_lines_tbl(rec$file_url)
    if (isTRUE(show_progress)) pb2$tick()
    if (is.character(tbl)) return(list(key = rec$uniprot, error = tbl, df = NULL, meta = rec))

    # robust mmCIF parse
    parsed <- tryCatch(.parse_mmcif_atom_site(tbl$X1), error = function(e) e)
    if (inherits(parsed, "error")) {
      return(list(key = rec$uniprot, error = paste0("Parse error: ", parsed$message), df = NULL, meta = rec))
    }

    # add metadata columns for convenience
    parsed$uniprot_id   <- rec$uniprot
    parsed$fragment     <- rec$fragment
    parsed$model_id     <- rec$model_id
    parsed$model_version<- rec$model_version
    parsed$file_url     <- rec$file_url

    list(key = rec$uniprot, error = NULL, df = parsed, meta = rec)
  })

  # report per-file errors
  got_err <- Filter(function(x) !is.null(x$error), fetched)
  if (length(got_err)) {
    err_tbl <- tibble::tibble(
      id    = vapply(got_err, function(x) x$meta$uniprot, ""),
      frag  = vapply(got_err, function(x) x$meta$fragment %||% "F?", ""),
      error = vapply(got_err, `[[`, "", "error")
    ) |> dplyr::distinct()
    message("The following downloads/parses failed:")
    message(paste0(utils::capture.output(err_tbl), collapse = "\n"))
  }

  # assemble output
  ok <- Filter(function(x) !is.null(x$df), fetched)
  if (!length(ok)) {
    message("No valid information could be retrieved!")
    return(invisible(NULL))
  }

  # by default: list per UniProt (binds fragments inside each ID)
  out_list <- lapply(split(ok, vapply(ok, `[[`, "", "key")), function(items) {
    dfs <- lapply(items, `[[`, "df")
    dplyr::bind_rows(dfs)
  })

  if (!return_data_frame) {
    # keep original shape: named list by UniProt ID
    ids <- vapply(split(ok, vapply(ok, `[[`, "", "key")), function(x) x[[1]]$key, "")
    names(out_list) <- ids
    return(out_list)
  } else {
    # single data frame of all rows
    return(dplyr::bind_rows(out_list))
  }
}

