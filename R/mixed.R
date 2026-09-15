#' Flag clusters whose LLM predictions were compound (mixed) labels
#'
#' When the LLM answers a cluster with a compound label such as
#' "cardiomyocyte|endothelial cell", OntoAnno splits the weight between the two
#' terms, so the final summary shows the first term as "selected" with the second
#' as a ~50% "alternative". That is not a choice between two candidates: it is the
#' model saying the marker list looks like two cell types at once. This function
#' counts, per cluster, how many runs returned a compound label and attaches a flag
#' so the user checks such clusters for doublets, under-clustering or ambient
#' contamination (higher resolution, marker co-expression, doublet scores) before
#' accepting the top label.
#'
#' @param annotation_summary One element of an `ontoanno()` result (a list with
#'   `combined_results` and `final_summary`).
#' @param min_fraction Fraction of runs with a compound label at or above which a
#'   cluster is flagged "mixed-signal" (default 0.5); any compound run below that
#'   gives "partly-mixed".
#' @return The same list with `final_summary` gaining the columns `compound_runs`,
#'   `n_runs`, `compound_fraction`, `compound_labels` and `flag`.
#' @export
flag_mixed_predictions <- function(annotation_summary, min_fraction = 0.5) {
  cr <- as.data.frame(annotation_summary$combined_results)
  fs <- as.data.frame(annotation_summary$final_summary)
  cols <- setdiff(colnames(cr), "run")
  info <- do.call(rbind, lapply(cols, function(k) {
    v <- as.character(cr[[k]]); v <- v[!is.na(v) & nzchar(v)]
    mx <- grepl("|", v, fixed = TRUE)
    data.frame(cluster = as.character(k), compound_runs = sum(mx), n_runs = length(v),
               compound_fraction = if (length(v)) round(sum(mx) / length(v), 2) else NA_real_,
               compound_labels = paste(unique(tolower(gsub("\\s*\\|\\s*", "|", v[mx]))), collapse = " / "),
               stringsAsFactors = FALSE)
  }))
  info$flag <- ifelse(is.na(info$compound_fraction), "",
                      ifelse(info$compound_fraction >= min_fraction, "mixed-signal",
                             ifelse(info$compound_fraction > 0, "partly-mixed", "")))
  fs$cluster <- as.character(fs$cluster)
  for (col in c("compound_runs", "n_runs", "compound_fraction", "compound_labels", "flag")) fs[[col]] <- NULL
  ord <- fs$cluster
  fs <- merge(fs, info, by = "cluster", all.x = TRUE, sort = FALSE)
  fs <- fs[match(ord, fs$cluster), ]
  rownames(fs) <- NULL
  annotation_summary$final_summary <- fs
  flagged <- fs$cluster[!is.na(fs$flag) & fs$flag == "mixed-signal"]
  if (length(flagged)) {
    message("[OntoAnno] ", length(flagged), " cluster(s) flagged 'mixed-signal' (compound label in >= ",
            round(100 * min_fraction), "% of runs): ", paste(flagged, collapse = ", "),
            ". A compound label means the markers look like two cell types at once, not that the model is choosing between two names. ",
            "Check these clusters for doublets, under-clustering or ambient contamination before accepting the top label.")
  }
  annotation_summary
}
