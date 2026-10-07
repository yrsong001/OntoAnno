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
#' @return The same list with (1) `final_summary` gaining the columns `compound_runs`,
#'   `n_runs`, `compound_fraction`, `compound_labels` and `flag`, and (2) a new element
#'   `mixed_summary`: one row per cluster that had at least one compound run, with the
#'   same columns, ordered by compound fraction. The console message reports the
#'   run counts so that a compound label seen in only a few runs is not over-read.
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
  # Separate table: only clusters with at least one compound run, most frequent first
  ms <- info[info$compound_runs > 0, c("cluster", "compound_runs", "n_runs", "compound_fraction", "compound_labels", "flag")]
  ms <- ms[order(-ms$compound_fraction, -ms$compound_runs), ]
  rownames(ms) <- NULL
  annotation_summary$mixed_summary <- ms
  if (nrow(ms)) {
    message("[OntoAnno] compound (multi-cell-type) predictions: ",
            paste(sprintf("cluster %s in %d/%d runs (%s)", ms$cluster, ms$compound_runs, ms$n_runs, ms$compound_labels), collapse = "; "),
            ". A compound label means the markers looked like two cell types at once in that run, not a choice between two names. ",
            "Clusters compound in >= ", round(100 * min_fraction), "% of runs are flagged 'mixed-signal': check for doublets, ",
            "under-clustering or ambient contamination before accepting the top label; a compound label in only a few runs may be noise.")
  }
  annotation_summary
}
