#' Score GPT Annotation Results Across Multiple Resolutions
#'
#' Calculates composite scores for annotation results at different resolutions,
#' considering ontology path length (if available) and annotation confidence metrics.
#'
#' @param annotation_result_list A named list of annotation results (e.g., output from `run_annotation_all_resolutions()`).
#' @param output_csv Optional. File path to write the summary table as CSV.
#'
#' @return A data.frame summarizing each resolution's metrics and a composite score (1 is ideal: short ontology distance, high percentage).
#'
#' @details
#' The composite score is calculated by combining normalized and scaled metrics:
#'
#' **With ontology distance** (when avg_distance is available):
#' - `avg_path_length`: Mean of ontology distances across all clusters
#'   (lower is better) → normalized to 0-1 across resolutions, inverted.
#'   NA distances are treated as 0
#' - `avg_max_percentage`: Mean of max vote percentages
#'   (higher is better) → scaled directly to 0-1 by dividing by 100
#' - `min_max_percentage`: Minimum max vote percentage
#'   (higher is better) → scaled directly to 0-1 by dividing by 100
#' - Composite = (norm_avg_dist + scaled_avg + scaled_min) / 3
#'
#' **Without ontology distance** (when avg_distance is not available):
#' - `avg_max_percentage`: Mean of max vote percentages
#'   (higher is better) → scaled directly to 0-1 by dividing by 100
#' - `min_max_percentage`: Minimum max vote percentage
#'   (higher is better) → scaled directly to 0-1 by dividing by 100
#' - Composite = (scaled_avg + scaled_min) / 2
#'
#' **Normalization vs Scaling:**
#' - Path length: Normalized across resolutions (min→0, max→1)
#'   to make distances comparable
#' - Percentages: Scaled directly (percentage/100) to preserve
#'   absolute quality differences between resolutions
#' @examples
#' # result_list <- list(res_01 = ..., res_02 = ...)  # Your annotation results
#' # score_annotation_resolutions(result_list)
#' @keywords internal
#' @noRd
.is_unknown_annotation_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- sub("^\\s*(-|[0-9]+\\.)\\s*", "", x, perl = TRUE)
  x <- sub("^[-\\s]+", "", x, perl = TRUE)
  x <- trimws(x)
  x[!grepl("[a-z]", x)] <- ""
  x %in% c("", "unknown")
}

.clusters_all_unknown <- function(result_summary) {
  combined_results <- tryCatch(result_summary$combined_results, error = function(e) NULL)
  if (is.null(combined_results) || !is.data.frame(combined_results) || nrow(combined_results) == 0) {
    return(character(0))
  }

  cluster_cols <- setdiff(colnames(combined_results), "run")
  if (length(cluster_cols) == 0) {
    return(character(0))
  }

  cluster_cols[vapply(cluster_cols, function(cluster_col) {
    cluster_values <- combined_results[[cluster_col]]
    length(cluster_values) > 0 && all(.is_unknown_annotation_value(cluster_values))
  }, logical(1))]
}

.filter_final_summary_for_scoring <- function(result_summary, final_summary = result_summary$final_summary) {
  if (is.null(final_summary) || !is.data.frame(final_summary) || nrow(final_summary) == 0) {
    return(final_summary)
  }
  if (!"cluster" %in% colnames(final_summary)) {
    return(final_summary)
  }

  excluded_clusters <- .clusters_all_unknown(result_summary)
  if (length(excluded_clusters) == 0) {
    return(final_summary)
  }

  final_summary[!(as.character(final_summary$cluster) %in% as.character(excluded_clusters)), , drop = FALSE]
}

.compute_resolution_score_metrics <- function(result_summary, include_distance = FALSE) {
  final_summary <- .filter_final_summary_for_scoring(result_summary)

  if (is.null(final_summary) || !is.data.frame(final_summary) || nrow(final_summary) == 0) {
    return(list(
      avg_path_length = NA_real_,
      avg_max_percentage = NA_real_,
      min_max_percentage = NA_real_
    ))
  }

  max_percentage <- final_summary$max_percentage
  avg_max_percentage <- mean(max_percentage, na.rm = TRUE)
  if (is.nan(avg_max_percentage)) {
    avg_max_percentage <- NA_real_
  }

  min_max_percentage <- suppressWarnings(min(max_percentage, na.rm = TRUE))
  if (is.infinite(min_max_percentage)) {
    min_max_percentage <- NA_real_
  }

  avg_path_length <- NA_real_
  if (include_distance && "avg_distance" %in% colnames(final_summary)) {
    distances <- final_summary$avg_distance
    distances[is.na(distances)] <- 0
    avg_path_length <- mean(distances)
    if (is.nan(avg_path_length)) {
      avg_path_length <- NA_real_
    }
  }

  list(
    avg_path_length = avg_path_length,
    avg_max_percentage = avg_max_percentage,
    min_max_percentage = min_max_percentage
  )
}

#' Score GPT Annotation Results Across Multiple Resolutions
#'
#' Calculates composite scores for annotation results at different resolutions,
#' considering ontology path length (if available) and annotation confidence metrics.
#'
#' @param annotation_result_list A named list of annotation results (e.g., output from `run_annotation_all_resolutions()`).
#' @param output_csv Optional. File path to write the summary table as CSV.
#'
#' @return A data.frame summarizing each resolution's metrics and a composite score.
#' @export
score_annotation_resolutions <- function(annotation_result_list, output_csv = NULL) {
  if (!is.list(annotation_result_list)) stop("Input must be a named list of results.")
  norm_vec <- function(x) {
    out <- rep(NA_real_, length(x))
    valid <- !is.na(x)
    if (!any(valid)) return(out)

    rng <- diff(range(x[valid], na.rm = TRUE))
    out[valid] <- if (rng == 0) {
      rep(1, sum(valid))
    } else {
      (x[valid] - min(x[valid], na.rm = TRUE)) / rng
    }
    out
  }

  # Check if avg_distance exists in any result
  has_distance <- any(sapply(annotation_result_list, function(res) {
    "avg_distance" %in% colnames(res$final_summary)
  }))

  metrics <- lapply(annotation_result_list, function(res) {
    .compute_resolution_score_metrics(res, include_distance = has_distance)
  })

  avg_path_length <- vapply(metrics, function(metric) metric$avg_path_length, numeric(1))
  avg_max_perc <- vapply(metrics, function(metric) metric$avg_max_percentage, numeric(1))
  min_max_perc <- vapply(metrics, function(metric) metric$min_max_percentage, numeric(1))

  # Composite score: 1 is ideal (shortest ontology distance, max % = 100)
  # Only normalize path length; percentages scaled directly to 0-1
  if (has_distance) {
    norm_avg_dist <- 1 - norm_vec(avg_path_length)  # smaller is better
    scaled_avg <- avg_max_perc / 100  # scale percentage to 0-1
    scaled_min <- min_max_perc / 100  # scale percentage to 0-1
    composite_score <- ifelse(
      is.na(norm_avg_dist) | is.na(scaled_avg) | is.na(scaled_min),
      NA_real_,
      (norm_avg_dist + scaled_avg + scaled_min) / 3
    )
  } else {
    scaled_avg <- avg_max_perc / 100  # scale percentage to 0-1
    scaled_min <- min_max_perc / 100  # scale percentage to 0-1
    composite_score <- ifelse(
      is.na(scaled_avg) | is.na(scaled_min),
      NA_real_,
      (scaled_avg + scaled_min) / 2
    )
  }

  n_mixed <- vapply(annotation_result_list, function(res) {
    fs <- res$final_summary
    if (!"flag" %in% colnames(fs)) fs <- flag_mixed_predictions(res)$final_summary
    sum(!is.na(fs$flag) & fs$flag == "mixed-signal")
  }, numeric(1))
  summary_table <- data.frame(
    resolution          = names(avg_path_length),
    avg_path_length     = avg_path_length,
    avg_max_percentage  = avg_max_perc,
    min_max_percentage  = min_max_perc,
    composite_score     = composite_score,
    n_mixed_clusters    = n_mixed,
    row.names           = NULL
  ) |>
    dplyr::arrange(dplyr::desc(composite_score))
  if (any(n_mixed > 0)) message("[OntoAnno] mixed-signal clusters per resolution: ",
    paste(names(n_mixed), n_mixed, sep = " = ", collapse = ", "),
    ". These clusters carried compound labels and should be checked (doublets / under-clustering) before the resolution is accepted.")
  if (!is.null(output_csv)) {
    utils::write.csv(summary_table, output_csv, row.names = FALSE)
    message("Score summary written to: ", output_csv)
  }
  return(summary_table)
}

#' Extract Synonyms from OBO-Style Strings
#'
#' Helper function to extract synonyms from CL ontology synonym strings.
#' @param x Character vector of OBO-style synonym fields.
#' @return Character vector of synonyms.
#' @examples
#' extract_synonyms("fibroblast")
#' @export
extract_synonyms <- function(x) {
  if (length(x) == 0) return(character(0))
  syns <- stringr::str_match(x, '^"([^"]+)"')[,2]
  syns[!is.na(syns)]
}

#' Build Name/Synonym-to-CL-ID Map from Ontology Object
#'
#' Constructs a lookup table mapping CL names and synonyms to CL IDs.
#'
#'
#' @param cl An ontologyIndex CL ontology object.
#' @param verbose Logical; print mapping stats.
#' @return Data.frame with columns: key, clid, cl_label.
#' @examples
#' # cl <- ontologyIndex::get_ontology(
#' #   "http://purl.obolibrary.org/obo/cl.obo",
#' #   extract_tags = "everything"
#' # )
#' # cl_term_map <- build_cl_term_map(cl)
#' @export
build_cl_term_map <- function(cl, verbose = TRUE) {
  all_ids <- names(cl$name)
  all_labels <- cl$name
  all_labels_lower <- tolower(all_labels)
  df <- data.frame(key = all_labels_lower, clid = all_ids, cl_label = all_labels, stringsAsFactors = FALSE)
  for (clid in names(cl$synonym)) {
    syn_vec <- extract_synonyms(cl$synonym[[clid]])
    if (length(syn_vec)) {
      syn_vec <- tolower(syn_vec)
      syn_vec <- setdiff(syn_vec, df$key)
      if (length(syn_vec)) {
        df <- rbind(
          df,
          data.frame(key = syn_vec, clid = clid, cl_label = cl$name[[clid]], stringsAsFactors = FALSE)
        )
      }
    }
  }
  if (verbose) message("Mapping includes ", nrow(df), " name/synonym entries.")
  df <- unique(df)
  rownames(df) <- NULL
  return(df)
}

#' Map Cell Type Names to CL IDs
#'
#' Maps character vector of cell type names to CL ontology IDs using a prebuilt mapping table.
#'
#'
#' @param terms Character vector of cell type names.
#' @param cl_term_map Data.frame; defaults to package built-in map, or returned by `build_cl_term_map()`.
#' @param verbose Print mapping results.
#' @return Data.frame with columns: term, clid, cl_label.
#' @examples
#' # cell_names <- c("fibroblast", "endothelial cell")
#' # map_celltypes_to_cl(cell_names, OntoAnno::cl_term_map)
#' @export
map_celltypes_to_cl <- function(terms, cl_term_map = OntoAnno::cl_term_map, verbose = TRUE) {
  terms <- as.character(terms)
  keys <- tolower(terms)
  match_idx <- match(keys, cl_term_map$key)
  clid <- cl_term_map$clid[match_idx]
  cl_label <- cl_term_map$cl_label[match_idx]
  if (verbose) {
    for (i in seq_along(terms)) {
      if (!is.na(clid[i])) {
        message("Mapped: '", terms[i], "' to ", clid[i], " (", cl_label[i], ")")
      } else {
        message("No match: '", terms[i], "'")
      }
    }
  }
  return(data.frame(term = terms, clid = clid, cl_label = cl_label, stringsAsFactors = FALSE))
}

#' Calculate Mean Shortest Path Distance in CL Ontology
#'
#' For each row of a CL ID dataframe, calculates the ontology graph shortest path distance between pairs.
#'
#' @param clid_df Data.frame with columns clid1 and clid2 (CL IDs).
#' @param graph igraph object representing CL ontology DAG.
#' @param verbose Print progress and stats.
#' @return List: mean_distance, inverse_mean, dist_vector.
#' @importFrom igraph V shortest_paths
#' @importFrom utils flush.console
#' @examples
#' # res <- calculate_mean_ontology_distance(clid_df, graph)
#' @export
calculate_mean_ontology_distance <- function(clid_df, graph, verbose = TRUE) {
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  n <- nrow(clid_df)
  dist_vector <- rep(NA_real_, n)
  last_pct <- -1
  for (i in seq_len(n)) {
    id1 <- clid_df$clid1[i]
    id2 <- clid_df$clid2[i]
    dist <- tryCatch({
      if (is.na(id1) || is.na(id2)) NA_real_
      else if (id1 %in% igraph::V(graph)$name && id2 %in% igraph::V(graph)$name) {
        sp <- igraph::shortest_paths(graph, from = id1, to = id2, mode = "all")$vpath[[1]]
        if (length(sp) > 0) length(sp) - 1 else NA_real_
      } else NA_real_
    }, error = function(e) NA_real_)
    dist_vector[i] <- dist
    if (verbose && n > 1) {
      pct <- floor(100 * i / n)
      if (pct != last_pct) {
        cat(sprintf("\rProgress: %d%%", pct))
        utils::flush.console()
        last_pct <- pct
      }
    }
  }
  if (verbose && n > 1) cat("\n")
  mean_distance <- mean(dist_vector, na.rm = TRUE)
  inverse_mean <- 1 / mean_distance
  if (verbose) {
    message("Mean ontology distance: ", mean_distance)
    message("Inverse (1 / mean): ", inverse_mean)
  }
  return(list(mean_distance = mean_distance, inverse_mean = inverse_mean, dist_vector = dist_vector))
}

#' Score Annotation Distance Based on Cell Ontology
#'
#' Computes the mean shortest-path distance between annotations in two columns using CL ontology.
#' (Distance-based complement to score_annotation_agreement_ontology.)
#'
#' @param manual_clid Character. Reference/manual CL ID.
#' @param predicted_clid Character. Predicted CL ID.
#' @param cl_ontology ontologyIndex CL ontology object.
#' @param graph Ontology graph built from CL object.
#' @return List with the mapped terms, relationship label, score, and ontology distance.
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' graph <- build_ontology_graph(cl)
#' rel <- check_cl_relationship("CL:0000066", "CL:0000312", cl, graph)
#' print(rel)  # Should show child relationship
#' }
#' @export
check_cl_relationship <- function(manual_clid, predicted_clid, cl_ontology, graph) {
  scoring_weights <- c("exact" = 1.0, "parent" = 0.75, "child" = 0.75, "sibling" = 0.5, "no_match" = 0.0)

  # Create single-row data frame for classification
  test_pair <- data.frame(
    manual_clid = manual_clid,
    predicted_clid = predicted_clid,
    stringsAsFactors = FALSE
  )

  # Add names
  test_pair <- .add_cl_names_to_combinations(test_pair, cl_ontology)

  # Classify
  result <- .classify_cl_matches(test_pair, cl_ontology, graph, scoring_weights)

  return(list(
    manual_term = paste(manual_clid, "-", result$manual_cl_name[1]),
    predicted_term = paste(predicted_clid, "-", result$predicted_cl_name[1]),
    relationship = result$match_type[1],
    score = result$agreement_score[1],
    distance = result$ontology_distance[1]
  ))
}

# INTERNAL HELPER FUNCTIONS

#' Add Cell Type Names to Unique Combinations
#'
#' @param unique_combinations Data.frame with unique manual_clid/predicted_clid pairs
#' @param cl_ontology ontologyIndex CL object
#' @return Data.frame with cell type names added for each unique combination
#' @keywords internal
.add_cl_names_to_combinations <- function(unique_combinations, cl_ontology) {
  # Add cell type names for better readability
  unique_combinations$manual_cl_name <- character(nrow(unique_combinations))
  unique_combinations$predicted_cl_name <- character(nrow(unique_combinations))

  for (i in 1:nrow(unique_combinations)) {
    manual_id <- unique_combinations$manual_clid[i]
    predicted_id <- unique_combinations$predicted_clid[i]

    # Get manual cell type name
    if (!is.na(manual_id) && manual_id %in% names(cl_ontology$name)) {
      unique_combinations$manual_cl_name[i] <- cl_ontology$name[[manual_id]]
    } else {
      unique_combinations$manual_cl_name[i] <- NA_character_
    }

    # Get predicted cell type name
    if (!is.na(predicted_id) && predicted_id %in% names(cl_ontology$name)) {
      unique_combinations$predicted_cl_name[i] <- cl_ontology$name[[predicted_id]]
    } else {
      unique_combinations$predicted_cl_name[i] <- NA_character_
    }
  }

  return(unique_combinations)
}

#' Get Direct Parents for CL IDs (Internal)
#'
#' @param clids Character vector of CL IDs
#' @param cl_ontology ontologyIndex CL object
#' @return Character vector of unique parent CL IDs
#' @keywords internal
.get_direct_parents <- function(clids, cl_ontology) {
  if (length(clids) == 0 || all(is.na(clids))) return(character(0))

  all_parents <- character(0)
  for (clid in clids) {
    if (!is.na(clid) && clid %in% names(cl_ontology$parents)) {
      all_parents <- c(all_parents, cl_ontology$parents[[clid]])
    }
  }
  return(unique(all_parents))
}

#' Get Direct Children for CL IDs (Internal)
#'
#' @param clids Character vector of CL IDs
#' @param cl_ontology ontologyIndex CL object
#' @return Character vector of unique child CL IDs
#' @keywords internal
.get_direct_children <- function(clids, cl_ontology) {
  if (length(clids) == 0 || all(is.na(clids))) return(character(0))

  all_children <- character(0)
  valid_clids <- clids[!is.na(clids)]

  for (clid in valid_clids) {
    children <- names(cl_ontology$parents)[
      vapply(cl_ontology$parents, function(parents) clid %in% parents, logical(1))
    ]
    all_children <- c(all_children, children)
  }
  return(unique(all_children))
}

#' Calculate Node Depths in Ontology Graph (Internal)
#'
#' @param clids Character vector of CL IDs
#' @param graph igraph object
#' @return Named numeric vector of depths from root
#' @keywords internal
.get_node_depths <- function(clids, graph) {
  depths <- stats::setNames(rep(NA_real_, length(clids)), clids)

  # Find root nodes (no incoming edges)
  in_degrees <- igraph::degree(graph, mode = "in")
  roots <- names(in_degrees[in_degrees == 0])

  if (length(roots) == 0) {
    roots <- names(sort(igraph::degree(graph), decreasing = TRUE))[1]
  }

  root <- roots[1]  # Use first root for consistency

  # Calculate distances from root
  valid_clids <- clids[clids %in% igraph::V(graph)$name & !is.na(clids)]

  if (length(valid_clids) > 0 && root %in% igraph::V(graph)$name) {
    distances <- igraph::distances(graph, v = root, to = valid_clids, mode = "out")
    depths[valid_clids] <- as.numeric(distances[1, ])
    depths[is.infinite(depths)] <- NA_real_
  }

  return(depths)
}

#' Calculate Pairwise Distances Between CL IDs (Internal)
#'
#' @param pairs Data.frame with 'from' and 'to' columns containing CL IDs
#' @param graph igraph object
#' @return Numeric vector of shortest path distances
#' @keywords internal
.get_pairwise_distances <- function(pairs, graph) {
  n_pairs <- nrow(pairs)
  distances <- rep(NA_real_, n_pairs)

  # Filter to valid pairs
  valid_mask <- !is.na(pairs$from) & !is.na(pairs$to) &
    pairs$from %in% igraph::V(graph)$name &
    pairs$to %in% igraph::V(graph)$name

  if (sum(valid_mask) > 0) {
    valid_pairs <- pairs[valid_mask, ]
    unique_from <- unique(valid_pairs$from)

    for (from_node in unique_from) {
      mask <- valid_pairs$from == from_node
      to_nodes <- valid_pairs$to[mask]

      if (length(to_nodes) > 0) {
        dist_matrix <- igraph::distances(graph, v = from_node, to = to_nodes, mode = "all")
        original_indices <- which(valid_mask)[mask]
        distances[original_indices] <- as.numeric(dist_matrix[1, ])
      }
    }
  }

  distances[is.infinite(distances)] <- NA_real_
  return(distances)
}

#' Classify Unique CL ID Combinations (Internal)
#'
#' @param unique_pairs Data.frame with manual_clid and predicted_clid columns
#' @param cl_ontology ontologyIndex CL object
#' @param graph igraph object
#' @param scoring_weights Named vector of match type weights
#' @return Data.frame with match classifications
#' @keywords internal
.classify_cl_matches <- function(unique_pairs, cl_ontology, graph, scoring_weights) {
  message("Classifying ", nrow(unique_pairs), " unique CL ID combinations...")

  result_df <- unique_pairs
  result_df$match_type <- "no_match"
  result_df$agreement_score <- scoring_weights["no_match"]
  result_df$ontology_distance <- NA_real_

  # Handle exact matches
  exact_mask <- !is.na(result_df$manual_clid) &
    !is.na(result_df$predicted_clid) &
    result_df$manual_clid == result_df$predicted_clid

  result_df$match_type[exact_mask] <- "exact"
  result_df$agreement_score[exact_mask] <- scoring_weights["exact"]
  result_df$ontology_distance[exact_mask] <- 0

  # Process non-exact matches
  non_exact_mask <- !exact_mask & !is.na(result_df$manual_clid) & !is.na(result_df$predicted_clid)

  if (sum(non_exact_mask) > 0) {
    non_exact_data <- result_df[non_exact_mask, ]

    # Calculate distances
    distance_pairs <- data.frame(
      from = non_exact_data$manual_clid,
      to = non_exact_data$predicted_clid
    )

    distances <- .get_pairwise_distances(distance_pairs, graph)
    result_df$ontology_distance[non_exact_mask] <- distances

    # Process distance = 1 (parent/child relationships)
    result_df <- .process_distance_one_relationships(result_df, cl_ontology, scoring_weights)

    # Process distance = 2 (potential siblings)
    result_df <- .process_distance_two_relationships(result_df, cl_ontology, graph, scoring_weights)
  }

  return(result_df)
}

#' Process Distance-1 Relationships (Parent/Child) (Internal)
#'
#' @param result_df Data.frame with classification results
#' @param cl_ontology ontologyIndex CL object
#' @param scoring_weights Named vector of scoring weights
#' @return Updated result_df
#' @keywords internal
.process_distance_one_relationships <- function(result_df, cl_ontology, scoring_weights) {
  dist_1_mask <- !is.na(result_df$ontology_distance) & result_df$ontology_distance == 1 &
    result_df$match_type == "no_match"

  if (sum(dist_1_mask) > 0) {
    dist_1_indices <- which(dist_1_mask)

    for (idx in dist_1_indices) {
      manual_id <- result_df$manual_clid[idx]
      predicted_id <- result_df$predicted_clid[idx]

      # Check parent relationship
      manual_parents <- .get_direct_parents(manual_id, cl_ontology)
      if (predicted_id %in% manual_parents) {
        result_df$match_type[idx] <- "parent"
        result_df$agreement_score[idx] <- scoring_weights["parent"]
        next
      }

      # Check child relationship
      manual_children <- .get_direct_children(manual_id, cl_ontology)
      if (predicted_id %in% manual_children) {
        result_df$match_type[idx] <- "child"
        result_df$agreement_score[idx] <- scoring_weights["child"]
      }
    }
  }

  return(result_df)
}

#' Process Distance-2 Relationships (Siblings) (Internal)
#'
#' @param result_df Data.frame with classification results
#' @param cl_ontology ontologyIndex CL object
#' @param graph igraph object
#' @param scoring_weights Named vector of scoring weights
#' @return Updated result_df
#' @keywords internal
.process_distance_two_relationships <- function(result_df, cl_ontology, graph, scoring_weights) {
  dist_2_mask <- !is.na(result_df$ontology_distance) & result_df$ontology_distance == 2 &
    result_df$match_type == "no_match"

  if (sum(dist_2_mask) > 0) {
    dist_2_indices <- which(dist_2_mask)

    # Calculate depths for sibling verification
    all_relevant_terms <- unique(c(
      result_df$manual_clid[dist_2_indices],
      result_df$predicted_clid[dist_2_indices]
    ))
    term_depths <- .get_node_depths(all_relevant_terms, graph)

    for (idx in dist_2_indices) {
      manual_id <- result_df$manual_clid[idx]
      predicted_id <- result_df$predicted_clid[idx]

      # Check same depth requirement for siblings
      manual_depth <- term_depths[manual_id]
      predicted_depth <- term_depths[predicted_id]

      if (!is.na(manual_depth) && !is.na(predicted_depth) && manual_depth == predicted_depth) {
        # Check if they share a parent
        manual_parents <- .get_direct_parents(manual_id, cl_ontology)
        if (length(manual_parents) > 0) {
          parent_children <- .get_direct_children(manual_parents, cl_ontology)
          if (predicted_id %in% parent_children) {
            result_df$match_type[idx] <- "sibling"
            result_df$agreement_score[idx] <- scoring_weights["sibling"]
          }
        }
      }
    }
  }

  return(result_df)
}

#' Map Classifications Back to Individual Cells (Internal)
#'
#' @param cell_results Data.frame with cell-level data
#' @param classified_combinations Data.frame with unique classifications
#' @return Updated cell_results with match classifications
#' @keywords internal
.map_classifications_to_cells <- function(cell_results, classified_combinations) {
  # Initialize result columns
  cell_results$match_type <- NA_character_
  cell_results$agreement_score <- NA_real_
  cell_results$ontology_distance <- NA_real_

  # Create lookup keys
  cell_keys <- paste(cell_results$manual_clid, cell_results$predicted_clid, sep = "|||")
  combo_keys <- paste(classified_combinations$manual_clid, classified_combinations$predicted_clid, sep = "|||")

  # Map results back
  match_indices <- match(cell_keys, combo_keys)

  valid_matches <- !is.na(match_indices)
  cell_results$match_type[valid_matches] <- classified_combinations$match_type[match_indices[valid_matches]]
  cell_results$agreement_score[valid_matches] <- classified_combinations$agreement_score[match_indices[valid_matches]]
  cell_results$ontology_distance[valid_matches] <- classified_combinations$ontology_distance[match_indices[valid_matches]]

  return(cell_results)
}

#' Calculate Summary Statistics (Internal)
#'
#' @param cell_results Data.frame with cell-level results
#' @param classified_combinations Data.frame with unique classifications
#' @param scoring_weights Named vector of scoring weights
#' @return List with summary statistics
#' @keywords internal
.calculate_summary_statistics <- function(cell_results, classified_combinations, scoring_weights) {
  valid_scores <- cell_results$agreement_score[!is.na(cell_results$agreement_score)]
  valid_types <- cell_results$match_type[!is.na(cell_results$match_type)]

  type_counts <- table(factor(valid_types, levels = c("exact", "parent", "child", "sibling", "no_match")))
  type_proportions <- prop.table(type_counts) * 100

  list(
    total_cells = nrow(cell_results),
    valid_cells = length(valid_scores),
    unique_combinations = nrow(classified_combinations),
    mean_score = round(mean(valid_scores), 4),
    median_score = round(stats::median(valid_scores), 4),

    # Counts
    n_exact = type_counts["exact"],
    n_parent = type_counts["parent"],
    n_child = type_counts["child"],
    n_sibling = type_counts["sibling"],
    n_no_match = type_counts["no_match"],

    # Percentages
    pct_exact = round(type_proportions["exact"], 2),
    pct_parent = round(type_proportions["parent"], 2),
    pct_child = round(type_proportions["child"], 2),
    pct_sibling = round(type_proportions["sibling"], 2),
    pct_no_match = round(type_proportions["no_match"], 2),

    scoring_weights = scoring_weights
  )
}

#' Print Results Summary (Internal)
#'
#' @param summary_stats List with summary statistics
#' @keywords internal
.print_results_summary <- function(summary_stats) {
  message("\n=== SCORING RESULTS ===")
  message("Total cells processed: ", summary_stats$total_cells)
  message("Valid annotations: ", summary_stats$valid_cells)
  message("Unique combinations: ", summary_stats$unique_combinations)
  message("Mean agreement score: ", summary_stats$mean_score)
  message("Median agreement score: ", summary_stats$median_score)
  message("\nMatch Type Distribution:")
  message("  Exact matches: ", summary_stats$n_exact, " (", summary_stats$pct_exact, "%)")
  message("  Parent matches: ", summary_stats$n_parent, " (", summary_stats$pct_parent, "%)")
  message("  Child matches: ", summary_stats$n_child, " (", summary_stats$pct_child, "%)")
  message("  Sibling matches: ", summary_stats$n_sibling, " (", summary_stats$pct_sibling, "%)")
  message("  No matches: ", summary_stats$n_no_match, " (", summary_stats$pct_no_match, "%)")
}

#' Score Annotation Distance Based on Cell Ontology
#'
#' Computes the mean shortest-path distance between annotations in two columns using CL ontology.
#' (Distance-based complement to score_annotation_agreement_ontology.)
#'
#' @param seurat_obj Seurat object.
#' @param manual_col Character. Metadata column name for reference/manual annotation.
#' @param predicted_col Character. Metadata column for predicted/auto annotation.
#' @param cl_term_map Data.frame; defaults to package built-in map, or returned by `build_cl_term_map()`.
#' @param graph Ontology graph built from CL object.
#' @param verbose Print mapping/progress.
#' @return List with mean distance and mapping dataframe result.
#' @examples
#' \dontrun{
#' # Build ontology and graph (one-time setup)
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo",
#'  extract_tags = "everything")
#' graph <- build_ontology_graph(cl)
#' cl_term_map <- build_cl_term_map(cl) # if need the up-to-date mapping
#' # Score annotation distance
#' dist_res <- score_annotation_distance_ontology(
#'   seurat_obj, manual_col = "celltype_manual", predicted_col = "celltype_parent",
#'   cl_term_map = cl_term_map, # leave it blank will use the built-in mapping
#'   graph = ontology_graph
#' )
#' print(dist_res$mean_distance)
#' }
#'
#' @export
#' @importFrom dplyr mutate recode transmute
#' @importFrom magrittr %>%
score_annotation_distance_ontology <- function(
    seurat_obj,
    manual_col = "manual_celltype",
    predicted_col = "predicted_celltype",
    cl_term_map = OntoAnno::cl_term_map,
    graph,
    verbose = TRUE
) {
  meta <- seurat_obj@meta.data
  meta <- meta[!is.na(meta[[manual_col]]) & !is.na(meta[[predicted_col]]), ]

  # Apply clean_and_match_annotation following score_old.R convention
  cl_map_col1 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[manual_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col1$label0 <- names(table(meta[[manual_col]]))
  cl_map_col2 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[predicted_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col2$label0 <- names(table(meta[[predicted_col]]))
  cl_map <- rbind(cl_map_col1, cl_map_col2)

  map_name <- stats::setNames(cl_map$cl_label, cl_map$label0)
  map_id <- stats::setNames(cl_map$clid, cl_map$label0)
  meta <- meta %>%
    mutate(col1_cl = recode(meta[[manual_col]], !!!map_name)) %>%
    mutate(col1_clid = recode(meta[[manual_col]], !!!map_id)) %>%
    mutate(col2_cl = recode(meta[[predicted_col]], !!!map_name)) %>%
    mutate(col2_clid = recode(meta[[predicted_col]], !!!map_id))

  clid_df <- meta %>%
    mutate(cell = row.names(.)) %>%
    transmute(cell, manual_term = meta[[manual_col]], predicted_term = meta[[predicted_col]],
              clid1 = col1_clid, clid2 = col2_clid,
              manual_label = col1_cl, predicted_label = col2_cl)

  dist_res <- calculate_mean_ontology_distance(clid_df, graph, verbose = verbose)
  dist_res$clid_df <- clid_df
  return(dist_res)
}

#' Score Agreement Based on Cell Ontology Ancestry
#'
#' Scoring strategy: Fully match includes same CL ID or manual is ancestor of predicted.
#' Partially match when manual is child of predicted.
#'
#' **NOTE:** This function uses revised scoring (1/0.5/0). For detailed 5-category
#' scoring, use \code{score_annotation_agreement_ontology_detailed()}.
#'
#' @param seurat_obj A Seurat object with both annotation columns.
#' @param manual_col Character. Metadata column for manual annotation.
#' @param predicted_col Character. Metadata column for predicted annotation.
#' @param cl_term_map Data.frame. Term name/synonym to CL ID; defaults to package built-in map, or returned by `build_cl_term_map()`.
#' @param ancestor_type_map Named list: CL ID to character vector of ancestors (including self).
#'        \strong{Run \code{ancestor_type_map <- build_ancestor_type_map(cl)}} to create the ancestor map before using this function.
#' @param output_csv Optional. Save the detailed score of each cells. Path for CSV export. if NULL, no export.
#' @return List with per-cell scores and summary.
#' @details
#' \strong{NOTE:} You need to first run \code{ancestor_type_map <- build_ancestor_type_map(cl)}
#' to create the ancestor mapping for your ontology.
#'
#' This function provides revised scoring approach:
#' \itemize{
#'   \item Score 1.0: Exact match (same CL ID or manual is ancestor of predicted)
#'   \item Score 0.5: Partial match (manual is child of predicted)
#'   \item Score 0.0: No match (unrelated terms)
#' }
#'
#' For more detailed analysis with 5 match types (exact/parent/child/sibling/no_match),
#' use \code{score_annotation_agreement_ontology_detailed()}.
#'
#' @export
#' @importFrom dplyr mutate recode rename filter arrange desc left_join
score_annotation_agreement_ontology <- function(
    seurat_obj,
    manual_col = "manual_celltype",
    predicted_col = "predicted_celltype",
    cl_term_map = OntoAnno::cl_term_map,
    ancestor_type_map,
    output_csv = NULL
) {
  meta <- seurat_obj@meta.data
  meta <- meta[!is.na(meta[[manual_col]]) & !is.na(meta[[predicted_col]]), ]
  total_cells <- nrow(meta)  # Store total cell count for percentage calculation

  # Apply clean_and_match_annotation following score_old.R convention
  cl_map_col1 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[manual_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col1$label0 <- names(table(meta[[manual_col]]))
  cl_map_col2 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[predicted_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col2$label0 <- names(table(meta[[predicted_col]]))
  cl_map <- rbind(cl_map_col1, cl_map_col2)

  map_name <- stats::setNames(cl_map$cl_label, cl_map$label0)
  map_id <- stats::setNames(cl_map$clid, cl_map$label0)
  meta <- meta %>%
    mutate(col1_cl = recode(meta[[manual_col]], !!!map_name)) %>%
    mutate(col1_clid = recode(meta[[manual_col]], !!!map_id)) %>%
    mutate(col2_cl = recode(meta[[predicted_col]], !!!map_name)) %>%
    mutate(col2_clid = recode(meta[[predicted_col]], !!!map_id))

  # Helper: check ancestry
  is_ancestor <- function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    a %in% ancestor_type_map[[b]]
  }

  # Revised agreement score strategy:
  # 1.0: same CL ID OR manual is ancestor of predicted
  # 0.5: manual is child of predicted
  # 0.0: unrelated
  scores <- mapply(function(c1, c2) {
    if (is.na(c1) || is.na(c2)) return(NA_real_)
    if (c1 == c2) return(1)  # Exact match
    if (is_ancestor(c1, c2)) return(1)  # Manual is ancestor of predicted
    if (is_ancestor(c2, c1)) return(0.5)  # Manual is child of predicted
    return(0)  # Unrelated
  }, meta$col1_clid, meta$col2_clid)

  result_df <- data.frame(
    cell = rownames(meta),
    manual_term = meta[[manual_col]],
    predicted_term = meta[[predicted_col]],
    manual_clid = meta$col1_clid,
    predicted_clid = meta$col2_clid,
    manual_label = meta$col1_cl,
    predicted_label = meta$col2_cl,
    agreement_score = scores,
    stringsAsFactors = FALSE
  )

  summary <- list(
    mean_score = mean(scores, na.rm = TRUE),
    n_cells = sum(!is.na(scores)),
    fully_match = mean(scores == 1, na.rm = TRUE),
    partial_match = mean(scores == 0.5, na.rm = TRUE),
    mismatch = mean(scores == 0, na.rm = TRUE),
    unscored = mean(is.na(scores))
  )

  if (!is.null(output_csv)) {
    write.csv(result_df, file = output_csv, row.names = FALSE)
    message("Agreement scores written to: ", output_csv)
  }

  # Create pairwise cross-tabulation with agreement scores
  pairwise_raw <- table(result_df$manual_label, result_df$predicted_label, useNA = "ifany")
  pairwise <- as.data.frame(pairwise_raw) |>
    dplyr::rename(manual_label = Var1, predicted_label = Var2, freq = Freq) |>
    dplyr::filter(freq > 0, as.character(manual_label) != as.character(predicted_label)) |>
    dplyr::mutate(percent = freq / total_cells * 100)  # Use total cell count

  # Add agreement scores to pairwise table
  # Calculate mean agreement score for each manual-predicted pair
  score_summary <- result_df |>
    dplyr::group_by(manual_label, predicted_label) |>
    dplyr::summarize(agreement_score = mean(agreement_score, na.rm = TRUE), .groups = "drop")

  pairwise <- pairwise |>
    dplyr::left_join(score_summary, by = c("manual_label", "predicted_label")) |>
    dplyr::arrange(dplyr::desc(freq))

  return(list(scores = result_df, summary = summary, pairwise = pairwise))
}

# DETAILED SCORING FUNCTIONS

#' Score Annotation Agreement Based on Cell Ontology with Detailed Match Types
#'
#' Enhanced version of \code{score_annotation_agreement_ontology()} with detailed
#' match classification using the same ancestor relationship strategy.
#'
#' @param seurat_obj A Seurat object with both annotation columns.
#' @param manual_col Character. Metadata column name for reference/manual annotation.
#' @param predicted_col Character. Metadata column name for predicted/auto annotation.
#' @param cl_term_map Data.frame with mapping from cell type names to CL IDs.
#'   Should have columns: key, clid, cl_label. Defaults to package built-in map.
#' @param cl_ontology The CL ontology object from \code{ontologyIndex::get_ontology()}.
#' @param graph igraph object representing CL ontology DAG from \code{build_ontology_graph()}.
#' @param ancestor_type_map Named list: CL ID to character vector of ancestors (including self).
#'        \strong{Run \code{ancestor_type_map <- build_ancestor_type_map(cl)}} to create the ancestor map.
#' @param scoring_weights Named numeric vector with weights for each match type.
#'   Default: c("exact" = 1.0, "parent" = 0.8, "child" = 0.7, "sibling" = 0.5, "no_match" = 0.0)
#' @param output_csv Optional character. File path to save detailed per-cell results as CSV.
#'
#' @return List with components:
#' \describe{
#'   \item{scores}{Data.frame with per-cell detailed results}
#'   \item{summary}{List with summary statistics and match type counts}
#'   \item{pairwise}{Data.frame with unique combination classifications including freq and percent}
#' }
#'
#' @details
#' Match types using ancestor relationships:
#' \itemize{
#'   \item \strong{Exact}: Same CL ID
#'   \item \strong{Parent}: Predicted is an ancestor of manual
#'   \item \strong{Child}: Manual is an ancestor of predicted
#'   \item \strong{Sibling}: Distance <= 3 and share common ancestors but neither is ancestor of the other
#'   \item \strong{No match}: None of the above criteria are met
#' }
#'
#' @examples
#' \dontrun{
#' # Setup ontology (one-time)
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo",
#'                                   extract_tags = "everything")
#' graph <- build_ontology_graph(cl)
#' cl_term_map <- build_cl_term_map(cl)
#' ancestor_type_map <- build_ancestor_type_map(cl)
#'
#' # Score annotations with detailed match types
#' results <- score_annotation_agreement_ontology_detailed(
#'   seurat_obj,
#'   manual_col = "celltype_manual",
#'   predicted_col = "celltype_predicted",
#'   cl_term_map = cl_term_map,
#'   cl_ontology = cl,
#'   graph = graph,
#'   ancestor_type_map = ancestor_type_map,
#'   output_csv = "detailed_scores.csv"
#' )
#'
#' # View summary
#' print(results$summary)
#' }
#'
#' @export
#' @importFrom utils write.csv
#' @importFrom dplyr mutate recode transmute rename filter arrange desc left_join group_by summarize
score_annotation_agreement_ontology_detailed <- function(
    seurat_obj,
    manual_col = "manual_celltype",
    predicted_col = "predicted_celltype",
    cl_term_map = OntoAnno::cl_term_map,
    cl_ontology,
    graph,
    ancestor_type_map,
    scoring_weights = c("exact" = 1.0, "parent" = 0.8, "child" = 0.7, "sibling" = 0.5, "no_match" = 0.0),
    output_csv = NULL
) {
  # Input validation
  stopifnot(requireNamespace("igraph", quietly = TRUE))
  if (!all(names(scoring_weights) %in% c("exact", "parent", "child", "sibling", "no_match"))) {
    stop("scoring_weights must contain keys: exact, parent, child, sibling, no_match")
  }
  if (!manual_col %in% colnames(seurat_obj@meta.data)) {
    stop("manual_col '", manual_col, "' not found in seurat_obj metadata")
  }
  if (!predicted_col %in% colnames(seurat_obj@meta.data)) {
    stop("predicted_col '", predicted_col, "' not found in seurat_obj metadata")
  }

  message("=== Detailed Cell Ontology Agreement Scoring ===")

  # Extract and validate data - following score_old.R convention
  meta <- seurat_obj@meta.data
  meta <- meta[!is.na(meta[[manual_col]]) & !is.na(meta[[predicted_col]]), ]
  total_cells <- nrow(meta)  # Store total cell count for percentage calculation

  if (nrow(meta) == 0) {
    stop("No valid overlapping annotations found in specified columns")
  }

  message("Processing ", nrow(meta), " cells...")

  # Map cell types to CL IDs using clean_and_match_annotation - following score_old.R convention
  cl_map_col1 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[manual_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col1$label0 <- names(table(meta[[manual_col]]))
  cl_map_col2 <- map_celltypes_to_cl(sapply(gsub("_", " ", names(table(meta[[predicted_col]]))), clean_and_match_annotation), cl_term_map, verbose = FALSE)
  cl_map_col2$label0 <- names(table(meta[[predicted_col]]))
  cl_map <- rbind(cl_map_col1, cl_map_col2)

  map_name <- stats::setNames(cl_map$cl_label, cl_map$label0)
  map_id <- stats::setNames(cl_map$clid, cl_map$label0)
  meta <- meta %>%
    mutate(col1_cl = recode(meta[[manual_col]], !!!map_name)) %>%
    mutate(col1_clid = recode(meta[[manual_col]], !!!map_id)) %>%
    mutate(col2_cl = recode(meta[[predicted_col]], !!!map_name)) %>%
    mutate(col2_clid = recode(meta[[predicted_col]], !!!map_id))

  # Create cell-level results dataframe
  cell_results <- data.frame(
    cell_id = rownames(meta),
    manual_term = meta[[manual_col]],
    predicted_term = meta[[predicted_col]],
    manual_clid = meta$col1_clid,
    predicted_clid = meta$col2_clid,
    manual_label = meta$col1_cl,
    predicted_label = meta$col2_cl,
    stringsAsFactors = FALSE
  )

  # Find unique combinations for efficient processing
  unique_combinations <- unique(cell_results[, c("manual_clid", "predicted_clid")])
  message("Found ", nrow(unique_combinations), " unique CL ID combinations")

  # Add cell type names to unique combinations
  unique_combinations_with_names <- .add_cl_names_to_combinations(unique_combinations, cl_ontology)

  # Helper function to check if a is ancestor of b (same as in score_annotation_agreement_ontology)
  is_ancestor <- function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    return(a %in% ancestor_type_map[[b]])
  }

  # Calculate distances for all pairs first
  distance_pairs <- data.frame(
    from = unique_combinations$manual_clid,
    to = unique_combinations$predicted_clid
  )
  distances <- .get_pairwise_distances(distance_pairs, graph)

  # Classify unique combinations directly in main function
  result_df <- unique_combinations_with_names
  result_df$match_type <- "no_match"
  result_df$agreement_score <- scoring_weights["no_match"]
  result_df$ontology_distance <- distances

  # Process each unique combination
  for (i in seq_len(nrow(result_df))) {
    manual_id <- result_df$manual_clid[i]
    predicted_id <- result_df$predicted_clid[i]
    distance <- result_df$ontology_distance[i]

    # Skip if either ID is NA
    if (is.na(manual_id) || is.na(predicted_id)) {
      next
    }

    # Check exact match
    if (manual_id == predicted_id) {
      result_df$match_type[i] <- "exact"
      result_df$agreement_score[i] <- scoring_weights["exact"]
      result_df$ontology_distance[i] <- 0
      next
    }

    # Check parent relationship: predicted is ancestor of manual
    if (is_ancestor(predicted_id, manual_id)) {
      result_df$match_type[i] <- "parent"
      result_df$agreement_score[i] <- scoring_weights["parent"]
      next
    }

    # Check child relationship: manual is ancestor of predicted
    if (is_ancestor(manual_id, predicted_id)) {
      result_df$match_type[i] <- "child"
      result_df$agreement_score[i] <- scoring_weights["child"]
      next
    }

    # Check sibling relationship: distance <= 3 AND they share common ancestors
    if (!is.na(distance) && distance <= 3) {
      if (manual_id %in% names(ancestor_type_map) && predicted_id %in% names(ancestor_type_map)) {
        manual_ancestors <- ancestor_type_map[[manual_id]]
        predicted_ancestors <- ancestor_type_map[[predicted_id]]

        # Find common ancestors (excluding very general root terms)
        common_ancestors <- intersect(manual_ancestors, predicted_ancestors)

        # Check if they have meaningful common ancestors (more than just root)
        if (length(common_ancestors) > 2) {
          result_df$match_type[i] <- "sibling"
          result_df$agreement_score[i] <- scoring_weights["sibling"]
        }
      }
    }
  }

  # Add frequency and percentage information to the classified combinations
  combination_counts <- cell_results |>
    dplyr::group_by(manual_clid, predicted_clid) |>
    dplyr::summarize(freq = n(), .groups = "drop") |>
    dplyr::mutate(percent = freq / total_cells * 100)

  # Merge with classified combinations to create final pairwise table
  pairwise_table <- result_df |>
    dplyr::left_join(combination_counts, by = c("manual_clid", "predicted_clid")) |>
    dplyr::arrange(dplyr::desc(freq))

  # Map results back to cells efficiently
  cell_results <- .map_classifications_to_cells(cell_results, result_df)

  # Calculate summary statistics
  summary_stats <- .calculate_summary_statistics(cell_results, result_df, scoring_weights)

  # Print summary
  .print_results_summary(summary_stats)

  # Export results if requested
  if (!is.null(output_csv)) {
    utils::write.csv(cell_results, output_csv, row.names = FALSE)
    message("\nDetailed results saved to: ", output_csv)
  }

  return(list(
    scores = cell_results,
    summary = summary_stats,
    pairwise = pairwise_table
  ))
}

#' Classify Unique CL ID Combinations Using Ancestor Relationships (Fixed Version)
#'
#' @param unique_pairs Data.frame with manual_clid and predicted_clid columns
#' @param cl_ontology ontologyIndex CL object
#' @param graph igraph object (for distance calculation)
#' @param ancestor_type_map Named list mapping CL IDs to their ancestors
#' @param scoring_weights Named vector of match type weights
#' @return Data.frame with match classifications
#' @keywords internal
.classify_cl_matches_with_ancestors <- function(unique_pairs, cl_ontology, graph, ancestor_type_map, scoring_weights) {
  message("Classifying ", nrow(unique_pairs), " unique CL ID combinations using ancestor relationships...")

  result_df <- unique_pairs
  result_df$match_type <- "no_match"
  result_df$agreement_score <- scoring_weights["no_match"]
  result_df$ontology_distance <- NA_real_

  # Helper function to check if a is ancestor of b (same as in score_annotation_agreement_ontology)
  is_ancestor <- function(a, b) {
    if (is.na(a) || is.na(b)) return(FALSE)
    return(a %in% ancestor_type_map[[b]])
  }

  # Calculate distances for all pairs first
  distance_pairs <- data.frame(
    from = unique_pairs$manual_clid,
    to = unique_pairs$predicted_clid
  )
  distances <- .get_pairwise_distances(distance_pairs, graph)
  result_df$ontology_distance <- distances

  # Process each unique combination
  for (i in seq_len(nrow(result_df))) {
    manual_id <- result_df$manual_clid[i]
    predicted_id <- result_df$predicted_clid[i]

    # Skip if either ID is NA
    if (is.na(manual_id) || is.na(predicted_id)) {
      next
    }

    # Check exact match
    if (manual_id == predicted_id) {
      result_df$match_type[i] <- "exact"
      result_df$agreement_score[i] <- scoring_weights["exact"]
      result_df$ontology_distance[i] <- 0
      next
    }

    # Check parent relationship: predicted is ancestor of manual
    if (is_ancestor(predicted_id, manual_id)) {
      result_df$match_type[i] <- "parent"
      result_df$agreement_score[i] <- scoring_weights["parent"]
      next
    }

    # Check child relationship: manual is ancestor of predicted
    if (is_ancestor(manual_id, predicted_id)) {
      result_df$match_type[i] <- "child"
      result_df$agreement_score[i] <- scoring_weights["child"]
      next
    }

    # Check sibling relationship: they share a common ancestor but neither is ancestor of the other
    if (manual_id %in% names(ancestor_type_map) && predicted_id %in% names(ancestor_type_map)) {
      manual_ancestors <- ancestor_type_map[[manual_id]]
      predicted_ancestors <- ancestor_type_map[[predicted_id]]

      # Find common ancestors (excluding very general root terms)
      common_ancestors <- intersect(manual_ancestors, predicted_ancestors)

      # Filter out root-level terms (heuristic: keep terms with reasonable specificity)
      if (length(common_ancestors) > 2) {  # More than just basic root terms
        result_df$match_type[i] <- "sibling"
        result_df$agreement_score[i] <- scoring_weights["sibling"]
      }
    }
  }

  return(result_df)
}

# PUBLIC UTILITY FUNCTIONS FOR ONTOLOGY EXPLORATION

#' Get Direct Parents for CL IDs
#'
#' Retrieve direct parent terms for given Cell Ontology IDs.
#'
#' @param clids Character vector of CL IDs (e.g., "CL:0000066")
#' @param cl_ontology ontologyIndex CL object from \code{ontologyIndex::get_ontology()}
#' @return Character vector of unique parent CL IDs
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' parents <- get_cl_direct_parents("CL:0000066", cl)  # epithelial cell parents
#' print(parents) # CL:0000255 eukaryotic cell
#' }
#' @export
get_cl_direct_parents <- function(clids, cl_ontology) {
  .get_direct_parents(clids, cl_ontology)
}

#' Get Direct Children for CL IDs
#'
#' Retrieve direct child terms for given Cell Ontology IDs.
#'
#' @param clids Character vector of CL IDs (e.g., "CL:0000066")
#' @param cl_ontology ontologyIndex CL object from \code{ontologyIndex::get_ontology()}
#' @return Character vector of unique child CL IDs
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' children <- get_cl_direct_children("CL:0000066", cl)  # epithelial cell children
#' print(children)
#' }
#' @export
get_cl_direct_children <- function(clids, cl_ontology) {
  .get_direct_children(clids, cl_ontology)
}

#' Calculate Node Depths in CL Ontology Graph
#'
#' Calculate the depth (distance from root) for Cell Ontology terms.
#'
#' @param clids Character vector of CL IDs
#' @param graph igraph object from \code{build_ontology_graph()}
#' @return Named numeric vector of depths from root (NA for unreachable nodes)
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' graph <- build_ontology_graph(cl)
#' depths <- get_cl_node_depths(c("CL:0000066", "CL:0000057"), graph)
#' print(depths)
#' }
#' @export
get_cl_node_depths <- function(clids, graph) {
  .get_node_depths(clids, graph)
}

#' Calculate Pairwise Distances Between CL IDs
#'
#' Calculate shortest path distances between pairs of Cell Ontology terms.
#'
#' @param pairs Data.frame with 'from' and 'to' columns containing CL IDs
#' @param graph igraph object from \code{build_ontology_graph()}
#' @return Numeric vector of shortest path distances (NA for unreachable pairs)
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' graph <- build_ontology_graph(cl)
#' pairs <- data.frame(from = "CL:0000066", to = "CL:0000312")
#' distances <- get_cl_pairwise_distances(pairs, graph)
#' print(distances)
#' }
#' @export
get_cl_pairwise_distances <- function(pairs, graph) {
  .get_pairwise_distances(pairs, graph)
}

#' Get Cell Type Names for CL IDs
#'
#' Retrieve official cell type names from Cell Ontology for given CL IDs.
#'
#' @param clids Character vector of CL IDs
#' @param cl_ontology ontologyIndex CL object from \code{ontologyIndex::get_ontology()}
#' @return Named character vector of cell type names
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#' names <- get_cl_names(c("CL:0000066", "CL:0000057"), cl)
#' print(names)
#' }
#' @export
get_cl_names <- function(clids, cl_ontology) {
  names_vec <- stats::setNames(rep(NA_character_, length(clids)), clids)

  for (i in seq_along(clids)) {
    clid <- clids[i]
    if (!is.na(clid) && clid %in% names(cl_ontology$name)) {
      names_vec[i] <- cl_ontology$name[[clid]]
    }
  }

  return(names_vec)
}
