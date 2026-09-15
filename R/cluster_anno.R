#' Parse Raw LLM Annotation Lines
#'
#' @param response_text Character. Raw response text returned by the LLM.
#'
#' @return Character vector of parsed annotation labels.
#' @keywords internal
#' @noRd
.parse_llm_annotation_lines <- function(response_text) {
  # Split into lines and remove empty/whitespace-only lines
  res_lines <- strsplit(response_text, "\n", fixed = TRUE)[[1]]
  res_lines <- res_lines[nzchar(trimws(res_lines))]

  # Remove leading numbering (e.g., '1. ') and bullets ('- ', '* ')
  res_tmp <- gsub("^\\s*[0-9]+\\.\\s*", "", res_lines)
  res_tmp <- gsub("^\\s*[-*]\\s+", "", res_tmp)
  res_tmp <- trimws(res_tmp)

  # Some providers return 'marker - celltype' pairs; keep the right-hand side.
  res_tmp <- vapply(res_tmp, function(line) {
    if (grepl(" - ", line, fixed = TRUE)) {
      parts <- strsplit(line, " - ", fixed = TRUE)[[1]]
      return(parts[length(parts)])
    }
    line
  }, character(1))

  # Drop only exact header-like lines, not valid predictions such as "mixed cell types".
  is_header_line <- grepl(
    "^(?:here are(?: the)?(?: predicted)? cell types:?|the cell types are:?|cell types:?|predictions?:?)$",
    res_tmp,
    ignore.case = TRUE,
    perl = TRUE
  )
  unname(trimws(res_tmp[!is_header_line]))
}

.normalize_marker_input <- function(input, topgenenumber) {
  if (is.list(input) && !is.data.frame(input)) {
    collapsed <- vapply(input, function(x) {
      if (length(x) == 0) {
        return(NA_character_)
      }
      paste(x, collapse = ",")
    }, character(1))

    if (is.null(names(collapsed))) {
      names(collapsed) <- as.character(seq_along(collapsed))
    }
  } else {
    filtered <- input[input$avg_log2FC > 1, , drop = FALSE]
    split_genes <- split(filtered$gene, filtered$cluster, drop = FALSE)

    collapsed <- vapply(split_genes, function(genes) {
      genes <- genes[seq_len(min(topgenenumber, length(genes)))]
      if (length(genes) == 0) {
        return(NA_character_)
      }
      paste0(genes, collapse = ",")
    }, character(1))
  }

  valid_mask <- !is.na(collapsed) & nzchar(trimws(collapsed))
  list(
    all_input = collapsed,
    valid_input = collapsed[valid_mask],
    empty_clusters = names(collapsed)[!valid_mask]
  )
}

.log_llm_mismatch <- function(response_text, batch_idx, parsed_count, expected_count,
                              cluster_names, prompt_text = NULL) {
  out_dir <- file.path(tempdir(), "gptanno_llm_raw")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  timestamp <- gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d_%H%M%OS3"))
  unique_tag <- paste0(timestamp, "_pid", Sys.getpid())
  cluster_stub <- paste(cluster_names, collapse = "-")
  cluster_stub <- gsub("[^A-Za-z0-9._-]", "_", cluster_stub)

  raw_file <- file.path(
    out_dir,
    paste0("raw_response_", unique_tag, "_batch_", batch_idx, "_clusters_", cluster_stub, ".txt")
  )
  aggregate_file <- file.path(out_dir, "raw_response_log.txt")

  tryCatch({
    writeLines(response_text, con = raw_file)
    message("Saved raw response to: ", raw_file)
  }, error = function(e) {
    message("Failed to write raw response to file: ", e$message)
  })

  log_lines <- c(
    paste0("=== ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ==="),
    paste0("batch: ", batch_idx),
    paste0("expected_count: ", expected_count),
    paste0("parsed_count: ", parsed_count),
    paste0("clusters: ", paste(cluster_names, collapse = ", ")),
    "raw_response:",
    response_text
  )

  if (!is.null(prompt_text)) {
    log_lines <- c(log_lines, "prompt:", prompt_text)
  }
  log_lines <- c(log_lines, "")

  tryCatch({
    write(log_lines, file = aggregate_file, append = TRUE)
    message("Appended mismatch log to: ", aggregate_file)
  }, error = function(e) {
    message("Failed to append mismatch log: ", e$message)
  })
}

.is_unknown_prediction_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- sub("^\\s*(-|[0-9]+\\.)\\s*", "", x, perl = TRUE)
  x <- sub("^[-\\s]+", "", x, perl = TRUE)
  x <- trimws(x)
  x[!grepl("[a-z]", x)] <- ""
  x %in% c("", "unknown")
}

.empty_failed_run_summary <- function() {
  data.frame(
    run = character(0),
    note = character(0),
    stringsAsFactors = FALSE
  )
}

#' Query LLM for Cell Type Annotation Using Marker Genes
#'
#' Calls an LLM to predict cell types for each cluster or subcluster using marker genes.
#' Supports multiple providers (OpenAI, Anthropic, Gemini, Ollama, and vLLM) via the ellmer package,
#' plus the Claude Code and Codex CLIs in headless mode (providers "claude_cli" / "codex_cli").
#' For each cluster, marker genes are first ranked by avg_log2FC, and the top genes are selected.
#' The prompt can optionally include a request for Cell Ontology prediction, restriction to a set of cell types,
#' and/or an explicit instruction to predict a child cell type of a parent.
#'
#' @param input A data frame or list of marker genes per cluster/subcluster.
#' @param tissue_name Character. Biological context (e.g., "human lung").
#' @param model Character. Model name (default: 'gpt-5'). Used if llm_config not provided.
#' @param topgenenumber Integer. Number of top marker genes per cluster (default: 10).
#' @param add_cl_prompt Logical. If TRUE, adds "Please try to predict cell types in Cell Ontology." to the prompt.
#' @param restrict_to Optional character vector. Restrict predictions to these cell types.
#' @param parent_celltype Optional character. Predict child cell type of this parent.
#' @param llm_config Optional list. LLM configuration with provider, model, params, api_key,
#'   system_prompt, and provider-specific fields. Overrides `model` when supplied.
#'
#' @return A named character vector of predicted cell types for each cluster.
#' @importFrom ellmer params
#' @importFrom magrittr %>%
#' @export
gptcelltype <- function(input, tissue_name = NULL, model = 'gpt-5', topgenenumber = 10,
                        add_cl_prompt = FALSE, restrict_to = NULL, parent_celltype = NULL,
                        llm_config = NULL) {
  # Prepare LLM configuration
  config <- prepare_config(model = model, llm_config = llm_config)

  # Normalize input to named character vector of markers per cluster
  normalized_input <- .normalize_marker_input(input, topgenenumber = topgenenumber)
  input <- normalized_input$all_input
  valid_input <- normalized_input$valid_input
  failure_notes <- character(0)

  if (length(normalized_input$empty_clusters) > 0) {
    total_clusters <- length(input)
    empty_count <- length(normalized_input$empty_clusters)
    valid_count <- length(valid_input)
    message(
      "Clusters with no usable marker genes after filtering: ",
      paste(normalized_input$empty_clusters, collapse = ", "),
      " (", empty_count, "/", total_clusters, " cluster(s)); ",
      valid_count, " cluster(s) will be sent to the LLM and the empty cluster(s) will remain 'unknown'."
    )
  }

  # Build prompt
  base_prompt <- paste0(
    'Identify cell types of ', tissue_name, ' cells using the following markers separately for each\nrow. ',
    'Only provide the cell type name. Do not show numbers before the name.\n'
  )
  if (is.null(parent_celltype)) {
    base_prompt <- paste0(
      base_prompt,
      'Some clusters can be a mixture of multiple cell types.\n',
      'If a cluster is a mixture, provide the full cell type name for each component and separate them only with |.\n',
      'Do not use parentheses, brackets, commas, slashes, semicolons, plus signs, ampersands, or explanations.\n',
      'If the cluster is not mixed, provide only one cell type name.\n'
    )
  } else {
    base_prompt <- paste0(
      base_prompt,
      "For each subcluster, predict a more specific (child) cell type of the parent cell type: ",
      parent_celltype, ".\n",
      'Do not use parentheses, brackets, commas, slashes, semicolons, plus signs, ampersands, or explanations.\n'
    )
  }
  if (add_cl_prompt) {
    base_prompt <- paste0(base_prompt, "Please try to predict cell types in Cell Ontology.\n")
  }
  if (!is.null(restrict_to) && length(restrict_to) > 0) {
    base_prompt <- paste0(
      base_prompt,
      "Restrict your predictions to the following cell types:\n",
      paste(restrict_to, collapse = ", "), "\n"
    )
  }

  # Initialize result vector with "unknown" as default for all clusters
  res <- rep("unknown", length(input))
  names(res) <- names(input)

  if (length(valid_input) == 0) {
    message("No usable marker genes remain after filtering. Returning 'unknown' for all clusters.")
    attr(res, "run_note") <- "no usable marker genes after filtering"
    return(res)
  }

  # Batch requests: process 50 clusters per API call
  batch_size <- 50
  num_batches <- ceiling(length(valid_input) / batch_size)
  batch_ids <- if (num_batches > 1) {
    as.numeric(cut(seq_along(valid_input), num_batches))
  } else {
    rep(1, length(valid_input))
  }

  allres <- sapply(seq_len(num_batches), function(batch_idx) {
    cluster_indices <- which(batch_ids == batch_idx)
    batch_input <- valid_input[cluster_indices]
    batch_cluster_names <- names(valid_input)[cluster_indices]

    # Build prompt for this batch
    batch_prompt <- paste0(base_prompt, paste(batch_input, collapse = '\n'))

    batch_res <- rep("unknown", length(cluster_indices))

    tryCatch({
      # Build params object from llm_config$params (preferred) or legacy fields
      params_obj <- NULL
      if (!is.null(config$params)) {
        params_obj <- config$params
      } else if (!is.null(config$temperature) || !is.null(config$max_tokens)) {
        params_list <- list()
        if (!is.null(config$temperature)) params_list$temperature <- config$temperature
        if (!is.null(config$max_tokens)) params_list$max_tokens <- config$max_tokens
        params_obj <- do.call(ellmer::params, params_list)
      }

      response_text <- call_llm(
        prompt = batch_prompt,
        provider = config$provider,
        model = config$model,
        params = params_obj,
        api_key = config$api_key,
        system_prompt = config$system_prompt,
        api_url = config$api_url
      )

      # Parse provider output into one candidate label per line.
      res_tmp <- .parse_llm_annotation_lines(response_text)
      # Line-count validation: ensure response matches number of clusters in batch
      if (length(res_tmp) == length(cluster_indices)) {
        batch_res <- res_tmp
      } else {
        failure_notes <- unique(c(
          failure_notes,
          sprintf("batch %d response line count mismatch", batch_idx)
        ))
        message(sprintf("Batch %d: Response lines (%d) != expected (%d). Using 'unknown' for unmatched clusters.",
                        batch_idx, length(res_tmp), length(cluster_indices)))
        # Keep batch_res as "unknown" for all
          message(sprintf("Batch %d: Response lines (%d) != expected (%d). Writing raw response to console and file.",
                          batch_idx, length(res_tmp), length(cluster_indices)))
          cat("\n--- RAW LLM RESPONSE (batch ", batch_idx, ") ---\n", sep = "")
          cat(response_text, "\n")
          cat("--- END RAW LLM RESPONSE ---\n")
          .log_llm_mismatch(
            response_text = response_text,
            batch_idx = batch_idx,
            parsed_count = length(res_tmp),
            expected_count = length(cluster_indices),
            cluster_names = batch_cluster_names,
            prompt_text = batch_prompt
          )
          # Keep batch_res as "unknown" for all
      }
    }, error = function(e) {
      err_msg <- conditionMessage(e)
      failure_notes <<- unique(c(
        failure_notes,
        sprintf("batch %d LLM call failed: %s", batch_idx, err_msg)
      ))
      message(sprintf("LLM call failed for batch %d: %s", batch_idx, err_msg))
      if (identical(config$provider, "openai") && grepl("HTTP 400", err_msg, fixed = TRUE)) {
        message("Hint: OpenAI returned HTTP 400, which often indicates model/parameter incompatibility.")
        message("Hint: If the message includes 'Unsupported parameter', remove that field from llm_config$params")
        message("      or switch to a model that supports it.")
      }
      message("Prompt was:\n", batch_prompt)
      # batch_res remains as "unknown" for all
    })

    names(batch_res) <- batch_cluster_names
    batch_res
  }, simplify = FALSE)

  # Merge all batch results
  res[names(unlist(allres))] <- unlist(allres)
  res <- gsub(',$', '', res)
  attr(res, "run_note") <- if (length(failure_notes) > 0) {
    paste(failure_notes, collapse = "; ")
  } else {
    NA_character_
  }

  return(res)
}

#' Summarize LLM Cell Type Annotations Across Multiple Runs
#'
#' Calls LLM multiple times for cell type annotation and summarizes results.
#'
#' @param markers Named character vector or list of marker genes per cluster.
#' @param model Character. Model to use (default: 'gpt-5').
#' @param tissue_name Character. Optional context for prompt.
#' @param n_runs Integer. Number of LLM calls to aggregate (default: 2).
#' @param topgenenumber Integer. Number of top marker genes per cluster (default: 10).
#' @param add_cl_prompt Logical. If TRUE, adds Cell Ontology request to prompt.
#' @param restrict_to Optional character vector. Restrict predictions to these cell types.
#' @param parent_celltype Optional character. Predict child cell type of this parent.
#' @param llm_config Optional list. LLM configuration passed to gptcelltype.
#'   Supports OpenAI, Anthropic, and Google Gemini providers.
#'
#' @return List with: \itemize{
#'   \item combined_results: raw predictions from successful runs,
#'   \item summary: weighted summary table,
#'   \item final_summary: most frequent annotation per cluster,
#'   \item run_summary: failed runs with note column.
#' }
#' @importFrom dplyr bind_rows mutate group_by ungroup summarize arrange
#' @importFrom tidyr pivot_longer unnest
#' @importFrom stringr str_to_lower str_remove_all str_replace str_trim str_extract
#' @export
summarize_gptcelltype <- function(markers, model = 'gpt-5', tissue_name = "", n_runs = 2, topgenenumber = 10,
                                  add_cl_prompt = FALSE, restrict_to = NULL, parent_celltype = NULL,
                                  llm_config = NULL) {
  results_list <- vector("list", n_runs)
  successful_run_ids <- character(0)
  run_summary_rows <- list()
  run_summary_idx <- 1L

  for (i in seq_len(n_runs)) {
    res <- gptcelltype(
      markers,
      model = model,
      tissue_name = tissue_name,
      topgenenumber = topgenenumber,
      add_cl_prompt = add_cl_prompt,
      restrict_to = restrict_to,
      parent_celltype = parent_celltype,
      llm_config = llm_config
    )
    run_failed <- length(res) > 0 && all(.is_unknown_prediction_label(res))

    if (run_failed) {
      note <- attr(res, "run_note", exact = TRUE)
      if (is.null(note) || is.na(note) || !nzchar(note)) {
        note <- "all clusters predicted as unknown"
      }
      message("Excluding failed run ", i, ": ", note)
      run_summary_rows[[run_summary_idx]] <- data.frame(
        run = as.character(i),
        note = note,
        stringsAsFactors = FALSE
      )
      run_summary_idx <- run_summary_idx + 1L
    } else {
      results_list[[length(successful_run_ids) + 1L]] <- res
      successful_run_ids <- c(successful_run_ids, as.character(i))
    }
  }

  run_summary <- if (length(run_summary_rows) > 0) {
    do.call(rbind, run_summary_rows)
  } else {
    .empty_failed_run_summary()
  }

  if (length(successful_run_ids) == 0) {
    return(list(
      combined_results = data.frame(run = character(0), stringsAsFactors = FALSE),
      summary = data.frame(
        cluster = character(0),
        annotation_split = character(0),
        total_weight = numeric(0),
        total = numeric(0),
        percentage = numeric(0),
        stringsAsFactors = FALSE
      ),
      final_summary = data.frame(
        cluster = character(0),
        most_frequent_annotation = character(0),
        max_percentage = numeric(0),
        other_annotations = character(0),
        stringsAsFactors = FALSE
      ),
      run_summary = run_summary
    ))
  }

  results_list <- results_list[seq_along(successful_run_ids)]
  names(results_list) <- successful_run_ids
  combined_results <- dplyr::bind_rows(results_list, .id = "run")
  split_results <- combined_results %>%
    tidyr::pivot_longer(cols = -run, names_to = "cluster", values_to = "annotation") %>%
    dplyr::mutate(
      annotation = stringr::str_to_lower(annotation),
      annotation = stringr::str_remove_all(annotation, "^\\s*(-|\\d+\\.)\\s*"),
      annotation = stringr::str_replace(annotation, "^[-\\s]+", ""),
      annotation = stringr::str_trim(annotation, side = "right"),
      annotation = stringr::str_extract(annotation, "[a-zA-Z].*"),
      standardized_annotation = sapply(annotation, function(x) {
        if (is.na(x)) {
          return(NA_character_)
        } else {
          clean_and_match_annotation(x)
        }
      })
    ) %>%
    dplyr::mutate(annotation_split = strsplit(as.character(standardized_annotation), "\\s*\\|\\s*")) %>%
    tidyr::unnest(annotation_split) %>%
    dplyr::group_by(cluster, run) %>%
    dplyr::mutate(weight = 1 / length(annotation_split)) %>%
    dplyr::ungroup()
  summary <- split_results %>%
    dplyr::group_by(cluster, annotation_split) %>%
    dplyr::summarize(total_weight = sum(weight), .groups = 'drop') %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(total = sum(total_weight),
                  percentage = (total_weight / total) * 100) %>%
    dplyr::arrange(cluster, dplyr::desc(percentage))
  final_summary <- summary %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      most_frequent_annotation = annotation_split[1],
      max_percentage = round(percentage[1], 2),
      other_annotations = if (n() > 1) {
        paste(annotation_split[-1], round(percentage[-1], 2), "%", collapse = ", ")
      } else { "" },
      .groups = 'drop'
    )
  return(list(
    combined_results = combined_results,
    summary = summary,
    final_summary = final_summary,
    run_summary = run_summary
  ))
}



#' Run Annotation Workflow for All Cluster Resolutions
#'
#' Applies the OntoAnno LLM annotation workflow for all specified clustering resolutions.
#' For each, retrieves markers, runs annotation, calculates ontology distances, and stores output.
#'
#' @param seurat_obj A Seurat object.
#' @param resolutions Numeric vector of resolutions to annotate.
#' @param cl The ontology object.
#' @param graph An igraph object representing the Cell Ontology DAG. This can be created using:
#'   \code{cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo", extract_tags = "everything")}
#'   and \code{graph <- build_ontology_graph(cl)}.
#' @param mapping_dict Named character vector or data.frame; mapping from GPT-predicted to CL names. Defaults to package data \code{GPTCelltype_mapping}.
#' @param model Character. Model to use (default: 'gpt-5').
#' @param tissue_name Character. Optional context for prompt.
#' @param n_runs Integer. Number of LLM calls to aggregate (default: 2).
#' @param topgenenumber Integer. Number of top marker genes per cluster to include (default: 10).
#'   This parameter is passed to \code{summarize_gptcelltype()} and \code{gptcelltype()} to control
#'   how many of the top marker genes are used for each cluster in the annotation process.
#' @param add_cl_prompt Logical. If TRUE, adds "Please try to predict cell types in Cell Ontology." to the prompt.
#' @param marker_dir Character. Directory containing marker gene files (default: "output/marker_genes").
#'   Marker files should be named "markers_res_<resolution>.rds".
#' @param save_plots Logical. Whether to save comparison plots (default: FALSE).
#' @param plot_dir Character. Directory to save plots (default: "./output/prediction").
#'   Only used if save_plots = TRUE.
#' @param llm_config Optional list. LLM configuration passed to summarize_gptcelltype.
#'   Supports OpenAI, Anthropic, and Google Gemini providers.
#'
#' @return A named list of annotation summary objects for each resolution.
#' @importFrom dplyr arrange
#' @section Renamed function:
#' `gptanno()` is retained as a deprecated compatibility wrapper. New code
#' should use `ontoanno()`.
#'
#' @export
ontoanno <- function(seurat_obj, resolutions, cl, graph,
                    mapping_dict = GPTCelltype_mapping,
                    model = 'gpt-5',
                    tissue_name = NULL,
                    n_runs = 2,
                    topgenenumber = 10,
                    add_cl_prompt = FALSE,
                    marker_dir = "output/marker_genes",
                    save_plots = FALSE,
                    plot_dir = "./output/prediction",
                    llm_config = NULL) {
  results_list <- list()

  # Create plot directory only if saving plots
  if (save_plots && !dir.exists(plot_dir)) {
    message("Creating plot directory at: ", plot_dir)
    dir.create(plot_dir, recursive = TRUE)
  }

  for (res in resolutions) {
    message("\nRunning annotation for resolution: ", res)
    col_name <- paste0("cluster_res.", res)
    if (!col_name %in% colnames(seurat_obj@meta.data)) {
      warning("Column ", col_name, " not found in metadata. Skipping.")
      next
    }
    Seurat::Idents(seurat_obj) <- col_name

    # Construct marker file path using marker_dir parameter
    marker_file <- file.path(marker_dir, paste0("markers_res_", res, ".rds"))
    if (!file.exists(marker_file)) {
      warning("Marker file not found: ", marker_file, ". Skipping.")
      next
    }

    markers <- readRDS(marker_file)
    annotation_summary <- summarize_gptcelltype(
      markers,
      model = model,
      tissue_name = tissue_name,
      n_runs = n_runs,
      topgenenumber = topgenenumber,
      add_cl_prompt = add_cl_prompt,
      llm_config = llm_config
    )

    all_clusters <- unique(seurat_obj@meta.data[[col_name]])
    annotated_clusters <- unique(annotation_summary$summary$cluster)
    missing_clusters <- setdiff(all_clusters, annotated_clusters)
    if (length(missing_clusters) > 0) {
      warning("The following clusters lack annotations and will be labeled 'unannotated': ",
              paste(missing_clusters, collapse = ", "))
    }

    annotation_summary <- calculate_ontology_distance(
      annotation_summary,
      ontology_graph = graph,
      cl_term_map = OntoAnno::cl_term_map
    )

    annotated_seurat <- assign_celltype(seurat_obj, annotation_summary)
    results_list[[paste0("res_", res)]] <- annotation_summary

    # Save plots only if requested
    if (save_plots) {
      grDevices::pdf(file = file.path(plot_dir, paste0("res_", res, ".pdf")),
          width = 30, height = 10)
      print(plot_celltype_comparison(annotated_seurat))
      grDevices::dev.off()
      message("Plot saved to: ", file.path(plot_dir, paste0("res_", res, ".pdf")))
    }
  }
  return(results_list)
}

#' @rdname ontoanno
#' @export
gptanno <- function(seurat_obj, resolutions, cl, graph,
                    mapping_dict = GPTCelltype_mapping,
                    model = 'gpt-5',
                    tissue_name = NULL,
                    n_runs = 2,
                    topgenenumber = 10,
                    add_cl_prompt = FALSE,
                    marker_dir = "output/marker_genes",
                    save_plots = FALSE,
                    plot_dir = "./output/prediction",
                    llm_config = NULL) {
  .Deprecated("ontoanno", package = "OntoAnno")
  ontoanno(
    seurat_obj = seurat_obj,
    resolutions = resolutions,
    cl = cl,
    graph = graph,
    mapping_dict = mapping_dict,
    model = model,
    tissue_name = tissue_name,
    n_runs = n_runs,
    topgenenumber = topgenenumber,
    add_cl_prompt = add_cl_prompt,
    marker_dir = marker_dir,
    save_plots = save_plots,
    plot_dir = plot_dir,
    llm_config = llm_config
  )
}

#' Assign Annotated Cell Types to Seurat Object Metadata
#'
#' Adds a new metadata column to a Seurat object with annotated cell types for each cluster.
#'
#' @param seurat_obj A Seurat object.
#' @param annotation_summary Annotation summary as returned by `summarize_gptcelltype`.
#' @param cluster_col Character. Column to use as cluster identity (default: NULL, use Idents).
#' @param new_celltype Character. Name of new metadata column (default: "annotated_celltype").
#'
#' @return The updated Seurat object with new cell type annotations.
#' @importFrom Seurat Idents
#' @export
assign_celltype <- function(seurat_obj, annotation_summary, cluster_col = NULL, new_celltype = "annotated_celltype") {
  if (!is.null(cluster_col)) {
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("assign_celltype(): column '", cluster_col, "' not found in metadata.")
    }
    Seurat::Idents(seurat_obj) <- cluster_col
  }
  current_idents <- as.character(Seurat::Idents(seurat_obj))
  final_summary <- annotation_summary$final_summary
  if (is.null(final_summary) || !"cluster" %in% colnames(final_summary))
    stop("assign_celltype(): annotation_summary$final_summary is missing or malformed.")
  cluster_annotations <- stats::setNames(final_summary$most_frequent_annotation, as.character(final_summary$cluster))
  all_clusters    <- unique(current_idents)
  missing_clusters <- setdiff(all_clusters, names(cluster_annotations))
  if (length(missing_clusters) > 0) {
    warning("Unannotated clusters: ", paste(missing_clusters, collapse = ", "),
            " labeled 'unannotated'.")
    cluster_annotations[missing_clusters] <- "unannotated"
  }
  seurat_obj@meta.data[[new_celltype]] <- cluster_annotations[current_idents]
  return(seurat_obj)
}

#' Calculate Average Ontology Distance Between Cluster Annotations (using CL term map)
#'
#' Maps annotation names (including synonyms) to CL IDs using `cl_term_map` and computes mean ontology distances for each cluster.
#'
#' @param result_summary List as returned by `summarize_gptcelltype` (must include `summary` and `final_summary`).
#' @param ontology_graph An igraph object representing the Cell Ontology DAG.
#' @param cl_term_map Data frame produced by `build_cl_term_map()`, with columns 'key', 'clid', 'cl_label'.
#'
#' @return The input `result_summary` list, with `final_summary` augmented by `avg_distance` per cluster.
#' @importFrom dplyr left_join mutate group_by ungroup summarise filter
#' @importFrom igraph V shortest_paths
#' @export
calculate_ontology_distance <- function(result_summary, ontology_graph, cl_term_map) {
  # Defensive checks
  if (!inherits(ontology_graph, "igraph")) {
    warning("Invalid ontology graph provided. Skipping distance calculation.")
    result_summary$final_summary$avg_distance <- NA
    return(result_summary)
  }
  if (is.null(result_summary$summary) || !"annotation_split" %in% colnames(result_summary$summary)) {
    stop("result_summary$summary is missing or does not contain 'annotation_split'.")
  }
  if (!"cluster" %in% colnames(result_summary$summary)) {
    stop("result_summary$summary does not contain a 'cluster' column.")
  }
  # Helper: map annotation (synonyms or names) to CL IDs using the term map
  name_to_id <- function(name) {
    i <- match(tolower(name), cl_term_map$key)
    if (is.na(i)) return(NA_character_)
    cl_term_map$clid[[i]]
  }
  mapped_summary <- dplyr::mutate(
    result_summary$summary,
    cl_term = sapply(annotation_split, name_to_id)
  )
  # Compute average ontology distance for each cluster
  distances <- mapped_summary %>%
    dplyr::filter(!is.na(cl_term)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      avg_distance = {
        terms <- unique(cl_term)
        if (length(terms) < 2) {
          NA_real_
        } else {
          dist_values <- utils::combn(terms, 2, function(x) {
            if (all(x %in% igraph::V(ontology_graph)$name)) {
              path_result <- tryCatch({
                igraph::shortest_paths(ontology_graph, from = x[1], to = x[2], mode = "all")$vpath[[1]]
              }, error = function(e) NULL)
              if (is.null(path_result) || length(path_result) == 0) return(Inf)
              return(length(path_result) - 1)
            } else {
              return(NA_real_)
            }
          }, simplify = TRUE)
          finite_vals <- dist_values[is.finite(dist_values)]
          if (length(finite_vals) == 0) NA_real_ else mean(finite_vals, na.rm = TRUE)
        }
      },
      .groups = "drop"
    )
  result_summary$final_summary <- dplyr::left_join(
    result_summary$final_summary,
    distances,
    by = "cluster"
  )
  return(result_summary)
}
