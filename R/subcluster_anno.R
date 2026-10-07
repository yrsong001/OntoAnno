# ============================================================================
# Helper Functions for Subcluster Annotation
# ============================================================================

#' Get All Descendant Names from an Ontology Term
#'
#' Given a term name, returns all descendant term names (recursive) from a provided ontology.
#' Uses `get_descendants()` from the \pkg{ontologyIndex} package.
#'
#' @param ontology An ontology object from \pkg{ontologyIndex}.
#' @param term_name Character. Ontology term name to find descendants for.
#'
#' @return Character vector of descendant term names (unique).
#' @importFrom ontologyIndex get_descendants
#' @export
get_all_descendant_names <- function(ontology, term_name) {
  clean_name <- gsub("_", " ", term_name)
  term_ids <- names(ontology$name[ontology$name == clean_name])

  if (length(term_ids) == 0) {
    stop("Term name '", term_name, "' not found in the ontology.")
  }

  all_descendant_names <- unlist(lapply(term_ids, function(term_id) {
    descendant_ids <- ontologyIndex::get_descendants(ontology, roots = term_id, exclude_roots = TRUE)
    ontology$name[descendant_ids]
  }))

  unique(unname(all_descendant_names))
}

#' Prepare Ontology-Based Annotation Strategy
#'
#' Configures parameters for ontology-based subcluster annotation.
#' Determines restriction list from ontology descendants, user input, or both.
#'
#' @param cl Ontology object from \pkg{ontologyIndex}
#' @param parent_celltype Character. Parent cell type name
#' @param user_restrict_to Optional character vector. User-specified cell type restrictions
#' @param combine_restrictions Logical. Controls how restrictions are determined:
#'   \itemize{
#'     \item If TRUE (default) and user_restrict_to is provided: Returns union of ontology descendants + user list
#'     \item If FALSE and user_restrict_to is provided: Returns ONLY user list (ignores descendants)
#'     \item If user_restrict_to is NULL: Always returns ontology descendants (combine_restrictions ignored)
#'   }
#'
#' @return List with strategy configuration: strategy, restrict_to, parent_celltype, tissue_context, marker_prep_fn
#' @export
#'
#' @examples
#' \dontrun{
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo")
#'
#' # Use only ontology descendants
#' config <- prepare_ontology_strategy(cl, "T cell")
#'
#' # Combine descendants with user list
#' config <- prepare_ontology_strategy(
#'   cl, "T cell",
#'   user_restrict_to = c("regulatory T cell", "gamma-delta T cell"),
#'   combine_restrictions = TRUE
#' )
#'
#' # Use ONLY user list (ignore descendants)
#' config <- prepare_ontology_strategy(
#'   cl, "T cell",
#'   user_restrict_to = c("regulatory T cell", "effector T cell"),
#'   combine_restrictions = FALSE
#' )
#' }
prepare_ontology_strategy <- function(cl,
                                      parent_celltype,
                                      user_restrict_to = NULL,
                                      combine_restrictions = TRUE) {

  # Get ontology descendants
  descendants <- tryCatch({
    get_all_descendant_names(cl, parent_celltype)
  }, error = function(e) {
    warning("Could not get descendants for ", parent_celltype, ": ", e$message)
    character(0)
  })

  # Clean and standardize user-provided cell type names
  if (!is.null(user_restrict_to) && length(user_restrict_to) > 0) {
    user_restrict_to_cleaned <- sapply(user_restrict_to, function(x) {
      clean_and_match_annotation(x)
    }, USE.NAMES = FALSE)

    # Report any changes from cleaning
    changed_idx <- which(user_restrict_to != user_restrict_to_cleaned)
    if (length(changed_idx) > 0) {
      message("  Standardized user-provided cell type names:")
      for (idx in changed_idx) {
        message("    '", user_restrict_to[idx], "' -> '", user_restrict_to_cleaned[idx], "'")
      }
    }
    user_restrict_to <- user_restrict_to_cleaned
  }

  # Determine final restrict_to list
  if (is.null(user_restrict_to)) {
    # Use only ontology descendants
    restrict_to <- descendants
    message("  Using ontology descendants only: ", length(descendants), " cell types")
  } else if (combine_restrictions && length(descendants) > 0) {
    # Combine: union of descendants and user list (unique values)
    restrict_to <- unique(c(descendants, user_restrict_to))
    message("  Combining restrictions (union):")
    message("    - Ontology descendants: ", length(descendants), " cell types")
    message("    - User-specified (cleaned): ", length(user_restrict_to), " cell types")
    message("    - Union (unique): ", length(restrict_to), " cell types")
  } else {
    # Use only user-specified list
    restrict_to <- user_restrict_to
    message("  Using user-specified restrictions only: ",
            if (!is.null(user_restrict_to)) length(user_restrict_to) else 0,
            " cell types")
  }

  message("  Final restriction list for ", parent_celltype, ": ",
          if (length(restrict_to) > 0) paste(length(restrict_to), "cell types") else "none")

  list(
    strategy = "ontology",
    restrict_to = if (length(restrict_to) > 0) restrict_to else NULL,
    parent_celltype = NULL,  # Not used in prompt for ontology strategy
    tissue_context = parent_celltype,
    marker_prep_fn = function(subcluster_markers, ...) {
      # Return subcluster markers as-is for ontology strategy
      subcluster_markers
    }
  )
}

#' Prepare Marker-Inheritance Annotation Strategy
#'
#' Configures parameters for marker-inheritance based subcluster annotation.
#' Combines parent cluster markers with subcluster markers.
#'
#' @param parent_marker_file Path to parent marker RDS file
#' @param parent_celltype Character. Parent cell type name
#' @param user_restrict_to Optional character vector. User-specified cell type restrictions
#' @param topgenenumber Integer. Number of top subcluster genes (default: 10)
#' @param parent_topgenenumber Integer. Number of top parent-cluster genes (default: same as topgenenumber)
#' @param min_log2fc Numeric or NULL. avg_log2FC a gene must exceed (default 1); NULL = no fold-change filter, genes taken in
#'   FindAllMarkers order
#' @param max_genes Integer or NULL. Cap on the combined (unique) gene list (default 20); NULL = no cap
#'
#' @return List with strategy configuration
#' @export
#'
#' @examples
#' \dontrun{
#' config <- prepare_marker_inheritance_strategy(
#'   parent_marker_file = "output/markers/markers_res_0.3.rds",
#'   parent_celltype = "T cell",
#'   topgenenumber = 10
#' )
#' }
prepare_marker_inheritance_strategy <- function(parent_marker_file,
                                                parent_celltype,
                                                user_restrict_to = NULL,
                                                topgenenumber = 10,
                                                parent_topgenenumber = topgenenumber,
                                                min_log2fc = 1,
                                                max_genes = 20) {

  if (!file.exists(parent_marker_file)) {
    stop("Parent marker file not found: ", parent_marker_file)
  }

  parent_markers <- readRDS(parent_marker_file)

  # Clean and standardize user-provided cell type names
  if (!is.null(user_restrict_to) && length(user_restrict_to) > 0) {
    user_restrict_to_cleaned <- sapply(user_restrict_to, function(x) {
      clean_and_match_annotation(x)
    }, USE.NAMES = FALSE)

    # Report any changes from cleaning
    changed_idx <- which(user_restrict_to != user_restrict_to_cleaned)
    if (length(changed_idx) > 0) {
      message("  Standardized user-provided cell type names:")
      for (idx in changed_idx) {
        message("    '", user_restrict_to[idx], "' -> '", user_restrict_to_cleaned[idx], "'")
      }
    }
    user_restrict_to <- user_restrict_to_cleaned
  }

  list(
    strategy = "marker_inheritance",
    restrict_to = user_restrict_to,  # Can be NULL or user-specified (cleaned)
    parent_celltype = parent_celltype,  # Used in GPT prompt
    tissue_context = NULL,
    parent_markers = parent_markers,  # Store for later use
    topgenenumber = topgenenumber,
    parent_topgenenumber = parent_topgenenumber,
    min_log2fc = min_log2fc,
    max_genes = max_genes,
    marker_prep_fn = function(subcluster_markers,
                              subset_seurat,
                              subcluster_col,
                              original_cluster_col) {
      # This will be called to combine markers
      combine_parent_subcluster_markers_internal(
        parent_markers = parent_markers,
        subcluster_markers = subcluster_markers,
        subset_seurat = subset_seurat,
        subcluster_col = subcluster_col,
        original_cluster_col = original_cluster_col,
        topgenenumber = topgenenumber,
        parent_topgenenumber = parent_topgenenumber,
        min_log2fc = min_log2fc,
        max_genes = max_genes
      )
    }
  )
}

#' Combine Parent and Subcluster Markers (Internal)
#'
#' @keywords internal
combine_parent_subcluster_markers_internal <- function(parent_markers,
                                                       subcluster_markers,
                                                       subset_seurat,
                                                       subcluster_col,
                                                       original_cluster_col,
                                                       topgenenumber = 10,
                                                       parent_topgenenumber = topgenenumber,
                                                       min_log2fc = 1,
                                                       max_genes = 20) {

  # Get unique subclusters
  sub_ids <- unique(subset_seurat@meta.data[[subcluster_col]])

  # Build combined marker list for each subcluster
  combined_list <- list()

  for (sub_id in sub_ids) {

    # Find cells in this subcluster
    cells_in_subcluster <- rownames(subset_seurat@meta.data)[
      subset_seurat@meta.data[[subcluster_col]] == sub_id
    ]

    # Find parent cluster mode for this subcluster
    parent_clusters <- subset_seurat@meta.data[cells_in_subcluster, original_cluster_col]
    parent_cluster_mode <- as.numeric(names(sort(table(parent_clusters), decreasing = TRUE))[1])

    # Get top genes from parent (optional avg_log2FC filter, keep original order = FindAllMarkers ranking)
    parent_top <- parent_markers %>%
      dplyr::filter(cluster == parent_cluster_mode, if (is.null(min_log2fc)) TRUE else avg_log2FC > min_log2fc) %>%
      dplyr::slice(seq_len(parent_topgenenumber)) %>%
      dplyr::pull(gene)

    # Get top genes from subcluster (optional avg_log2FC filter, keep original order = FindAllMarkers ranking)
    sub_top <- subcluster_markers %>%
      dplyr::filter(cluster == sub_id, if (is.null(min_log2fc)) TRUE else avg_log2FC > min_log2fc) %>%
      dplyr::slice(seq_len(topgenenumber)) %>%
      dplyr::pull(gene)

    # Combine, keeping unique, prioritizing subcluster markers
    combined <- unique(c(sub_top, parent_top))
    if (!is.null(max_genes)) combined <- utils::head(combined, max_genes)  # cap (default 20); NULL = keep all

    combined_list[[as.character(sub_id)]] <- paste(combined, collapse = ",")
  }

  return(combined_list)
}

#' Annotate Subclusters with Unified Strategy
#'
#' Core annotation function that works for both ontology and marker-inheritance strategies.
#' Uses gptcelltype() and summarize_gptcelltype() from cluster_anno.R.
#'
#' @param seurat_obj Subcluster Seurat object
#' @param strategy_config List returned by prepare_*_strategy() functions
#' @param resolutions Numeric vector of resolutions to annotate
#' @param marker_dir Directory containing marker RDS files
#' @param model Character. GPT model (default: "gpt-5")
#' @param tissue_name Character. Base tissue context
#' @param n_runs Integer. Number of GPT runs to aggregate (default: 2)
#' @param topgenenumber Integer. Number of top marker genes per cluster (default: 10)
#' @param add_cl_prompt Logical. Add Cell Ontology prompt (default: FALSE)
#' @param ontology_graph Optional igraph for distance calculation. If provided, avg_distance
#'   will be calculated for use in scoring. Pass the result of build_ontology_graph(cl)
#' @param subcluster_prefix Character. Metadata column prefix (default: "subcluster_res.")
#' @param original_cluster_col Character. Column for parent cluster assignment (needed for marker_inheritance)
#' @param save_dir Optional directory to save results
#' @param save_plots Logical. Save comparison plots? (default: TRUE)
#' @param llm_config Optional list. LLM configuration forwarded to `summarize_gptcelltype()`.
#' @param min_log2fc Numeric or NULL. Fold-change filter for data-frame marker input (default 1); NULL disables it.
#'
#' @return List with annotation results per resolution
#' @importFrom Seurat Idents
#' @export
annotate_subclusters <- function(seurat_obj,
                                 strategy_config,
                                 resolutions,
                                 marker_dir,
                                 model = "gpt-4o",
                                 tissue_name = NULL,
                                 n_runs = 2,
                                 topgenenumber = 10,
                                 add_cl_prompt = FALSE,
                                 ontology_graph = NULL,
                                 subcluster_prefix = "subcluster_res.",
                                 original_cluster_col = NULL,
                                 save_dir = NULL,
                                 save_plots = TRUE,
                                 llm_config = NULL,
                                 min_log2fc = 1) {

  strategy <- strategy_config$strategy
  results_list <- list()

  for (res in resolutions) {
    message("\n--- Annotating resolution: ", res, " ---")

    # Construct column name and marker file path
    col_name <- paste0(subcluster_prefix, format(res, nsmall = 1))
    marker_file <- file.path(marker_dir, paste0("markers_res_", res, ".rds"))

    if (!file.exists(marker_file)) {
      warning("Marker file not found: ", marker_file, ". Skipping.")
      next
    }

    # Set identities
    if (!col_name %in% colnames(seurat_obj@meta.data)) {
      warning("Column ", col_name, " not found in metadata. Skipping.")
      next
    }
    Seurat::Idents(seurat_obj) <- col_name

    # Load markers
    markers <- readRDS(marker_file)

    # Prepare markers based on strategy
    if (strategy == "marker_inheritance") {
      # Need to combine parent + subcluster markers
      if (is.null(original_cluster_col)) {
        stop("original_cluster_col required for marker_inheritance strategy")
      }

      prepared_markers <- strategy_config$marker_prep_fn(
        subcluster_markers = markers,
        subset_seurat = seurat_obj,
        subcluster_col = col_name,
        original_cluster_col = original_cluster_col
      )

      parent_celltype <- strategy_config$parent_celltype
      full_tissue_name <- paste(parent_celltype, "in", tissue_name)

    } else {
      # Ontology strategy - use markers as-is
      prepared_markers <- markers
      parent_celltype <- NULL

      tissue_context <- strategy_config$tissue_context
      full_tissue_name <- paste(tissue_context, "in", tissue_name)
    }

    # Call unified annotation from cluster_anno.R
    annotation_summary <- summarize_gptcelltype(
      markers = prepared_markers,
      model = model,
      tissue_name = full_tissue_name,
      n_runs = n_runs,
      topgenenumber = topgenenumber,
      add_cl_prompt = add_cl_prompt,
      restrict_to = strategy_config$restrict_to,
      parent_celltype = parent_celltype,
      llm_config = llm_config,
      min_log2fc = min_log2fc
    )

    # Add ontology distance if graph provided
    if (!is.null(ontology_graph)) {
      annotation_summary <- calculate_ontology_distance(
        annotation_summary,
        ontology_graph = ontology_graph,
        cl_term_map = OntoAnno::cl_term_map
      )
    }

    # Assign to Seurat object
    annotation_summary$summary$cluster <- as.character(annotation_summary$summary$cluster)
    Seurat::Idents(seurat_obj) <- as.character(Seurat::Idents(seurat_obj))

    annotated_seurat <- assign_celltype(
      seurat_obj,
      annotation_summary,
      new_celltype = paste0("annotated_sub_", res)
    )

    # Save results if requested
    if (!is.null(save_dir)) {
      dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

      if (save_plots) {
        grDevices::pdf(file.path(save_dir, paste0("annotation_plot_res_", res, ".pdf")),
            width = 18, height = 10)
        print(plot_celltype_comparison(
          annotated_seurat,
          original_col = col_name,
          annotation_col = paste0("annotated_sub_", res)
        ))
        grDevices::dev.off()
      }

      saveRDS(annotation_summary,
              file.path(save_dir, paste0("annotation_summary_res_", res, ".rds")))
    }

    results_list[[paste0("res_", res)]] <- list(
      combined_results = annotation_summary$combined_results,
      summary = annotation_summary$summary,
      final_summary = annotation_summary$final_summary,
      seurat = annotated_seurat
    )
  }

  return(results_list)
}

# ============================================================================
# Main Workflow Functions
# ============================================================================

#' Subcluster and Find Markers for Each Major Cell Type (with Ontology Check)
#'
#' Performs subclustering for each predicted cell type in a Seurat object and saves marker gene files per subcluster.
#' Only subclusters cell types that have descendant terms in the Cell Ontology.
#'
#' @param seurat_obj A Seurat object.
#' @param cl Cell Ontology object for checking descendants.
#' @param predicted_celltype_column Metadata column with parent cell type labels.
#' @param cluster_col Metadata column with cluster assignment (default: "seurat_clusters").
#' @param output_dir Output directory for subcluster markers.
#' @param resolutions Numeric vector of resolutions to use for subclustering.
#' @param dims Numeric vector of dimensions for clustering.
#' @param assay Which assay to use (default: NULL).
#' @param reduction Character. Reduction to use for graph construction during
#'   subclustering (default: "pca").
#' @param min_cell_count Minimum number of cells to trigger subclustering (default: 10000).
#' @param celltypes_to_subcluster Character vector of specific cell types to subcluster (default: NULL).
#'   If provided, these cell types will be subclustered regardless of whether they have ontology descendants.
#'   If NULL, automatic selection is used with descendant checking.
#'
#' @return List of subcluster results per cell type, or NULL if no subclustering performed.
#' @importFrom dplyr group_by summarize filter n_distinct n arrange slice pull
#' @importFrom rlang sym
#' @export
#' @examples
#' \dontrun{
#' # Load your ontology
#' cl <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl.obo",
#'                                   extract_tags = "everything")
#' # Run subclustering (only for cell types with descendants)
#' results <- subcluster_and_find_markers(
#'   seurat_obj = my_seurat,
#'   cl = cl,
#'   predicted_celltype_column = "celltype_parent",
#'   resolutions = c(0.1, 0.2, 0.3),
#'   reduction = "pca"
#' )
#' }
subcluster_and_find_markers <- function(seurat_obj,
                                        cl,
                                        predicted_celltype_column = "new_celltype",
                                        cluster_col = "seurat_clusters",
                                        output_dir = "results/subcluster_markers",
                                        resolutions = c(0.1, 0.2, 0.3),
                                        dims = 1:30,
                                        assay = NULL,
                                        reduction = "pca",
                                        min_cell_count = 10000,
                                        celltypes_to_subcluster = NULL) {

  if (missing(cl) || is.null(cl)) {
    stop("Error: 'cl' argument (Cell Ontology object) is required but not provided.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- seurat_obj@meta.data

  # --- Selection logic ---
  if (!is.null(celltypes_to_subcluster)) {
    # User-specified cell types
    celltype_stats <- meta %>%
      dplyr::group_by(!!sym(predicted_celltype_column)) %>%
      dplyr::summarize(
        cluster_count = dplyr::n_distinct(!!sym(cluster_col)),
        cell_count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(!!sym(predicted_celltype_column) %in% celltypes_to_subcluster)
    message("Cell types selected by user: ", paste(celltype_stats[[predicted_celltype_column]], collapse = ", "))
  } else {
    # Automatic selection by basic criteria
    celltype_stats <- meta %>%
      dplyr::group_by(!!sym(predicted_celltype_column)) %>%
      dplyr::summarize(
        cluster_count = dplyr::n_distinct(!!sym(cluster_col)),
        cell_count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(cluster_count > 11 | cell_count > min_cell_count)
    message("Cell types meeting basic criteria: ", paste(celltype_stats[[predicted_celltype_column]], collapse = ", "))
  }

  if (nrow(celltype_stats) == 0) {
    message("No cell types selected for subclustering.")
    return(NULL)
  }

  # Filter by ontology descendants ONLY if user did NOT specify celltypes
  if (is.null(celltypes_to_subcluster)) {
    # Automatic selection - apply descendant check
    celltype_with_descendants <- character(0)
    celltype_without_descendants <- character(0)

    for (ct in celltype_stats[[predicted_celltype_column]]) {
      tryCatch({
        descendants <- get_all_descendant_names(cl, ct)
        if (length(descendants) > 0) {
          celltype_with_descendants <- c(celltype_with_descendants, ct)
          message(ct, " has ", length(descendants), " descendant terms will be subclustered")
        } else {
          celltype_without_descendants <- c(celltype_without_descendants, ct)
          message(ct, " has no descendant terms skipping subclustering")
        }
      }, error = function(e) {
        celltype_without_descendants <- c(celltype_without_descendants, ct)
        message(ct, " error checking descendants (", e$message, ") skipping subclustering")
      })
    }

    if (length(celltype_with_descendants) == 0) {
      message("No cell types have descendant terms. No subclustering will be performed.")
      return(NULL)
    }

    celltypes_to_process <- celltype_with_descendants
  } else {
    # User-specified - skip descendant check
    celltypes_to_process <- celltype_stats[[predicted_celltype_column]]
    message("User-specified cell types will be subclustered regardless of descendants")
  }

  message("\nCell types to be subclustered: ", paste(celltypes_to_process, collapse = ", "))

  results <- list()
  for (ct in celltypes_to_process) {
    message("\nProcessing: ", ct)
    ct_safe <- gsub("[^[:alnum:]_]+", "_", ct)
    ct_dir <- file.path(output_dir, ct_safe)
    dir.create(ct_dir, showWarnings = FALSE)

    subset_path <- file.path(ct_dir, "seurat_subset.rds")
    if (file.exists(subset_path)) {
      message(" Loading existing subset: ", subset_path)
      subset_obj <- readRDS(subset_path)
    } else {
      message(" Creating subset...")
      subset_obj <- subset(seurat_obj, subset = !!sym(predicted_celltype_column) == ct)
    }

    existing_marker_files <- list.files(ct_dir, pattern = "^markers_res_.*\\.rds$", full.names = TRUE)
    existing_resolutions <- gsub("markers_res_(.*)\\.rds", "\\1", basename(existing_marker_files))
    resolutions_to_run <- resolutions[!as.character(resolutions) %in% existing_resolutions]

    if (length(resolutions_to_run) == 0) {
      message(" Marker files exist for all resolutions. Will update metadata only.")
    } else {
      message(" Running marker detection for resolutions: ", paste(resolutions_to_run, collapse = ", "))
    }

    subcluster_results <- run_multi_resolution_clustering(
      seurat_obj = subset_obj,
      resolutions = resolutions_to_run,
      result_dir = ct_dir,
      dims = dims,
      assay = assay,
      reduction = reduction,
      group.by = "seurat_clusters"
    )

    # Update metadata with subcluster information
    all_resolutions <- union(names(subcluster_results$clusters), existing_resolutions)
    for (res in all_resolutions) {
      colname <- paste0("subcluster_res.", res)
      if (res %in% names(subcluster_results$clusters)) {
        subset_obj@meta.data[[colname]] <- subcluster_results$clusters[[res]]
      } else if (!colname %in% colnames(subset_obj@meta.data)) {
        warning(" Resolution ", res, " exists as a marker file but not in subset metadata. Skipping column.")
      }
    }

    saveRDS(subset_obj, file = subset_path)
    results[[ct]] <- subcluster_results
  }

  return(results)
}

# ============================================================================
# New Unified Workflow Functions
# ============================================================================

#' Run Complete Subcluster Annotation Workflow
#'
#' Master function that orchestrates the entire subcluster annotation process
#' for both ontology and marker-inheritance strategies.
#'
#' @param base_dir Directory with subcluster Seurat objects (each celltype in subdirectory)
#' @param strategy Character. "ontology" or "marker_inheritance"
#' @param cl Ontology object (required for both strategies). Used to calculate ontology distance
#'   for scoring, and to get descendants for ontology strategy.
#' @param parent_marker_root Directory for parent markers (required for marker_inheritance)
#' @param parent_res Character. Resolution of parent clustering (for marker_inheritance)
#' @param parent_cluster_col Character. Column name for parent clusters in subset metadata
#' @param user_restrict_to Named list or character vector. Cell type restrictions.
#'   \itemize{
#'     \item NULL: No user restrictions (uses ontology descendants for ontology strategy)
#'     \item Character vector: Applied to ALL celltypes being processed
#'     \item Named list: Per-celltype restrictions (RECOMMENDED for multiple celltypes)
#'   }
#'   Names are automatically cleaned/standardized using \code{clean_and_match_annotation()}.
#'   For ontology strategy, see \code{combine_restrictions} to control behavior.
#'   Example: \code{list("T cell" = c("regulatory T cell", "effector T cell"),
#'                       "B cell" = c("memory B cell", "naive B cell"))}
#' @param combine_restrictions Logical. For ontology strategy only. If TRUE (default),
#'   combines (union) user_restrict_to with ontology descendants. If FALSE, uses ONLY
#'   user_restrict_to and ignores ontology descendants. Has no effect on marker_inheritance strategy.
#' @param model Character. GPT model (default: "gpt-5")
#' @param tissue_name Character. Tissue context
#' @param resolutions Numeric vector. Resolutions to annotate
#' @param n_runs Integer. Number of GPT runs (default: 2)
#' @param topgenenumber Integer. Top genes per cluster (default: 10)
#' @param add_cl_prompt Logical. Add CL prompt (default: FALSE)
#' @param ontology_graph Optional igraph for distance calculation. If NULL and cl is provided,
#'   the graph will be automatically built using build_ontology_graph(cl)
#' @param save_results Logical. Save intermediate results? (default: TRUE)
#' @param select_best Logical. Auto-select best resolution using score_annotation_resolutions()? (default: TRUE)
#' @param output_dir Optional directory for outputs
#' @param celltypes_to_subcluster Optional character vector. Specific parent celltypes to process.
#'   If NULL, processes all directories in base_dir. Overrides any celltype filtering.
#' @param llm_config Optional list. LLM configuration forwarded to `annotate_subclusters()`.
#' @param min_log2fc Numeric or NULL. Fold-change filter on markers (default 1); NULL disables it.
#' @param parent_topgenenumber Integer or NULL. Parent-cluster genes for marker_inheritance (default: topgenenumber).
#' @param max_genes Integer or NULL. Cap on the combined gene list for marker_inheritance (default 20); NULL = no cap.
#'
#' @return List with: results (all resolutions for all celltypes), scores, best, strategy
#' @export
#'
#' @examples
#' \dontrun{
#' # Example 1: Per-celltype restrictions (RECOMMENDED for multiple celltypes)
#' workflow_results <- run_subcluster_annotation_workflow(
#'   base_dir = "output/subclusters",
#'   strategy = "ontology",
#'   cl = cl,
#'   user_restrict_to = list(
#'     "T cell" = c("regulatory T cell", "effector T cell", "memory T cell"),
#'     "B cell" = c("memory B cell", "naive B cell", "plasma cell"),
#'     "endothelial cell" = c("vascular endothelial cell", "lymphatic endothelial cell")
#'   ),
#'   combine_restrictions = TRUE,  # Union with ontology descendants
#'   tissue_name = "mouse heart",
#'   select_best = TRUE
#' )
#'
#' # Example 2: Single celltype with user restrictions only (no descendants)
#' workflow_results <- run_subcluster_annotation_workflow(
#'   base_dir = "output/subclusters",
#'   strategy = "ontology",
#'   cl = cl,
#'   celltypes_to_subcluster = "T cell",
#'   user_restrict_to = c("regulatory T cell", "effector T cell"),
#'   combine_restrictions = FALSE,  # Use ONLY specified types
#'   tissue_name = "mouse heart"
#' )
#'
#' # Example 3: Same restrictions for all celltypes (character vector)
#' workflow_results <- run_subcluster_annotation_workflow(
#'   base_dir = "output/subclusters",
#'   strategy = "ontology",
#'   cl = cl,
#'   user_restrict_to = c("fibroblast", "myofibroblast", "pericyte"),
#'   combine_restrictions = TRUE,
#'   tissue_name = "mouse heart"
#' )
#'
#' # Example 4: Marker inheritance strategy with per-celltype restrictions
#' workflow_results <- run_subcluster_annotation_workflow(
#'   base_dir = "output/subclusters",
#'   strategy = "marker_inheritance",
#'   cl = cl,  # Required for ontology distance calculation
#'   parent_marker_root = "output/markers",
#'   parent_res = "0.3",
#'   parent_cluster_col = "cluster_res.0.3",
#'   user_restrict_to = list(
#'     "T cell" = c("Treg", "Helper T cells", "cytotoxic T cell")
#'   ),
#'   tissue_name = "mouse heart"
#' )
#' # Note: Names like "Treg", "Helper T cells" are auto-cleaned to standard ontology terms
#' }
run_subcluster_annotation_workflow <- function(
    base_dir = "output/subclusters",
    strategy = c("ontology", "marker_inheritance"),
    cl = NULL,
    parent_marker_root = NULL,
    parent_res = "0.3",
    parent_cluster_col = "cluster_res.0.3",
    user_restrict_to = NULL,
    combine_restrictions = TRUE,
    model = "gpt-4o",
    tissue_name = NULL,
    resolutions = c(0.1, 0.3, 0.5),
    n_runs = 2,
    topgenenumber = 10,
    add_cl_prompt = FALSE,
    ontology_graph = NULL,
    save_results = TRUE,
    select_best = TRUE,
    output_dir = NULL,
    celltypes_to_subcluster = NULL,
    llm_config = NULL,
    min_log2fc = 1,
    parent_topgenenumber = NULL,
    max_genes = 20) {

  strategy <- match.arg(strategy)

  # Validate inputs
  if (is.null(cl)) {
    stop("Ontology object 'cl' is required for both strategies (used for ontology distance calculation)")
  }

  if (strategy == "marker_inheritance") {
    if (is.null(parent_marker_root)) {
      stop("'parent_marker_root' required for marker_inheritance strategy")
    }
    if (is.null(parent_cluster_col)) {
      stop("'parent_cluster_col' required for marker_inheritance strategy")
    }
  }

  # Build ontology graph if cl is provided and graph is not
  if (!is.null(cl) && is.null(ontology_graph)) {
    message("Building ontology graph from cl object...")
    ontology_graph <- build_ontology_graph(cl)
    message("Ontology graph built successfully")
  }

  # Find celltype directories
  all_celltype_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

  # Filter by celltypes_to_subcluster if provided
  if (!is.null(celltypes_to_subcluster)) {
    # Convert user-specified celltypes to safe names (spaces to underscores)
    safe_celltypes <- gsub("[^[:alnum:]_]+", "_", celltypes_to_subcluster)
    dir_basenames <- basename(all_celltype_dirs)

    # Match either original name or safe name
    matched_indices <- which(dir_basenames %in% celltypes_to_subcluster |
                               dir_basenames %in% safe_celltypes)

    if (length(matched_indices) == 0) {
      stop("None of the specified celltypes found in base_dir.\n",
           "  Looking for: ", paste(celltypes_to_subcluster, collapse = ", "), "\n",
           "  Available directories: ", paste(dir_basenames, collapse = ", "))
    }

    celltype_dirs <- all_celltype_dirs[matched_indices]
    message("\nProcessing user-specified celltypes: ",
            paste(basename(celltype_dirs), collapse = ", "))
  } else {
    celltype_dirs <- all_celltype_dirs
    message("\nProcessing all celltypes in base_dir: ", length(celltype_dirs), " found")
  }
  all_results <- list()
  all_scores <- list()
  best_selections <- list()

  for (ct_dir in celltype_dirs) {
    celltype <- basename(ct_dir)

    message("\n", strrep("=", 60))
    message("Processing: ", celltype)
    message("Strategy: ", strategy)
    message(strrep("=", 60))

    # Load subset
    seurat_path <- file.path(ct_dir, "seurat_subset.rds")
    if (!file.exists(seurat_path)) {
      warning("No Seurat subset for ", celltype, ". Skipping.")
      next
    }
    seurat_subset <- readRDS(seurat_path)

    # Get user restrictions for this celltype
    # Try to match celltype name (handle both "endothelial cell" and "endothelial_cell")
    celltype_restrict <- if (is.list(user_restrict_to)) {
      # Try exact match first
      if (celltype %in% names(user_restrict_to)) {
        user_restrict_to[[celltype]]
      } else {
        # Try converting celltype to original form (underscore to space)
        original_name <- gsub("_", " ", celltype)
        if (original_name %in% names(user_restrict_to)) {
          user_restrict_to[[original_name]]
        } else {
          NULL
        }
      }
    } else {
      user_restrict_to
    }

    # Prepare strategy config
    if (strategy == "ontology") {
      strategy_config <- prepare_ontology_strategy(
        cl = cl,
        parent_celltype = celltype,
        user_restrict_to = celltype_restrict,
        combine_restrictions = combine_restrictions
      )
      original_cluster_col <- NULL
    } else {
      # marker_inheritance
      parent_marker_file <- file.path(
        parent_marker_root,
        paste0("markers_res_", parent_res, ".rds")
      )
      strategy_config <- prepare_marker_inheritance_strategy(
        parent_marker_file = parent_marker_file,
        parent_celltype = celltype,
        user_restrict_to = celltype_restrict,
        topgenenumber = topgenenumber,
        parent_topgenenumber = if (is.null(parent_topgenenumber)) topgenenumber else parent_topgenenumber,
        min_log2fc = min_log2fc,
        max_genes = max_genes
      )
      original_cluster_col <- parent_cluster_col
    }

    # Set save directory
    if (save_results) {
      celltype_save_dir <- file.path(
        if (!is.null(output_dir)) output_dir else ct_dir,
        paste0("annotation_", strategy)
      )
    } else {
      celltype_save_dir <- NULL
    }

    # Run annotation
    annotation_results <- annotate_subclusters(
      seurat_obj = seurat_subset,
      strategy_config = strategy_config,
      resolutions = resolutions,
      marker_dir = ct_dir,
      model = model,
      tissue_name = tissue_name,
      n_runs = n_runs,
      topgenenumber = topgenenumber,
      add_cl_prompt = add_cl_prompt,
      ontology_graph = ontology_graph,
      original_cluster_col = original_cluster_col,
      save_dir = celltype_save_dir,
      llm_config = llm_config,
      min_log2fc = min_log2fc
    )

    # Store ALL resolution results (like parent workflow)
    all_results[[celltype]] <- annotation_results

    # Score resolutions if multiple exist
    if (select_best && length(annotation_results) > 1) {
      score_table <- score_annotation_resolutions(annotation_results)
      message("\nResolution scores for ", celltype, ":")
      print(score_table)

      best_res <- score_table$resolution[1]
      message("Best resolution: ", best_res,
              " (composite score: ", round(score_table$composite_score[1], 3), ")")

      all_scores[[celltype]] <- score_table
      best_selections[[celltype]] <- list(
        resolution = best_res,
        composite_score = score_table$composite_score[1]
      )
    } else {
      all_scores[[celltype]] <- NULL
      best_selections[[celltype]] <- NULL
    }
  }

  # Return structure aligned with parent workflow
  return(list(
    results = all_results,           # ALL resolutions for ALL celltypes
    scores = all_scores,             # Score tables per celltype
    best = best_selections,          # Best resolution metadata per celltype
    strategy = strategy
  ))
}

#' Assign Subcluster Annotations to Full Seurat Object
#'
#' Unified function that works for both ontology and marker-inheritance strategies.
#' Takes subcluster annotation results and assigns them back to the full Seurat object.
#'
#' @param full_seurat Full Seurat object (not subset)
#' @param workflow_results List returned by run_subcluster_annotation_workflow()
#'   Contains: results, scores, best, strategy
#' @param use_best Logical. If TRUE, use best resolution per celltype. If FALSE, user must specify resolution_map
#' @param resolution_map Optional named vector. Specify resolution per celltype (e.g., c("T cell" = "res_0.1"))
#' @param parent_column Character. Metadata column with parent cell type annotation
#' @param final_colname Character. Name for new annotation column
#' @param annotation_col_pattern Character. Pattern to find annotation column in subset (default: "^annotated_sub_")
#' @param clean_labels Logical. Remove " (not in GO)" suffix? (default: TRUE)
#'
#' @return Updated full Seurat object with subcluster annotations
#' @export
#'
#' @examples
#' \dontrun{
#' # Using auto-selected best resolutions
#' full_seurat <- assign_subcluster_annotations_to_full(
#'   full_seurat = my_seurat,
#'   workflow_results = workflow_results,
#'   use_best = TRUE,
#'   parent_column = "celltype_parent"
#' )
#'
#' # Using manual resolution selection
#' full_seurat <- assign_subcluster_annotations_to_full(
#'   full_seurat = my_seurat,
#'   workflow_results = workflow_results,
#'   use_best = FALSE,
#'   resolution_map = c("T cell" = "res_0.2", "B cell" = "res_0.1"),
#'   parent_column = "celltype_parent"
#' )
#' }
assign_subcluster_annotations_to_full <- function(
    full_seurat,
    workflow_results,
    use_best = TRUE,
    resolution_map = NULL,
    parent_column = "celltype_parent",
    final_colname = "celltype_subcluster",
    annotation_col_pattern = "^annotated_sub_",
    clean_labels = TRUE) {

  # Extract components from workflow results
  annotation_results <- workflow_results$results
  best_selections <- workflow_results$best

  # Initialize with parent annotations
  if (!parent_column %in% colnames(full_seurat@meta.data)) {
    stop("Parent column '", parent_column, "' not found in full Seurat object")
  }

  full_seurat@meta.data[[final_colname]] <- full_seurat@meta.data[[parent_column]]

  # Process each cell type
  for (celltype in names(annotation_results)) {

    celltype_results <- annotation_results[[celltype]]

    # Determine which resolution to use
    if (!is.null(resolution_map) && celltype %in% names(resolution_map)) {
      # User-specified resolution
      best_res <- resolution_map[[celltype]]
      message("Using user-specified resolution ", best_res, " for ", celltype)
    } else if (use_best && !is.null(best_selections[[celltype]])) {
      # Auto-selected best resolution
      best_res <- best_selections[[celltype]]$resolution
      score <- best_selections[[celltype]]$composite_score
      message("Using best resolution ", best_res, " for ", celltype,
              " (score: ", round(score, 3), ")")
    } else {
      # Use first available resolution
      best_res <- names(celltype_results)[1]
      message("Using first resolution ", best_res, " for ", celltype)
    }

    # Get the annotation object
    res_data <- celltype_results[[best_res]]

    if (is.null(res_data) || is.null(res_data$seurat)) {
      warning("No Seurat object found for ", celltype, " at ", best_res, ". Skipping.")
      next
    }

    # Find annotation column in subset
    subset_seurat <- res_data$seurat
    annot_cols <- grep(annotation_col_pattern,
                       colnames(subset_seurat@meta.data),
                       value = TRUE)

    if (length(annot_cols) == 0) {
      warning("No annotation column found for ", celltype, ". Skipping.")
      next
    }

    if (length(annot_cols) > 1) {
      warning("Multiple annotation columns found for ", celltype,
              ". Using first: ", annot_cols[1])
    }

    annot_col <- annot_cols[1]

    # Extract annotations
    cell_names <- colnames(subset_seurat)
    annotations <- as.character(subset_seurat@meta.data[[annot_col]])

    # Clean labels if requested
    if (clean_labels) {
      annotations <- gsub(" \\(not in GO\\)", "", annotations)
    }

    # Assign to full object (matching by cell names)
    matching_cells <- intersect(cell_names, colnames(full_seurat))

    if (length(matching_cells) == 0) {
      warning("No matching cells found for ", celltype, ". Skipping.")
      next
    }

    full_seurat@meta.data[matching_cells, final_colname] <- annotations[match(matching_cells, cell_names)]

    message("Assigned ", length(matching_cells), " cells for ", celltype)
  }

  # Clean parent labels if requested
  if (clean_labels) {
    full_seurat@meta.data[[final_colname]] <- gsub(
      " \\(not in GO\\)",
      "",
      full_seurat@meta.data[[final_colname]]
    )
  }

  return(full_seurat)
}

# ============================================================================
# DEPRECATED FUNCTIONS - For Backward Compatibility Only
# ============================================================================
#
# The functions below are deprecated and maintained only for backward compatibility.
# For new code, use the unified workflow functions above:
#   - run_subcluster_annotation_workflow() instead of anno_subcluster_*
#   - gptcelltype() instead of gptcelltype_sub()
#   - summarize_gptcelltype() instead of summarize_gptcelltype_sub()
#   - assign_subcluster_annotations_to_full() instead of assign_*_subcluster_annotations()
#
# ============================================================================

#' Summarize GPT Subcluster Annotations
#'
#' **DEPRECATED:** Use \code{\link{summarize_gptcelltype}} instead.
#'
#' Runs GPT-based annotation (optionally restricting to ontology descendants) multiple times and summarizes results.
#'
#' @param markers Named character vector or list of marker genes per subcluster.
#' @param model Character. GPT model to use.
#' @param tissue_name Character. Tissue context for prompt.
#' @param n_runs Integer. Number of times to query GPT.
#' @param restrict_to Character vector of allowable cell types (default: NULL).
#'
#' @return List with combined results, summary, and final_summary.
#' @keywords internal deprecated
#' @importFrom dplyr bind_rows mutate group_by ungroup summarize arrange
#' @importFrom tidyr pivot_longer unnest
#' @importFrom stringr str_to_lower str_remove_all str_replace str_trim str_extract
#' @export
summarize_gptcelltype_sub <- function(markers,
                                      model = 'gpt-5',
                                      tissue_name = "",
                                      n_runs = 2,
                                      restrict_to = NULL) {
  results_list <- vector("list", n_runs)
  for (i in seq_len(n_runs)) {
    res <- gptcelltype_sub(
      markers,
      model = model,
      tissue_name = tissue_name,
      restrict_to = restrict_to
    )
    results_list[[i]] <- res
  }
  combined_results <- dplyr::bind_rows(results_list, .id = "run")
  split_results <- combined_results %>%
    tidyr::pivot_longer(cols = -run,
                        names_to = "cluster",
                        values_to = "annotation") %>%
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
          # Only pass the annotation string, no extra arguments
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
      } else {
        ""
      },
      .groups = 'drop'
    )
  return(
    list(
      combined_results = combined_results,
      summary = summary,
      final_summary = final_summary
    )
  )
}

#' Query GPT for Cell Type Annotation of Subclusters Using Marker Genes
#'
#' **DEPRECATED:** Use \code{\link{gptcelltype}} instead.
#'
#' Calls OpenAI's GPT model to predict cell types for each subcluster using a marker-gene table.
#'
#' @param input A **data frame** of marker genes with (at least) the columns:
#'   \code{cluster}, \code{gene}, and \code{avg_log2FC}. Example columns:
#'   \code{p_val}, \code{avg_log2FC}, \code{pct.1}, \code{pct.2}, \code{p_val_adj}, \code{cluster}, \code{gene}.
#'   Rows should be the marker genes (from \code{FindAllMarkers}). Here, \code{cluster}
#'   refers to the subcluster label you want annotated.
#' @param tissue_name Character. **Highly recommended**: provide clear biological context
#'   (e.g., \emph{"human aging heart"}, \emph{"mouse postnatal heart day 7"}, \emph{"human lung"}).
#' @param model Character. GPT model to use (default: 'gpt-5').
#' @param topgenenumber Integer. Number of top marker genes per subcluster to include (default: 10).
#' @param restrict_to Character vector of allowed cell types to constrain GPT output (optional).
#'
#' @return A named character vector of predicted cell types for each subcluster.
#' @keywords internal deprecated
#' @export
gptcelltype_sub <- function(input,
                            tissue_name = NULL,
                            model = 'gpt-5',
                            topgenenumber = 10,
                            restrict_to = NULL) {
  if (is.list(input) && !is.data.frame(input)) {
    collapsed <- sapply(input, paste, collapse = ',')
  } else {
    # Validate input columns.
    req_cols <- c("cluster", "gene", "avg_log2FC")
    if (!all(req_cols %in% colnames(input))) {
      stop(sprintf("`input` must be a data.frame with columns: %s",
                   paste(req_cols, collapse = ", ")))
    }

    # Build per-subcluster marker strings using genes with avg_log2FC > 1.
    df <- input[input$avg_log2FC > 1, c("cluster", "gene", "avg_log2FC"), drop = FALSE]
    if (nrow(df) == 0) stop("No rows with avg_log2FC > 1 after filtering.")

    df <- df[order(df$cluster, -df$avg_log2FC), ]
    collapsed <- tapply(seq_len(nrow(df)), df$cluster, function(idx) {
      genes <- df$gene[idx]
      paste0(utils::head(genes, topgenenumber), collapse = ",")
    })
  }

  # Build prompt strings.
  base_prompt <- paste0(
    "Identify cell types of ", tissue_name, " using the following markers separately for each\n",
    "row. Only provide the cell type name. Do not show numbers before the name.\n",
    "Some can be a mixture of multiple cell types.\n"
  )
  if (!is.null(restrict_to)) {
    base_prompt <- paste0(
      base_prompt,
      "Restrict your predictions to the following cell types:\n",
      paste(restrict_to, collapse = ", "), "\n"
    )
  }
  debug_prompt <- paste0(
    base_prompt,
    paste0(names(collapsed), ": ", unlist(collapsed), collapse = "\n")
  )

  # Check API key and fail fast with the prompt for debugging.
  OPENAI_API_KEY <- Sys.getenv("OPENAI_API_KEY")
  if (OPENAI_API_KEY == "") {
    message("OPENAI_API_KEY missing. Prompt that would have been sent:\n", debug_prompt)
    stop("Error: OpenAI API key not found. Please set OPENAI_API_KEY.")
  }

  # Batch into chunks of 30 lines.
  cutnum <- ceiling(length(collapsed) / 30)
  cid <- if (cutnum > 1) as.numeric(cut(seq_along(collapsed), cutnum)) else rep(1, length(collapsed))

  allres <- sapply(seq_len(cutnum), function(i) {
    id <- which(cid == i)
    prompt <- paste0(base_prompt, paste(unlist(collapsed[id]), collapse = "\n"))

    res <- rep("unknown", length(id))
    attempt <- 1
    success <- FALSE
    while (attempt <= 3 && !success) {
      tryCatch({
        k <- openai::create_chat_completion(
          model = model,
          message = list(list("role" = "user", "content" = prompt))
        )
        res_tmp <- strsplit(k$choices[, "message.content"], "\n")[[1]]
        if (length(res_tmp) == length(id)) {
          res <- res_tmp
          success <- TRUE
        } else {
          message(sprintf("Response lines (%d) != expected (%d) on attempt %d.",
                          length(res_tmp), length(id), attempt))
          Sys.sleep(1)
        }
      }, error = function(e) {
        message(sprintf("API call failed on attempt %d: %s", attempt, e$message))
        Sys.sleep(1)
      })
      attempt <- attempt + 1
    }

    names(res) <- names(collapsed)[id]
    res
  }, simplify = FALSE)

  return(gsub(",$", "", unlist(allres)))
}


# prev: run_annotation_on_subclusters
# For all the remaining high-level workflow and integration functions, here's the style:

#' Run Ontology-based Subcluster Annotation Workflow
#'
#' **DEPRECATED:** Use \code{\link{run_subcluster_annotation_workflow}} with \code{strategy = "ontology"} instead.
#'
#' For each subcluster (directory in base_dir), runs GPT annotation restricting output to ontology descendants.
#'
#' @param base_dir Path to directory containing subclustered Seurat objects.
#' @param cl Cell ontology object.
#' @param mapping_dict Mapping dictionary for annotation cleaning.
#' @param model GPT model (default: 'gpt-5').
#' @param tissue_name Character. Tissue context.
#' @param n_runs Number of GPT calls (default: 2).
#' @param resolutions Numeric vector of subcluster resolutions.
#' @param ontology_graph igraph object for ontology.
#' @param subcluster_prefix Metadata column prefix for subcluster IDs (default: 'subcluster_res.').
#'
#' @return Named list of annotation results per subcluster.
#' @keywords internal deprecated
#' @export
anno_subcluster_ontology <- function(base_dir = "output/subclusters",
                                     cl,
                                     mapping_dict,
                                     model = 'gpt-5',
                                     tissue_name = NULL,
                                     n_runs = 2,
                                     resolutions = c(0.1, 0.3, 0.5),
                                     ontology_graph = NULL,
                                     subcluster_prefix = "subcluster_res.") {
  subtypes <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  all_results <- list()
  for (ct_dir in subtypes) {
    ct_name <- basename(ct_dir)
    message("\nProcessing annotations for: ", ct_name)
    seurat_path <- file.path(ct_dir, "seurat_subset.rds")
    if (!file.exists(seurat_path)) {
      warning("Subset Seurat object missing for: ", ct_name)
      next
    }
    seurat_obj <- readRDS(seurat_path)
    tissue_name <- paste(gsub("_", " ", ct_name) , "in", tissue_name )
    children_names <- get_all_descendant_names(cl, ct_name)
    annotation_results <- gptanno_sub(
      seurat_obj = seurat_obj,
      resolutions = resolutions,
      cl = cl,
      mapping_dict = mapping_dict,
      marker_path = ct_dir,
      model = model,
      tissue_name = tissue_name,
      n_runs = n_runs,
      save_dir = file.path(ct_dir, "annotation"),
      save_plots = TRUE,
      save_objects = FALSE,
      ontology_graph = ontology_graph,
      subcluster_prefix = subcluster_prefix,
      children_names = children_names
    )
    all_results[[ct_name]] <- annotation_results
  }
  return(all_results)
}

# previous name: annotate_subclusters_with_inherited_parent_markers
#' Run Marker-Inheritance Subcluster Annotation Workflow
#'
#' **DEPRECATED:** Use \code{\link{run_subcluster_annotation_workflow}} with \code{strategy = "marker_inheritance"} instead.
#'
#' For each subcluster, combines top marker genes from parent and subcluster, runs GPT annotation, and returns results.
#'
#' @param sub_seurat Subcluster Seurat object.
#' @param parent_marker_file Path to parent cluster marker file.
#' @param subcluster_res Metadata column for subcluster assignment.
#' @param original_cluster_col Metadata column for parent cluster assignment.
#' @param celltype_col Metadata column for parent cell type annotation.
#' @param top_n Number of top markers per source (default: 10).
#' @param model Character. GPT model (default: 'gpt-5').
#' @param tissue_name Character. Tissue context.
#' @param n_runs Number of GPT calls.
#' @param marker_dir Directory with subcluster markers.
#' @param ontology_graph igraph object (optional).
#' @param ontology ontologyIndex object (optional).
#' @param save_dir Directory for outputs (optional).
#' @param save_plots Whether to save plots.
#' @param save_objects Whether to save annotated objects.
#'
#' @return List with `summary` (annotation summary) and `seurat` (annotated Seurat).
#' @keywords internal deprecated
#' @export
anno_subcluster_inherit <- function(
    sub_seurat,
    parent_marker_file,
    subcluster_res = "subcluster_res.0.1",
    original_cluster_col = "cluster_res.0.1",
    celltype_col = "celltype_res.0.1",
    top_n = 10,
    model = "gpt-5",
    tissue_name = "Your tissue name",
    n_runs = 1,
    marker_dir = "output/subclusters",
    ontology_graph = NULL,
    ontology = NULL,
    save_dir = NULL,
    save_plots = TRUE,
    save_objects = FALSE
) {
  res_value <- gsub("subcluster_res\\.", "", subcluster_res)
  marker_filename <- paste0("markers_res_", res_value, ".rds")
  sub_markers_path <- file.path(marker_dir, marker_filename)
  if (!file.exists(sub_markers_path)) stop("Missing subcluster marker file: ", sub_markers_path)
  sub_markers <- readRDS(sub_markers_path)
  if (!file.exists(parent_marker_file)) stop("Missing parent marker file: ", parent_marker_file)
  parent_markers <- readRDS(parent_marker_file)
  Seurat::Idents(sub_seurat) <- subcluster_res
  sub_ids <- unique(Seurat::Idents(sub_seurat))
  base_celltype <- gsub(" \\(not in GO\\)", "", unique(sub_seurat[[celltype_col]][,1]))
  gpt_input <- list()
  for (sub_id in sub_ids) {
    sub_cells <- Seurat::WhichCells(sub_seurat, idents = sub_id)
    parent_cluster_mode <- sub_seurat@meta.data[sub_cells, original_cluster_col] |>
      as.character() |> as.numeric() |> stats::na.omit() |> as.integer() |>
      as.data.frame() |> stats::setNames("cluster") |>
      dplyr::count(cluster) |> dplyr::arrange(desc(n)) |>
      dplyr::slice(1) |> dplyr::pull(cluster)
    parent_top_genes <- parent_markers |>
      dplyr::filter(cluster == parent_cluster_mode, avg_log2FC > 1) |>
      dplyr::arrange(desc(avg_log2FC)) |>
      dplyr::slice(1:top_n) |>
      dplyr::pull(gene)
    sub_top_genes <- sub_markers |>
      dplyr::filter(cluster == sub_id, avg_log2FC > 1) |>
      dplyr::arrange(desc(avg_log2FC)) |>
      dplyr::slice(1:top_n) |>
      dplyr::pull(gene)
    combined_genes <- unique(c(parent_top_genes, sub_top_genes))[1:min(20, length(c(parent_top_genes, sub_top_genes)))]
    gpt_input[[as.character(sub_id)]] <- paste(combined_genes, collapse = ", ")
  }
  annotation_summary <- summarize_gptcelltype_sub(
    markers = gpt_input,
    model = model,
    tissue_name = paste0(base_celltype, " in ", tissue_name),
    n_runs = n_runs
  )
  if (!is.null(ontology_graph)) {
    annotation_summary <- calculate_ontology_distance(
      annotation_summary,
      ontology_graph = ontology_graph,
      cl_term_map =  OntoAnno::cl_term_map
    )
  }
  annotation_summary$cluster <- as.character(annotation_summary$cluster)
  Seurat::Idents(sub_seurat) <- as.character(Seurat::Idents(sub_seurat))
  annotated_seurat <- assign_celltype(
    sub_seurat,
    annotation_summary,
    new_celltype = paste0("annotated_sub_", res_value)
  )
  if (!is.null(save_dir)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    if (save_plots) {
      grDevices::pdf(file.path(save_dir, paste0("annotation_plot_res_", res_value, ".pdf")), width = 18, height = 10)
      print(plot_celltype_comparison(
        annotated_seurat,
        original_col = subcluster_res,
        annotation_col = paste0("annotated_sub_", res_value)
      ))
      grDevices::dev.off()
    }
    if (save_objects) {
      saveRDS(annotated_seurat, file.path(save_dir, paste0("annotated_seurat_res_", res_value, ".rds")))
    }
    saveRDS(annotation_summary, file.path(save_dir, paste0("annotation_summary_res_", res_value, ".rds")))
  }
  return(list(
    summary = annotation_summary,
    seurat = annotated_seurat
  ))
}

# previous name: run_annotation_on_subclusters_v2
#' Subcluster Annotation Master Workflow (Ontology or Inheritance)
#'
#' **DEPRECATED:** Use \code{\link{run_subcluster_annotation_workflow}} instead.
#'
#' Main function to run subcluster annotation using either ontology-based or marker-inheritance strategy.
#'
#' @param base_dir Directory containing subcluster Seurat objects.
#' @param cl Ontology object.
#' @param mapping_dict Dictionary for annotation cleaning.
#' @param model GPT model.
#' @param tissue_name Context string.
#' @param n_runs Number of GPT runs.
#' @param resolutions Numeric vector of resolutions.
#' @param ontology_graph igraph object.
#' @param subcluster_prefix Metadata column prefix.
#' @param strategy "ontology" or "marker_inheritance".
#' @param parent_marker_root Directory for parent marker files (marker_inheritance).
#' @param parent_res_val Parent cluster resolution value (marker_inheritance).
#' @param parent_celltype_col Column for parent cell type.
#'
#' @return List of annotation results.
#' @keywords internal deprecated
#' @export
anno_subcluster <- function( base_dir = "output/subclusters", cl, mapping_dict, model = "gpt-5",
                             tissue_name = NULL,
                             n_runs = 2,
                             resolutions = c(0.1, 0.3, 0.5),
                             ontology_graph        = NULL,
                             subcluster_prefix     = "subcluster_res.",
                             strategy              = c("ontology", "marker_inheritance"),
                             parent_marker_root    = NULL,
                             parent_res_val        = "0.3",
                             parent_celltype_col   = "celltype_parent"
) {
  strategy  <- match.arg(strategy)
  subtypes  <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  all_out   <- list()
  for (ct_dir in subtypes) {
    ct_name  <- basename(ct_dir)
    message("\nProcessing ", ct_name, " (", strategy, ")")
    seurat_path <- file.path(ct_dir, "seurat_subset.rds")
    if (!file.exists(seurat_path)) {
      warning("No Seurat subset in ", ct_dir)
      next
    }
    sub_obj <- readRDS(seurat_path)
    tissue_string <- paste(gsub("_", " ", ct_name), "in", tissue_name)
    if (strategy == "ontology") {
      children <- get_all_descendant_names(cl, ct_name)
      annot <- gptanno_sub(
        seurat_obj        = sub_obj,
        resolutions       = resolutions,
        cl                = cl,
        mapping_dict      = mapping_dict,
        marker_path       = ct_dir,
        model             = model,
        tissue_name        = tissue_string,
        n_runs            = n_runs,
        save_dir          = file.path(ct_dir, "annotation_ontology"),
        save_plots        = TRUE,
        save_objects      = FALSE,
        ontology_graph    = ontology_graph,
        subcluster_prefix = subcluster_prefix,
        children_names    = children
      )
      all_out[[ct_name]] <- annot
      next
    }
    original_cluster_col <- paste0("cluster_res.", parent_res_val)
    parent_marker_file <- file.path(parent_marker_root, paste0("markers_res_", parent_res_val, ".rds"))
    message(paste("Access Parent marker in", parent_marker_file))
    annot_list <- list()
    for (res in resolutions) {
      subcluster_res <- paste0(subcluster_prefix, format(res, nsmall = 1))
      annot_save_dir <- file.path(ct_dir, paste0("annotation_marker_inheritance_res_", res))
      annot <- anno_subcluster_inherit(
        sub_seurat           = sub_obj,
        parent_marker_file   = parent_marker_file,
        subcluster_res       = subcluster_res,
        original_cluster_col = original_cluster_col,
        celltype_col         = parent_celltype_col,
        top_n                = 10,
        model                = model,
        tissue_name           = tissue_string,
        n_runs               = n_runs,
        ontology_graph       = ontology_graph,
        ontology             = cl,
        save_dir             = annot_save_dir,
        save_plots           = TRUE,
        save_objects         = FALSE,
        marker_dir           = ct_dir
      )
      annot_list[[paste0("res_", res)]] <- annot
    }
    score_tbl <- summarize_annotation_scores_subclusters(list(tmp = annot_list))
    best_row  <- score_tbl |>
      dplyr::slice_max(composite_score, n = 1, with_ties = FALSE)
    best_res  <- best_row$resolution
    message("   best subcluster resolution: ", best_res,
            " (score ", round(best_row$composite_score, 3), ")")
    all_out[[ct_name]] <- annot_list[best_res]
  }
  return(all_out)
}

#' GPT-based Annotation for All Subclusters at Multiple Resolutions
#'
#' **DEPRECATED:** Use \code{\link{run_subcluster_annotation_workflow}} instead.
#'
#' Calls GPT-based annotation for each subcluster and resolution.
#'
#' @param seurat_obj Seurat object.
#' @param resolutions Numeric vector.
#' @param cl Ontology.
#' @param mapping_dict Mapping dictionary.
#' @param marker_path Marker file directory.
#' @param model GPT model.
#' @param tissue_name Tissue context.
#' @param n_runs Number of runs.
#' @param save_dir Where to save outputs.
#' @param save_plots Save plot PDFs?
#' @param save_objects Save annotated objects?
#' @param ontology_graph igraph object.
#' @param subcluster_prefix Prefix for subcluster metadata column.
#' @param children_names Optional restriction of predictions.
#'
#' @return List of annotation summaries and annotated Seurat objects.
#' @keywords internal deprecated
#' @export
gptanno_sub <- function(seurat_obj, resolutions, cl,  mapping_dict,
                        marker_path = "output/marker_genes", model = 'gpt-5',
                        tissue_name = NULL, n_runs = 2,
                        save_dir = NULL, save_plots = TRUE,
                        save_objects = FALSE, ontology_graph = NULL,
                        subcluster_prefix = "subcluster_res.", children_names = NULL) {
  results_list <- list()
  for (res in resolutions) {
    message("\nRunning annotation for resolution: ", res)
    col_name <- paste0(subcluster_prefix, format(res, nsmall = 1))
    if (!col_name %in% colnames(seurat_obj@meta.data)) {
      warning("Column ", col_name, " not found in metadata. Skipping.")
      next
    } else {
      message("Setting Idents to ", col_name)
      Seurat::Idents(seurat_obj) <- col_name
    }
    marker_file <- file.path(marker_path, paste0("markers_res_", res, ".rds"))
    if (!file.exists(marker_file)) {
      warning("Marker file not found for resolution ", res, ". Skipping.")
      next
    }
    markers <- readRDS(marker_file)
    annotation_summary <- summarize_gptcelltype_sub(markers, model = model, tissue_name = tissue_name, n_runs = n_runs, restrict_to = children_names)
    if (!is.null(ontology_graph)) {
      annotation_summary <- calculate_ontology_distance(annotation_summary, ontology_graph = ontology_graph, cl_term_map)
    }
    annotation_summary$cluster <- as.character(annotation_summary$cluster)
    Seurat::Idents(seurat_obj) <- as.character(Seurat::Idents(seurat_obj))
    annotated_seurat <- assign_celltype(seurat_obj, annotation_summary, new_celltype = "annotated_sub_celltype")
    if (!is.null(save_dir)) {
      dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
      if (save_plots) {
        grDevices::pdf(file.path(save_dir, paste0("annotation_plot_res_", res, ".pdf")), width = 18, height = 10)
        print(plot_celltype_comparison(annotated_seurat, original_col = paste0("subcluster_res.", res), annotation_col = "annotated_sub_celltype"))
        grDevices::dev.off()
      }
      if (save_objects) {
        saveRDS(annotated_seurat, file.path(save_dir, paste0("annotated_seurat_res_", res, ".rds")))
      }
    }
    results_list[[paste0("res_", res)]] <- list(
      summary = annotation_summary,
      seurat = annotated_seurat
    )
  }
  return(results_list)
}

## need to check in new seaneario
#' Summarize Annotation Scores for Subclusters Across Resolutions
#'
#' **DEPRECATED:** Scores are now automatically calculated in \code{\link{run_subcluster_annotation_workflow}}.
#'
#' Computes composite annotation quality scores across subcluster resolutions for each parent cell type.
#'
#' @param annotation_sub Nested list as returned by subcluster annotation workflows.
#'
#' @return Data.frame with celltype, resolution, sum_path_length, avg_max_percentage, and composite_score.
#' @keywords internal deprecated
#' @export
summarize_annotation_scores_subclusters <- function(annotation_sub) {
  all_scores <- list()
  for (celltype in names(annotation_sub)) {
    celltype_entry <- annotation_sub[[celltype]]
    res_names <- names(celltype_entry)
    sum_path_length <- numeric(0)
    avg_max_perc <- numeric(0)
    valid_res <- character(0)
    for (res in res_names) {
      res_obj <- celltype_entry[[res]]
      fs <- tryCatch(res_obj$summary$final_summary, error = function(e) NULL)
      if (is.null(fs) || !is.data.frame(fs) || nrow(fs) == 0) next

      fs <- .filter_final_summary_for_scoring(res_obj$summary, fs)
      if (nrow(fs) == 0) next
      if (!"max_percentage" %in% names(fs) || all(is.na(fs$max_percentage))) next

      sum_path_length[res] <- if ("avg_distance" %in% names(fs)) {
        sum(fs$avg_distance, na.rm = TRUE)
      } else {
        NA_real_
      }
      avg_max_perc[res] <- mean(fs$max_percentage, na.rm = TRUE)
      valid_res <- c(valid_res, res)
    }
    if (length(valid_res) == 0) {
      dummy_res <- if (length(res_names) > 0) res_names[1] else NA_character_
      score_df <- data.frame(
        celltype = celltype,
        resolution = dummy_res,
        sum_path_length = NA_real_,
        avg_max_percentage = NA_real_,
        composite_score = NA_real_,
        stringsAsFactors = FALSE
      )
      all_scores[[celltype]] <- score_df
      next
    }
    # --- New: robust normalization for path (smaller is better), percent (bigger is better) ---
    if (all(is.na(sum_path_length))) {
      path_score <- rep(NA_real_, length(sum_path_length))
    } else {
      max_possible_path <- max(sum_path_length, na.rm = TRUE)
      # If max_possible_path == 0 (ideal), set all path_scores = 1
      if (max_possible_path == 0) {
        path_score <- rep(1, length(sum_path_length))
      } else {
        path_score <- 1 - (sum_path_length / max_possible_path)
      }
    }
    percent_score <- avg_max_perc / 100
    composite <- ifelse(
      is.na(path_score) | is.na(percent_score),
      NA_real_,
      (path_score + percent_score) / 2
    )
    # If both ideal, force composite to 1
    is_ideal <- (sum_path_length == 0) & (avg_max_perc == 100)
    if (any(is_ideal)) composite[is_ideal] <- 1

    score_df <- data.frame(
      celltype = rep(celltype, length(valid_res)),
      resolution = valid_res,
      sum_path_length = sum_path_length,
      avg_max_percentage = avg_max_perc,
      composite_score = composite,
      stringsAsFactors = FALSE
    )
    all_scores[[celltype]] <- score_df
  }
  if (length(all_scores) == 0) return(data.frame())
  final_table <- do.call(rbind, all_scores)
  final_table <- dplyr::arrange(final_table, celltype, dplyr::desc(composite_score))
  return(final_table)
}

#' Assign Best Subcluster Annotations to Full Seurat Object
#'
#' **DEPRECATED:** Use \code{\link{assign_subcluster_annotations_to_full}} instead.
#'
#' Updates a full Seurat object with the best subcluster annotation for each parent cell type.
#'
#' @param seurat_obj Full Seurat object.
#' @param summary_scores_sub Data.frame as returned by summarize_annotation_scores_subclusters().
#' @param annotation_sub Nested list of subcluster annotation objects.
#' @param parent_column Metadata column with parent annotation.
#' @param subcluster_prefix Prefix for subcluster metadata columns.
#' @param final_colname Name for new column with best annotation.
#'
#' @return The updated Seurat object with subcluster annotation column.
#' @keywords internal deprecated
#' @export
assign_best_subcluster_annotations <- function(seurat_obj,
                                               summary_scores_sub,
                                               annotation_sub,
                                               parent_column = "celltype_parent",
                                               subcluster_prefix = "subcluster_res.",
                                               final_colname = "celltype_final") {
  seurat_obj@meta.data[[final_colname]] <- seurat_obj@meta.data[[parent_column]]
  best_res_df <- summary_scores_sub |>
    dplyr::group_by(celltype) |>
    dplyr::slice_max(composite_score, n = 1, with_ties = FALSE) |>
    dplyr::ungroup()
  for (i in seq_len(nrow(best_res_df))) {
    celltype <- best_res_df$celltype[i]
    best_res <- best_res_df$resolution[i]
    best_obj <- annotation_sub[[celltype]][[best_res]]
    if (is.null(best_obj)) {
      warning("Skipping ", celltype, " no annotation found for ", best_res)
      next
    }
    annotated_seurat <- best_obj$seurat
    sub_annot_col <- grep("^annotated_sub_", colnames(annotated_seurat@meta.data), value = TRUE)
    if (length(sub_annot_col) != 1) {
      warning("Skipping ", celltype, " could not uniquely identify subcluster annotation column.")
      next
    }
    annotations <- annotated_seurat@meta.data[[sub_annot_col]]
    cell_names <- colnames(annotated_seurat)
    seurat_obj@meta.data[cell_names, final_colname] <- annotations
  }
  return(seurat_obj)
}

#' Assign Inherited Subcluster Annotations
#'
#' **DEPRECATED:** Use \code{\link{assign_subcluster_annotations_to_full}} instead.
#'
#' Combines subcluster and parent annotations to create a new celltype_subcluster column.
#'
#' @param full_seurat Seurat object for full dataset.
#' @param annotation_sub_inherited Nested list of inherited subcluster annotation objects.
#' @param subcluster_annotation_col_prefix Prefix for subcluster annotation columns.
#' @param final_colname Name of new metadata column.
#' @param parent_colname Column for parent annotation.
#' @param clean_labels Whether to clean annotation labels of " (not in GO)" suffix.
#'
#' @return The updated Seurat object with inherited subcluster annotation column.
#' @keywords internal deprecated
#' @export
assign_inherited_subcluster_annotations <- function(
    full_seurat,
    annotation_sub_inherited,
    subcluster_annotation_col_prefix = "annotated_sub_",
    final_colname = "celltype_subcluster",
    parent_colname = "celltype_parent",
    clean_labels = TRUE
) {
  annotations <- rep(NA, ncol(full_seurat))
  names(annotations) <- colnames(full_seurat)
  for (ct in names(annotation_sub_inherited)) {
    entry <- annotation_sub_inherited[[ct]]
    if (!is.list(entry) || length(entry) == 0) {
      warning("Skipping ", ct, ": invalid annotation entry.")
      next
    }
    best_res <- names(entry)[[1]]
    res_data <- entry[[best_res]]
    if (is.null(res_data$seurat)) {
      warning("Skipping ", ct, ": missing Seurat object in ", best_res)
      next
    }
    seurat_obj <- res_data$seurat
    annot_col <- grep(paste0("^", subcluster_annotation_col_prefix), colnames(seurat_obj@meta.data), value = TRUE)
    if (length(annot_col) != 1) {
      warning("Skipping ", ct, ": expected one annotation column, found ", length(annot_col), " in ", best_res)
      next
    }
    cells <- colnames(seurat_obj)
    annots <- as.character(seurat_obj[[annot_col]][, 1])
    if (clean_labels) {
      annots <- gsub(" \\(not in GO\\)", "", annots)
    }
    annotations[cells] <- annots
  }
  if (!parent_colname %in% colnames(full_seurat@meta.data)) {
    stop("`parent_colname` (", parent_colname, ") not found in Seurat metadata.")
  }
  parent_annots <- as.character(full_seurat[[parent_colname]][, 1])
  if (clean_labels) {
    parent_annots <- gsub(" \\(not in GO\\)", "", parent_annots)
  }
  missing_idx <- which(is.na(annotations))
  annotations[missing_idx] <- parent_annots[missing_idx]
  full_seurat[[final_colname]] <- annotations
  return(full_seurat)
}
