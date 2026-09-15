#' Preprocess a Seurat Object for Downstream Analysis
#'
#' This function normalizes, identifies variable features, scales, and performs PCA on a Seurat object.
#' Optionally saves the processed object to disk if a save path is provided.
#'
#' @param seurat_obj A Seurat object.
#' @param assay Character. The assay to use (default: "RNA").
#' @param nfeatures Integer. Number of highly variable features to keep (default: 3000).
#' @param scale_factor Numeric. Scale factor for normalization (default: 10000).
#' @param npcs Integer. Number of principal components to compute (default: 30).
#' @param save_path Character. Path to save the processed Seurat object as RDS (default: NULL). If NULL, the object is not saved.
#'
#' @return The processed Seurat object.
#' @importFrom Seurat DefaultAssay NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP VariableFeatures
#' @export
#' @examples
#' # seurat_obj <- preprocess_seurat_object(seurat_obj, save_path = "preprocessed_seurat.rds")
preprocess_seurat_object <- function(seurat_obj,
                                     assay = "RNA",
                                     nfeatures = 3000,
                                     scale_factor = 10000,
                                     npcs = 30,
                                     save_path = NULL) {
  Seurat::DefaultAssay(seurat_obj) <- assay
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = scale_factor)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = Seurat::VariableFeatures(seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(seurat_obj), npcs = npcs)

  # Run UMAP if not already present
  if (!"umap" %in% names(seurat_obj@reductions)) {
    message("Running UMAP...")
    seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:npcs, reduction = "pca")
  } else {
    message("UMAP reduction already exists, skipping RunUMAP")
  }

  # Save the object only if save_path is provided
  if (!is.null(save_path)) {
    saveRDS(seurat_obj, file = save_path)
  }

  return(seurat_obj)
}

#' Run Multi-Resolution Clustering and Marker Detection
#'
#' This function runs graph-based clustering at multiple resolutions and finds marker genes for each, saving results to disk.
#'
#' @param seurat_obj A Seurat object.
#' @param resolutions Numeric vector of resolution parameters (e.g., c(0.1, 0.2, 0.3)).
#' @param result_dir Character. Directory to save marker RDS files (default: "results/markers").
#' @param dims Numeric vector. Which dimensions to use for clustering (default: 1:30).
#' @param assay Character. Assay to use, or NULL for current (default: NULL).
#' @param group.by Character. Metadata column to use for group assignment (default: "seurat_clusters").
#' @param reduction Character. Reduction to use for graph construction (default: "pca").
#' @param use_existing_neighbors Logical. Use an existing neighbor graph if present (default: TRUE).
#' @param algorithm Integer passed to Seurat::FindClusters: 1 Louvain (default), 2 Louvain with
#'   multilevel refinement, 3 SLM, 4 Leiden.
#' @param method Character passed to Seurat::FindClusters ("matrix" default; "igraph" recommended
#'   for Leiden on large datasets).
#' @param graph.name Character. Graph to cluster on (default: the SNN graph of the active assay,
#'   "<assay>_snn"). On multimodal objects this avoids silently clustering on a WNN or kNN graph.
#' @param ... Further arguments passed to Seurat::FindClusters (e.g. n.iter, random.seed).
#'
#' @return A list with: updated Seurat object, marker lists, and cluster assignments for each resolution.
#' @importFrom Seurat DefaultAssay Assays FindNeighbors FindClusters VariableFeatures Idents FindAllMarkers
#' @export
#' @examples
#' # result <- run_multi_resolution_clustering(seurat_obj, c(0.1, 0.2, 0.3))
run_multi_resolution_clustering <- function(seurat_obj,
                                            resolutions,
                                            result_dir = "output/markers",
                                            dims = 1:30,
                                            assay = NULL,
                                            group.by = "seurat_clusters",
                                            reduction = "pca",
                                            use_existing_neighbors = TRUE,
                                            algorithm = 1,
                                            method = "matrix",
                                            graph.name = NULL,
                                            ...) {
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }
  if (!is.null(assay)) {
    if (!(assay %in% Seurat::Assays(seurat_obj))) {
      stop(paste("Assay", assay, "not found. Available assays:", paste(Seurat::Assays(seurat_obj), collapse = ", ")))
    }
    Seurat::DefaultAssay(seurat_obj) <- assay
  }
  # Graph selection. If graph.name is given it must exist (or be created below); otherwise the
  # SNN graph of the active assay is used, never the first graph on the object (on a multimodal
  # object that may be a WNN graph, and FindNeighbors lists the kNN graph before the SNN graph).
  snn_default <- paste0(Seurat::DefaultAssay(seurat_obj), "_snn")
  if (use_existing_neighbors && length(names(seurat_obj@graphs)) > 0 &&
      (is.null(graph.name) && snn_default %in% names(seurat_obj@graphs) || !is.null(graph.name) && graph.name %in% names(seurat_obj@graphs))) {
    if (is.null(graph.name)) graph.name <- snn_default
    message("Using existing graph for clustering: ", graph.name)
  } else {
    message("Running FindNeighbors to generate neighbor graph...")
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, reduction = reduction, dims = dims, assay = assay, verbose = TRUE)
    if (is.null(graph.name)) graph.name <- snn_default
    if (!graph.name %in% names(seurat_obj@graphs)) stop("Graph '", graph.name, "' not found after FindNeighbors; available: ", paste(names(seurat_obj@graphs), collapse = ", "))
  }
  message("Clustering with algorithm = ", algorithm, " (1 Louvain, 2 Louvain refined, 3 SLM, 4 Leiden), method = ", method)
  all_markers_list <- list()
  cluster_assignments <- list()
  for (res in resolutions) {
    cat("Running clustering at resolution:", res, "\n")
    seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = res, graph.name = graph.name,
                                       algorithm = algorithm, method = method, ...)
    cluster_col <- paste0("cluster_res.", res)
    seurat_obj@meta.data[[cluster_col]] <- seurat_obj@meta.data$seurat_clusters
    cluster_assignments[[as.character(res)]] <- seurat_obj@meta.data[[cluster_col]]
    Seurat::Idents(seurat_obj) <- seurat_obj@meta.data[[cluster_col]]
    markers <- Seurat::FindAllMarkers(seurat_obj, assay = assay, only.pos = TRUE, group.by = group.by, min.pct = 0.25, logfc.threshold = 0.25)
    all_markers_list[[paste0("res_", res)]] <- markers
    saveRDS(markers, file = file.path(result_dir, paste0("markers_res_", res, ".rds")))
  }
  return(list(
    seurat_obj = seurat_obj,
    markers = all_markers_list,
    clusters = cluster_assignments
  ))
}
