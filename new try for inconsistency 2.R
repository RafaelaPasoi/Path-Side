# ????????????
library(nmadb)
library(dplyr)
library(netmeta)
library(writexl)
library(igraph)
library(MASS)

# --- ?????????????? dataset ---
bid <- 473552

netb <- readByID(bid)

# --- ???????? 1: Pairwise (?????????? OR ?????? binary) ---
pw <- pairwise(
  treat   = t,
  event   = r,
  n       = n,
  studlab = id,
  data    = netb$data,
  sm      = "OR"
)

# --- ???????? 2: Network meta-analysis ---
nm <- netmeta(pw, common = FALSE, random = TRUE)

# --- ???????? 3: Side-splitting ---
ns <- netsplit(nm)

side_df <- ns$compare.random |>
  as.data.frame() |>
  dplyr::transmute(
    comparison,
    Q_split = z^2,
    p_side_split = p,
    df_split = 1L
  ) |>
  dplyr::left_join(
    dplyr::select(as.data.frame(ns$direct.random),  comparison, TE_dir  = TE),
    by = "comparison"
  ) |>
  dplyr::left_join(
    dplyr::select(as.data.frame(ns$indirect.random), comparison, TE_indir = TE),
    by = "comparison"
  ) |>
  dplyr::select(comparison, TE_dir, TE_indir, Q_split, p_side_split, df_split)

# -----------------------------
# Teacher path-based method as a function for any pair (a, b)
# -----------------------------
path_stats_onepair <- function(nm, a, b, tol = 1e-8) {
  # 1) Hat matrix (common-effect hat per teacher???s code)
  Hc <- hatmatrix(nm, method = "Davies", type = "full")$common
  row_lab <- paste0(a, ":", b)
  
  # Build directed graph from hat row (teacher function)
  build_directed_network_from_hat_row <- function(hat_matrix, row_label) {
    row_vals <- hat_matrix[row_label, , drop = FALSE]
    edge_list <- character()
    for (col_label in colnames(row_vals)) {
      nodes <- unlist(strsplit(col_label, ":"))
      if (length(nodes) == 2) {
        i <- nodes[1]; j <- nodes[2]
        val <- row_vals[1, col_label]
        if (val > 0) edge_list <- c(edge_list, i, j)
        else if (val < 0) edge_list <- c(edge_list, j, i)
      }
    }
    graph(edges = edge_list, directed = TRUE)
  }
  
  g <- build_directed_network_from_hat_row(Hc, row_lab)
  
  # If no path exists, return NA row
  if (!(a %in% V(g)$name) || !(b %in% V(g)$name)) {
    return(tibble(
      comparison = paste(a, b, sep=":"),
      Q_path = NA_real_, p_path = NA_real_,
      n_paths = 0L, df_path = NA_integer_
    ))
  }
  
  # 2) All simple paths a -> b
  all_paths <- all_simple_paths(g, from = a, to = b)
  if (length(all_paths) == 0) {
    return(tibble(
      comparison = paste(a, b, sep=":"),
      Q_path = NA_real_, p_path = NA_real_,
      n_paths = 0L, df_path = NA_integer_
    ))
  }
  path_list <- lapply(all_paths, function(p) V(g)[p]$name)
  
  # 3) Path adjacency matrix (counts of shared edges; diag = path length)
  edgeify <- function(path_nodes) {
    nodes <- as.vector(path_nodes)
    if (length(nodes) < 2) return(character(0))
    # undirected label for matching comparisons regardless of order
    mapply(function(x,y) paste(sort(c(x,y)), collapse=":"), nodes[-length(nodes)], nodes[-1])
  }
  edge_sets <- lapply(path_list, edgeify)
  nP <- length(edge_sets)
  A <- matrix(0, nP, nP)
  diag(A) <- sapply(edge_sets, length)
  for (i in seq_len(nP)) for (j in seq_len(nP)) if (i != j)
    A[i, j] <- length(intersect(edge_sets[[i]], edge_sets[[j]]))
  rownames(A) <- colnames(A) <- paste0("Path", seq_len(nP))
  
  # 4) Reduce to independent paths (QR)
  qra <- qr(A, tol = tol)
  keep_idx <- qra$pivot[seq_len(qra$rank)]
  A_red <- A[keep_idx, keep_idx, drop = FALSE]
  kept_paths <- path_list[keep_idx]
  Pprime <- length(kept_paths)
  
  # 5) Build C and V aligned to comparisons used on kept paths
  # Network comparisons (standardized "A:B" sorted)
  comp_all <- nm$comparisons
  comp_std <- sapply(strsplit(gsub(":", "-", comp_all), "-"),
                     function(x) paste(sort(x), collapse=":"))
  
  # edges actually present on kept paths
  kept_edges <- sort(unique(unlist(lapply(kept_paths, edgeify))))
  # keep only edges that are true network comparisons
  valid_edges <- kept_edges[kept_edges %in% comp_std]
  
  # C: rows = kept paths; cols = valid_edges; 1 if edge on path
  C <- matrix(0, nrow = Pprime, ncol = length(valid_edges))
  rownames(C) <- paste0("Path", seq_len(Pprime))
  colnames(C) <- valid_edges
  for (i in seq_len(Pprime)) {
    pe <- edgeify(kept_paths[[i]])
    C[i, colnames(C) %in% pe] <- 1
  }
  
  # V diagonal from (seTE.direct.common)^2
  V_diag_lookup <- (nm$seTE.direct.common)^2
  getV <- function(pair) {
    tr <- unlist(strsplit(pair, ":", fixed=TRUE))
    V_diag_lookup[tr[1], tr[2]]
  }
  V <- diag(sapply(valid_edges, getV))
  rownames(V) <- colnames(V) <- valid_edges
  
  # 6) S = C V C^T  (+ tiny ridge; symmetric)
  S <- C %*% V %*% t(C)
  if (!isSymmetric(S)) S <- (S + t(S))/2
  eps <- 1e-10 * mean(diag(S))
  S_spd <- S + diag(eps, nrow(S))
  S_inv <- tryCatch(solve(S_spd), error = function(e) MASS::ginv(S_spd))
  
  # 7) ??_path and ????_ab
  TEdir <- nm$TE.direct.common
  theta_path <- vapply(kept_paths, function(nodes) {
    nodes <- as.vector(nodes)
    if (length(nodes) < 2) return(0)
    sum(mapply(function(x,y) TEdir[x,y], nodes[-length(nodes)], nodes[-1]), na.rm = TRUE)
  }, numeric(1))
  
  theta_ab <- nm$TE.common[a, b]
  
  # 8) Q, df, p
  theta_diff <- theta_path - theta_ab
  Q <- as.numeric(t(theta_diff) %*% S_inv %*% theta_diff)
  df <- length(theta_path) - 1L
  p  <- 1 - pchisq(Q, df = df)
  
  tibble(
    comparison = paste(a, b, sep=":"),
    Q_path = Q,
    p_path = p,
    n_paths = length(theta_path),
    df_path = df
  )
}

# --- ???????? 4: Path-based inconsistency (?????????????? ????????????????, ????????????????????) ---

path_stats_onepair <- function(nm, a, b, tol = 1e-8) {
  # 1) Hat matrix (Davies, full) ?????? ?????????????? ???????????? ????????????
  Hc_full <- hatmatrix(nm, method = "Davies", type = "full")$common
  Hc_mat  <- if (is.list(Hc_full) && "common" %in% names(Hc_full)) Hc_full$common else Hc_full
  
  row_lab <- paste0(a, ":", b)
  if (!(row_lab %in% rownames(Hc_mat))) {
    # ????????????????: ???????????????????????? ?????????? "A:B" (????????????????????)
    row_lab <- paste(sort(c(a, b)), collapse = ":")
  }
  
  # Helper: ?????????????? ?????? ???? ???????????? ?????? hat matrix (???????????????? igraph API)
  build_directed_network_from_hat_row <- function(hat_matrix, row_label, all_trts = NULL) {
    row_vals <- hat_matrix[row_label, , drop = FALSE]
    edges <- list()
    for (col_label in colnames(row_vals)) {
      nodes <- strsplit(col_label, ":", fixed = TRUE)[[1]]
      if (length(nodes) != 2) next
      val <- row_vals[1, col_label]
      if (is.na(val) || val == 0) next
      if (val > 0) edges[[length(edges) + 1]] <- c(nodes[1], nodes[2]) else edges[[length(edges) + 1]] <- c(nodes[2], nodes[1])
    }
    if (length(edges) == 0) {
      g <- igraph::make_empty_graph(directed = TRUE)
      if (!is.null(all_trts)) g <- igraph::add_vertices(g, nv = length(all_trts), name = all_trts)
      return(g)
    }
    el <- do.call(rbind, edges)
    g <- igraph::graph_from_edgelist(el, directed = TRUE)
    if (!is.null(all_trts)) {
      missing <- setdiff(all_trts, igraph::V(g)$name)
      if (length(missing) > 0) g <- igraph::add_vertices(g, length(missing), name = missing)
    }
    g
  }
  
  g <- build_directed_network_from_hat_row(Hc_mat, row_lab, all_trts = nm$trts)
  
  # ???? ?????? ?????????????? ???????????????? a->b, ?????????? NA
  if (!(a %in% igraph::V(g)$name) || !(b %in% igraph::V(g)$name)) {
    return(tibble::tibble(
      comparison = paste(a, b, sep=":"),
      Q_path = NA_real_, p_path = NA_real_,
      n_paths = 0L, df_path = NA_integer_
    ))
  }
  all_paths <- igraph::all_simple_paths(g, from = a, to = b)
  if (length(all_paths) == 0) {
    return(tibble::tibble(
      comparison = paste(a, b, sep=":"),
      Q_path = NA_real_, p_path = NA_real_,
      n_paths = 0L, df_path = NA_integer_
    ))
  }
  path_list <- lapply(all_paths, function(p) igraph::V(g)[p]$name)
  
  # 2) Path-adjacency A (?????????????????? = ?????????? ????????????????????)
  edgeify <- function(path_nodes) {
    nodes <- as.vector(path_nodes)
    if (length(nodes) < 2) return(character(0))
    mapply(function(x,y) paste(sort(c(x,y)), collapse=":"), nodes[-length(nodes)], nodes[-1])
  }
  edge_sets <- lapply(path_list, edgeify)
  nP <- length(edge_sets)
  A <- matrix(0, nP, nP)
  diag(A) <- sapply(edge_sets, length)
  for (i in seq_len(nP)) for (j in seq_len(nP)) if (i != j)
    A[i, j] <- length(intersect(edge_sets[[i]], edge_sets[[j]]))
  rownames(A) <- colnames(A) <- paste0("Path", seq_len(nP))
  
  # 3) ???????????? ???? ???????????????????? ?????????????????? (QR)
  qra <- qr(A, tol = tol)
  keep_idx   <- qra$pivot[seq_len(qra$rank)]
  kept_paths <- path_list[keep_idx]
  Pprime     <- length(kept_paths)
  
  # 4) ?????????????? ???????????????????? ???????????? ?????? ???? (a:b)
  #    (i) ???????????????????? ?????? ??????????????, (ii) edges ???? kept_paths, (iii) edges ???? ???? ???????????????? hat-weight ?????? ???????????????????????? row
  comp_all <- nm$comparisons
  comp_std <- sapply(strsplit(gsub(":", "-", comp_all), "-"),
                     function(x) paste(sort(x), collapse=":"))
  kept_edges <- sort(unique(unlist(lapply(kept_paths, edgeify))))
  edge_cols_raw <- names(which(abs(Hc_mat[row_lab, ]) > 1e-12))
  edge_cols_std <- sapply(strsplit(gsub(":", "-", edge_cols_raw), "-"),
                          function(x) paste(sort(x), collapse=":"))
  valid_edges <- Reduce(intersect, list(kept_edges, comp_std, edge_cols_std))
  if (length(valid_edges) == 0L) {
    return(tibble::tibble(
      comparison = paste(a, b, sep=":"),
      Q_path = NA_real_, p_path = NA_real_,
      n_paths = Pprime, df_path = max(Pprime - 1L, 0L)
    ))
  }
  
  # 5) C (paths ?? edges) ?????? V (?????????????????? ?????? variances ?????? direct)
  C <- matrix(0, nrow = Pprime, ncol = length(valid_edges))
  rownames(C) <- paste0("Path", seq_len(Pprime))
  colnames(C) <- valid_edges
  for (i in seq_len(Pprime)) {
    pe <- edgeify(kept_paths[[i]])
    C[i, colnames(C) %in% pe] <- 1
  }
  
  V_lookup <- (nm$seTE.direct.common)^2
  getV <- function(pair) { tr <- strsplit(pair, ":", fixed = TRUE)[[1]]; V_lookup[tr[1], tr[2]] }
  V <- diag(sapply(valid_edges, getV))
  rownames(V) <- colnames(V) <- valid_edges
  
  # 6) S = C V C??? (???????????????????????????????? + ridge ?????? PD) ?????? S?????
  S <- C %*% V %*% t(C)
  if (!isSymmetric(S)) S <- (S + t(S)) / 2
  eps <- 1e-10 * mean(diag(S))
  S   <- S + diag(eps, nrow(S))
  S_inv <- tryCatch(solve(S), error = function(e) MASS::ginv(S))
  
  # 7) ??_path ?????? ????_ab
  TEdir <- nm$TE.direct.common
  theta_path <- vapply(kept_paths, function(nodes) {
    nodes <- as.vector(nodes)
    if (length(nodes) < 2) return(0)
    sum(mapply(function(x,y) TEdir[x,y], nodes[-length(nodes)], nodes[-1]), na.rm = TRUE)
  }, numeric(1))
  theta_ab <- nm$TE.common[a, b]
  
  # 8) Q, df, p
  theta_diff <- theta_path - theta_ab
  Q  <- as.numeric(t(theta_diff) %*% S_inv %*% theta_diff)
  df <- length(theta_path) - 1L
  p  <- 1 - pchisq(Q, df = df)
  
  tibble::tibble(
    comparison = paste(a, b, sep=":"),
    Q_path = Q,
    p_path = p,
    n_paths = length(theta_path),
    df_path = df
  )
}

# ?????????????? ?????? ?????? ???? ?????????? ??????????????????
trts  <- nm$trts
pairs <- t(combn(trts, 2))

path_list <- lapply(seq_len(nrow(pairs)), function(i) {
  a <- as.character(pairs[i, 1])
  b <- as.character(pairs[i, 2])
  tryCatch(
    path_stats_onepair(nm, a, b),
    error = function(e) {
      tibble::tibble(
        comparison = paste(a, b, sep=":"),
        Q_path = NA_real_, p_path = NA_real_,
        n_paths = NA_integer_, df_path = NA_integer_
      )
    }
  )
})

path_df <- dplyr::bind_rows(path_list)


# --- ???????? 5: ???????????????????? ???? Excel ---
write_xlsx(
  list(
    "side_splitting" = side_df,
    "path_based"     = path_df
  ),
  "test_inconsistency2.xlsx"
)

