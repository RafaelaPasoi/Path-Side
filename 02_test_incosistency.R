# Πακέτα
library(nmadb)
library(dplyr)
library(netmeta)
library(writexl)

# --- Διάλεξε dataset ---
bid <- 473269
netb <- readByID(bid)

# --- Βήμα 1: Pairwise (πάντα OR για binary) ---
pw <- pairwise(
  treat   = t,
  event   = r,
  n       = n,
  studlab = id,
  data    = netb$data,
  sm      = "OR"
)

# --- Βήμα 2: Network meta-analysis ---
nm <- netmeta(pw, common = FALSE, random = TRUE)

# --- Βήμα 3: Side-splitting ---
ns <- netsplit(nm)

side_df <- ns$compare.random %>%
  select(comparison, p) %>%
  rename(p_side_split = p) %>%
  left_join(
    ns$direct.random %>% select(comparison, TE_dir = TE, Q_split = Q),
    by = "comparison"
  ) %>%
  left_join(
    ns$indirect.random %>% select(comparison, TE_indir = TE),
    by = "comparison"
  ) %>%
  mutate(df_split = 1L) %>%
  select(comparison, TE_dir, TE_indir, Q_split, p_side_split, df_split)

# --- Βήμα 4: Path-based inconsistency ---
trts <- nm$trts
pairs <- t(combn(trts, 2))

path_list <- lapply(1:nrow(pairs), function(i) {
  a <- as.character(pairs[i,1])
  b <- as.character(pairs[i,2])
  out <- tryCatch({
    np <- netpath(nma_obj = nm, node1 = a, node2 = b)
    tibble(
      comparison = np$Comparison,
      Q_path     = np$Q,
      p_path     = np$p_value,
      n_paths    = np$`No. of independent paths`,
      df_path    = np$`No. of independent paths` - 1
    )
  }, error = function(e) NULL)
  return(out)
})

path_df <- bind_rows(path_list)

# --- Βήμα 5: Αποθήκευση σε Excel ---
write_xlsx(
  list(
    "side_splitting" = side_df,
    "path_based"     = path_df
  ),
  "test_inconsistency.xlsx"
)

message("✅ Έτοιμο! Δες το αρχείο test_inconsistency.xlsx")

