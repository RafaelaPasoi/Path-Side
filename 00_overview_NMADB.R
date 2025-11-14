sink("overview_log.txt", split = TRUE)


library(nmadb)
library(dplyr)

# 1. Φέρε τον πλήρη κατάλογο
nmas <- getNMADB()
all_ids <- nmas$Record.ID

# --- Μπορείς να το περιορίσεις αν θες
# all_ids <- head(all_ids, 200)   # π.χ. για δοκιμή με 200 μόνο

# 2. Συνάρτηση που επιστρέφει type + effect + status
get_type_effect <- function(bid) {
  out <- tryCatch({
    netb <- readByID(bid)
    tibble(
      id     = bid,
      type   = netb$type,
      effect = netb$effect,
      status = "OK"
    )
  }, error = function(e) {
    tibble(
      id     = bid,
      type   = NA,
      effect = NA,
      status = "FAIL"
    )
  })
  return(out)
}

# 3. Εφάρμοσε σε όλα τα IDs
overview_types <- bind_rows(lapply(all_ids, get_type_effect))

# 4. Γενική σύνοψη
cat("Σύνολο datasets που προσπάθησα:", length(all_ids), "\n")
cat("Πόσα διαβάστηκαν σωστά:", sum(overview_types$status == "OK"), "\n")
cat("Πόσα απέτυχαν:", sum(overview_types$status == "FAIL"), "\n\n")

cat("Κατανομή τύπων (μόνο για OK):\n")
print(table(overview_types$type, useNA = "ifany"))

cat("\nΚατανομή μέτρων αποτελέσματος (μόνο για OK):\n")
print(table(overview_types$effect, useNA = "ifany"))

# 5. Πιο ωραία με dplyr
overview_summary <- overview_types %>%
  filter(status == "OK") %>%
  count(type, effect, sort = TRUE)

overview_summary

# 6. (Προαιρετικά) δες και μερικά IDs που απέτυχαν
head(overview_types %>% filter(status == "FAIL"))

sink()
