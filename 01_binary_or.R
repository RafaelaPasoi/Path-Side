
library(nmadb)
library(dplyr)
library(netmeta)
library(writexl)

# 1. Get all IDs from NMADB ---------------------------------------------------

nmas <- getNMADB()
all_ids <- nmas$Record.ID

# 2. Function to create a clean binary OR dataset for a single ID -------------

process_binary <- function(bid, outdir = "binary_datasets") {
  
  # Create output folder if it does not exist
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  
  # Progress message (this produced the "getting dataset: xxxx" output)
  print(c("getting dataset:", bid))
  
  out <- tryCatch({
    netb <- readByID(bid)
    
    # Keep only binary datasets
    if (netb$type != "binary") {
      return(NULL)
    }
    
    dat <- netb$data
    
    # Pairwise transformation to Odds Ratios (OR)
    pw <- pairwise(
      treat   = t,
      event   = r,
      n       = n,
      studlab = id,
      data    = dat,
      sm      = "OR"   # enforce Odds Ratio
    )
    
    # Netmeta object (for later side-splitting / path-based analyses)
    nm <- netmeta(pw, common = FALSE, random = TRUE)
    
    # Prepare a clean dataset
    clean_dat <- pw %>%
      dplyr::select(studlab, treat1, treat2, event1, n1, event2, n2, TE, seTE) %>%
      dplyr::mutate(dataset_id = bid)
    
    # Save to Excel
    outfile <- file.path(outdir, paste0("Binary_OR_", bid, ".xlsx"))
    write_xlsx(clean_dat, outfile)
    
    message("Binary dataset saved: ", outfile)
    
    # Return summary info for this dataset
    list(
      id       = bid,
      nstudies = netb$nstudies,
      ntreat   = netb$ntreat,
      file     = outfile
    )
    
  }, error = function(e) {
    message("Error in ", bid, ": ", e$message)
    NULL
  })
  
  out
}

# 3. Apply the function to all IDs --------------------------------------------

binary_info <- lapply(all_ids, process_binary)

# 4. Combine information for all successfully processed datasets --------------

binary_summary <- bind_rows(binary_info)

# 5. Save the summary table ---------------------------------------------------

if (!dir.exists("binary_datasets")) {
  dir.create("binary_datasets")
}
write_xlsx(binary_summary, "binary_datasets/Binary_summary.xlsx")

# Show the summary in the knitted document
binary_summary

