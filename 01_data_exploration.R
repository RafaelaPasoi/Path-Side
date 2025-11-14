
# Load required packages
suppressPackageStartupMessages({
  library(nmadb)
  library(dplyr)
  library(purrr)
  library(writexl)
})

# Main function: data exploration with log and Excel export
run_nmadb_exploration <- function(max_ids = 20,
                                  outdir   = "datasets",
                                  log_file = "data_exploration_log.txt") {
  
  # Inner function that does the actual work and prints progress
  core <- function() {
    cat("=== NMADB data exploration ===\n")
    cat("Date-time:", as.character(Sys.time()), "\n\n")
    
    # 1. Load the full NMADB catalogue
    nmas <- getNMADB()
    cat("Total NMAs in the catalogue:", nrow(nmas), "\n")
    
    # Use the first `max_ids` IDs for exploration
    all_ids  <- nmas$Record.ID
    some_ids <- head(all_ids, max_ids)
    cat("IDs selected for testing (first", max_ids, "):\n")
    print(some_ids)
    cat("\n")
    
    # 2. Test function: check if a dataset can be read and display basic info
    test_one <- function(bid) {
      message("=== Trying dataset ID: ", bid, " ===")
      out <- tryCatch({
        netb <- readByID(bid)
        
        cat("Studies:",   netb$nstudies,
            "| Treatments:", netb$ntreat,
            "| Type:",       netb$type,
            "| Effect:",     netb$effect, "\n")
        
        cat("Column names:\n")
        print(colnames(netb$data))
        cat("First rows:\n")
        print(utils::head(netb$data))
        cat("\n")
        
        TRUE
      }, error = function(e) {
        message("  -> Error: ", e$message)
        FALSE
      })
      out
    }
    
    # 3. Apply the test to the selected IDs
    ok_flags  <- purrr::map_lgl(some_ids, test_one)
    valid_ids <- some_ids[ok_flags]
    
    cat("\nSummary of ID testing:\n")
    cat("  Total IDs tested:", length(some_ids), "\n")
    cat("  Successful (OK):",  sum(ok_flags), "\n")
    cat("  Failed:",           sum(!ok_flags), "\n")
    cat("  Valid IDs:\n")
    print(valid_ids)
    cat("\n")
    
    # 4. Optional: light column normalization helper (if needed later)
    normalize_columns <- function(dat) {
      cn <- tolower(colnames(dat))
      colnames(dat) <- cn
      dat %>%
        rename(
          study_id = dplyr::any_of(c("study id", "study", "studlab", "trial", "id", "studyid")),
          outcome  = dplyr::any_of(c("outcome", "outcome:", "endpoint")),
          effect   = dplyr::any_of(c("effect", "effect size", "yi", "estimate")),
          n        = dplyr::any_of(c("n", "sample size", "total", "number of patients", "ntotal"))
        )
    }
    
    # 5. Build a compact overview for each valid dataset
    summarize_dataset <- function(bid) {
      out <- tryCatch({
        netb <- readByID(bid)
        dat  <- normalize_columns(netb$data)
        
        tibble(
          id     = bid,
          type   = netb$type,
          effect = netb$effect,
          cols   = paste(colnames(netb$data), collapse = ", ")
        )
      }, error = function(e) {
        NULL
      })
      out
    }
    
    overview_df <- dplyr::bind_rows(lapply(valid_ids, summarize_dataset))
    
    cat("Overview table (first rows):\n")
    print(utils::head(overview_df, 20))
    cat("\n")
    
    # 6. Save each valid dataset as an Excel file
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    save_dataset <- function(bid, outdir = outdir) {
      out <- tryCatch({
        netb <- readByID(bid)
        dat  <- netb$data
        outfile <- file.path(outdir, paste0("NMA_", bid, ".xlsx"))
        writexl::write_xlsx(dat, outfile)
        message("Saved dataset to: ", outfile)
        outfile
      }, error = function(e) {
        message("Error saving ID ", bid, ": ", e$message)
        NULL
      })
      out
    }
    
    exported_files <- lapply(valid_ids, save_dataset)
    
    # 7. Save the overview table as Excel
    overview_path <- file.path(outdir, "NMAs_overview.xlsx")
    writexl::write_xlsx(overview_df, overview_path)
    cat("Overview Excel saved to:\n")
    cat(overview_path, "\n\n")
    
    # Return objects to R
    list(
      nmas        = nmas,
      tested_ids  = some_ids,
      ok_flags    = ok_flags,
      valid_ids   = valid_ids,
      overview_df = overview_df,
      exported    = exported_files,
      overview_xlsx = overview_path
    )
  }
  
  # 8. Capture all printed output + messages into a log
  log_lines <- character(0)
  result    <- NULL
  
  log_lines <- withCallingHandlers(
    capture.output({
      result <- core()
    }, type = "output"),
    message = function(m) {
      # Redirect messages into the same stream so they appear in the log
      cat(conditionMessage(m), "\n")
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }
  )
  
  # Write log file
  writeLines(log_lines, con = log_file)
  
  message("Log saved to: ", normalizePath(log_file))
  message("Datasets and overview saved under: ", normalizePath(outdir))
  
  invisible(c(result, list(log_file = log_file)))
}

# ---- Run the exploration ----
res <- run_nmadb_exploration(
  max_ids  = 20,                          # change to all_ids length for full run
  outdir   = "datasets",                  # folder where Excel files will be saved
  log_file = "data_exploration_log.txt"   # log file with all printed output
)

# You can now use:
# res$valid_ids
# res$overview_df
# View(res$overview_df)
