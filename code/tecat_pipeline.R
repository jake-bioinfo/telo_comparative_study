#!/usr/bin/env Rscript

# Load the required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2) 
  library(cowplot)
  library(parallel)
  library(optparse)})
library(TECAT)


# Import options
option_list <- list(make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
                                help = "Print extra output [default]"),
                    make_option(c("-t", "--threads"), type = "integer", default = NULL,
                                help = "Number of threads to use"),
                    make_option(c("-o", "--output"), type = "character", default = file.path(getwd(), "output"),
                                help = "Output directory"),
                    make_option(c("-m", "--temporary"), type = "character", default = file.path(tempdir()),
                                help = "Temporary directory"),
                    make_option(c("-p", "--sample_name"), type = "character", default = "output",
                                help = "Output prefix"),
                    make_option(c("-i", "--input"), type = "character", default = NULL,
                                help = "Input file"),
                    make_option(c("-l", "--length"), type = "integer", default = 6,
                                help = "Length of telomere motif"),
                    make_option(c("-n", "--number"), type = "integer", default = 10,
                                help = "Number of repeats"),
                    make_option(c("-s", "--start"), type = "character", default = "Auto",
                                help = "Start threshold"),
                    make_option(c("-e", "--end"), type = "character", default = "Auto",
                                help = "End threshold"),
                    make_option(c("-d", "--downsample_ratio"), type = "numeric", default = 0.15,
                                help = "Ratio of total reads to to downsample for threshold determination"),
                    make_option(c("-r", "--results"), type = "character", default = file.path(getwd(), "output", "results.csv"),
                                help = "Results file"),
                    make_option(c("-f", "--platform"), type = "character", default = "pb",),
                    make_option(c("-x", "--reference"), type = "character", default = NULL,
                                help = "Reference file"),
                    make_option(c("-z", "--meme_bin"), type = "character", default = "~/meme/bin",
                                help = "Path to MEME bin directory"),
                    make_option(c("-T", "--test"), action = "store_true", default = FALSE,
                                help = "Run test mode")
)

# Start time
start_time <- Sys.time()
cat("\nStarting TECAT pipeline at: ", start_time, "\n")

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#opt$test <- TRUE

# Set test options
if ( opt$test ) {
  # Load paths
  opt$input <- file.path("/home/jake/tmp/comparative_study_temp/HG00732/lr/ERR4987503-5.fastq")
  opt$reference <- file.path("/home/jake/science_projects/telo/shared_data/hs/ref/chm13v2.0.fa.gz")
  opt$results <- file.path("/home/jake/tmp/comparative_study_temp/test_out/results.csv")
  opt$temporary <-file.path("~/tmp/comparative_study_temp/tmp") 
  opt$platform <- "ont"
  opt$sample_name <- "HG00732"
  opt$length <- 6
  opt$number <- 6
  opt$start <- "Auto"
  opt$end <- "Auto"
  opt$output <- file.path("/home/jake/tmp/comparative_study_temp/test_out")
}

# print some progress messages to stderr if \"quietly\" wasn't requested
if ( opt$verbose ) {
  write(paste("\nStarting to determine telomere read lengths at:",
              Sys.time(), collapse = ""), stderr())
}

# Determining number of threads
if ( is.null(opt$threads) ) {
  cat("\nNumber of parallel processors automatically set to max-2.\n")
  threads <- detectCores() - 2
} else {
  threads <- opt$threads
}

# Check output directory
if ( !dir.exists(opt$output) ) {
  dir.create(opt$output)
}

# Confirm all input options
if ( opt$verbose ) {
  write(paste("\nOptions:\n",
              "Threads: ", threads, "\n",
              "Output: ", opt$output, "\n",
              "Temporary: ", opt$temporary, "\n",
              "Sample_name: ", opt$sample_name, "\n",
              "Meme bin: ", opt$meme_bin, "\n",
              "Input: ", opt$input, "\n",
              "Reference: ", opt$reference, "\n",
              "Length: ", opt$length, "\n",
              "Number: ", opt$number, "\n",
              "Start: ", opt$start, "\n",
              "End: ", opt$end, "\n",
              "Platform: ", opt$platform, "\n",
              "Downsample ratio: ", opt$downsample_ratio, "\n",
              "Results: ", opt$results, "\n",
              sep = ""), stderr())
}

# Split fastq
if ( opt$verbose ) {
  write(paste("\nSplitting fastq file at:",
              Sys.time(), collapse = ""), stderr())
}

# Split fastq file
split_stats <- TECAT::auto_split(in_fastq = opt$input,
                  out_dir = file.path(opt$temporary, "split"),
                  threads = threads,
                  target_ram_size_gb = 1)

# Process reference file
if ( opt$verbose ) {
  write(paste("\nProcessing reference file at:",
              Sys.time(), collapse = ""), stderr())
}

# Assign meme_bin to options
options(meme_bin = opt$meme_bin)

# Parse platform
if (opt$platform == "pb") {
  plat <- "PB"
  opt$number = 6
} else {
  plat <- "ONT"
}

#plat <- ifelse(opt$platform == "pb", "PB", "ONT")

motif_analysis <- TECAT::telomere_motif(reference_file = opt$reference,
                                        reference_telo_length = 2000,
                                        telo_motif_length = opt$length,
                                        number_of_motifs = 5,
                                        number_of_repeats = opt$number,
                                        threads = threads, 
                                        platform = plat)

# Print out motif information
if ( opt$verbose ) {
  cat("\nMotif information:\n")
  cat("\nNumber of motifs: ", length(motif_analysis$motifs), "\n")
  cat("\nGrep list: ", "\n")
  print(motif_analysis$grep_list)
}

# Save motif information
saveRDS(motif_analysis, file = file.path(opt$output, "motif_analysis.rds"))

# Pattern searching
if ( opt$verbose ) {
  write(paste("\nPattern searching at:",
              Sys.time(), collapse = ""), stderr())
}

# Input files
input_files <- list.files(file.path(opt$temporary, "split"), pattern = "fastq", full.names = TRUE)

out <- file.path(opt$temporary, "telo")
dir.create(out, showWarnings = FALSE)
search_res <- TECAT::telo_search(fastq_files = input_files,
                                 grep_list = motif_analysis$grep_list,
                                 out_dir = out,
                                 threads = threads,
                                 return_telomeres = FALSE,
                                 verbose = opt$verbose,
                                 progress = TRUE)

# Print stats 
knitr::kable(head(search_res$telomere_stats))

if( opt$verbose ) {
  cat("\nNumber of telomeres: ", sum(search_res$telomere_stats[, "telomere_count"]), "\n")
}

# Save telomere stats
saveRDS(search_res, file = file.path(opt$output, "telomere_stats.rds"))

# Cut into sliding windows
if ( opt$verbose ) {
  write(paste("\nCutting into sliding windows at:",
              Sys.time(), collapse = ""), stderr())
}

# Determine environment
if ( opt$verbose ) {
  write(paste("\nDetermining environment at:",
              Sys.time(), collapse = ""), stderr())
}

# Determine environment
environ <- .Platform$OS.type
environ <- ifelse(environ == "unix", "linux", "windows")

# Cut into sliding windows
windows <- TECAT::sliding_window_parallel(telomere_file = search_res$combined_fasta,
                                          window_length = 200,
                                          step = 100, 
                                          environment = environ,
                                          threads = threads)

# Save windows
saveRDS(windows, file = file.path(opt$output, "windows.rds"))

# Calculate frequencies 
if ( opt$verbose ) {
  write(paste("\nCalculating frequencies at:",
              Sys.time(), collapse = ""), stderr())
}

# Calculate frequencies
freq_out <- file.path(opt$output, "frequencies")
dir.create(freq_out, showWarnings = FALSE)
freq <- TECAT::frequencies(windows = windows,
                           motifs = motif_analysis$motifs,
                           out_dir = freq_out,
                           parallel = TRUE,
                           threads = threads, 
                           verbose = opt$verbose, 
                           progress = TRUE,
                           environment = environ, 
                           save_files = TRUE) 

# Save frequencies
saveRDS(freq, file = file.path(opt$output, "freqs.rds"))

# Print head of frequencies
knitr::kable(head(freq[[1]]))

# Determine thresholds
if ( opt$verbose ) {
  write(paste("\nDetermining thresholds at:",
              Sys.time(), collapse = ""), stderr())
}

# Determine thresholds
thresholds_dataframe <- TECAT::determine_threshold(telomere_list = freq,
                                                   threads = threads,
                                                   sample_ratio = opt$downsample_ratio)

# Determine optimal thresholds
thresholds <- TECAT::optimal_thresholds(threshold_dataframe = thresholds_dataframe)

# Plot thresholds
plot_thresholds <- TECAT::plot_thresholds(threshold_dataframe = thresholds_dataframe,
                                          optimal_thresholds = thresholds)
                                              
# Save thresholds plot
ggsave(file = file.path(opt$output, "thresholds_plot.png"), plot = plot_thresholds, width = 12, height = 6)

# Determine telomere length
if ( opt$verbose ) {
  write(paste("\nDetermining telomere length at:",
              Sys.time(), collapse = ""), stderr())
}

telomere_lengths <- parallel::mclapply(X = freq,
                                       FUN = TECAT::findTelLength,
                                       st.thresh = (thresholds$sensitivity * 2),
                                       en.thresh = thresholds$telomere_length,
                                       mc.cores = threads)
telomere_lengths <- do.call(rbind, telomere_lengths)

# Print telomere lengths
knitr::kable(head(telomere_lengths))
cat("\nNumber of NAs: ", sum(is.na(telomere_lengths$telomere_end)), "\n")
# Remove NAs
telomere_lengths <- telomere_lengths %>% filter(!is.na(telomere_end))

# Save telomere lengths
saveRDS(telomere_lengths, file = file.path(opt$output, "telomere_lengths.rds"))

# Truncate and map reads
if ( opt$verbose ) {
  write(paste("\nTruncating and mapping reads at:",
              Sys.time(), collapse = ""), stderr())
}

# Truncate and map reads
trunc_dir <- file.path(opt$output, "truncated")
dir.create(trunc_dir, showWarnings = FALSE)

# Truncate and map reads
truncated <- TECAT::truncate_file(combined_telomere_file = search_res$combined_fasta,
                                  results_data_frame = telomere_lengths,
                                  write_file = TRUE,
                                  out_dir = trunc_dir,
                                  return = FALSE)

# Mapped directory
mapped_dir <- file.path(opt$output, "mapped")
dir.create(mapped_dir, showWarnings = FALSE)

# Map reads
mapped <- TECAT::map(fasta = truncated$truncated_files_list,
                     reference = opt$reference,
                     results_data_frame = telomere_lengths,
                     threads = threads,
                     out_dir = mapped_dir, 
                     preset_string = "map-ont",
                     verbose = opt$verbose)
mapped$results$sample <- opt$sample_name

# Print head of mapped
knitr::kable(head(mapped$results))

# Save mapped
saveRDS(mapped, file = file.path(opt$output, "mapped.rds"))

# Plot telomere lengths
if ( opt$verbose ) {
  write(paste("\nPlotting telomere lengths at:",
              Sys.time(), collapse = ""), stderr())
}

# Plot telomere lengths
plots <- TECAT::tecat_plot(mapped_output = mapped, 
                           out_dir = opt$output,
                           prefix = opt$sample_name, 
                           save_plots = TRUE,
                           return = TRUE)

# End Time
end_time <- Sys.time()
cat("\nTECAT pipeline completed at: ", end_time, "\n")
cat("\nTotal time: ", end_time - start_time, "\n")
