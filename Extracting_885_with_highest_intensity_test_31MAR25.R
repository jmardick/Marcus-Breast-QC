# Load necessary package
library(MSnbase)
library(ggplot2)

# Step 1: Read the mzML file
mzml_file <- "C:/Users/Jacob/Documents/Breast-BernieMarcus/OR Data/BC_004_05DEC24_IntraOP.mzML"
file_name <- tools::file_path_sans_ext(basename(mzml_file))
raw_data <- readMSData(mzml_file, mode = "onDisk")
results_dir <- "C:/Users/Jacob/Documents/Breast-BernieMarcus/OR Data/Extracted Scans"

# Number of spectra
length(raw_data)

# Retention times (first 10)
head(rtime(raw_data), 10)

# Plot a quick TIC (Total Ion Chromatogram)
plot(rtime(raw_data), tic(raw_data), type = "l",
     xlab = "Retention Time (s)", ylab = "TIC",
     main = "Total Ion Chromatogram")


# Parameters
target_mz <- 885.55
ppm <- 10
scan_min <- 1
scan_max <- as.numeric(length(raw_data))

# Subset the scan range
scan_range <- raw_data[scan_min:scan_max]

# Function to get intensity for target m/z
get_ion_intensity <- function(spectrum, target_mz, ppm) {
  mzs <- spectrum@mz
  intensities <- spectrum@intensity
  tol <- target_mz * ppm / 1e6
  idx <- which(abs(mzs - target_mz) <= tol)
  if (length(idx) > 0) return(sum(intensities[idx]))
  return(0)
}

# Loop through spectra and collect intensity
results <- lapply(seq_along(scan_range), function(i) {
  s <- scan_range[[i]]
  data.frame(scan_num = scan_min + i - 1,
             intensity = get_ion_intensity(s, target_mz, ppm))
})

# Combine into dataframe
results_df <- do.call(rbind, results)

# Find the scan with highest intensity
top_hit <- results_df[which.max(results_df$intensity), ]

# Print result
print(top_hit)


plot(results_df$scan_num, results_df$intensity, type = "l",
     xlab = "Scan Number", ylab = "Intensity",
     main = paste("EIC for m/z", target_mz))
abline(v = top_hit$scan_num, col = "red", lty = 2)

best_spectrum <- raw_data[[top_hit$scan_num]]
plot(best_spectrum@mz, best_spectrum@intensity, type = "h", xlab = "m/z", ylab = "Intensity", main = paste("Full MS Spectrum - Scan", top_hit$scan_num))


# Extract the spectrum
best_spectrum <- raw_data[[top_hit$scan_num]]

# Plot it
plot(best_spectrum@mz, best_spectrum@intensity, type = "h",
     xlab = "m/z", ylab = "Intensity",
     main = paste("Full MS Spectrum - Scan", top_hit$scan_num))

spec_df <- data.frame(mz = best_spectrum@mz,
                      intensity = best_spectrum@intensity)

write.csv(spec_df, file = file.path(results_dir, paste0(file_name, "spectrum_scan_", top_hit$scan_num, ".csv")), row.names = FALSE)

# Get best scan spectrum
best_spectrum <- raw_data[[top_hit$scan_num]]

# Total Ion Current (TIC)
tic_val <- sum(best_spectrum@intensity)

# Target ion intensity (re-use same function)
target_mz <- 885.55
ppm <- 10
tol <- target_mz * ppm / 1e6

# Find peak index
idx <- which(abs(best_spectrum@mz - target_mz) <= tol)
ion_intensity <- sum(best_spectrum@intensity[idx])

# Normalize to TIC
normalized_intensity <- ion_intensity / tic_val

cat("Raw intensity:", ion_intensity, "\n")
cat("TIC:", tic_val, "\n")
cat("TIC-normalized intensity:", normalized_intensity, "\n")

# Define noise window (e.g., ±0.5 Da around the ion but excluding ±ppm)
mzs <- best_spectrum@mz
ints <- best_spectrum@intensity

noise_idx <- which(mzs >= (target_mz - 0.5) & mzs <= (target_mz + 0.5) &
                     (abs(mzs - target_mz) > tol))  # exclude the actual peak

noise_intensities <- ints[noise_idx]
noise_baseline <- median(noise_intensities)  # or use mean()

# SNR estimate
snr <- ion_intensity / noise_baseline

cat("Estimated SNR:", snr, "\n")

# ---- EIC Plot ----
eic_plot <- ggplot(results_df, aes(x = scan_num, y = intensity)) +
  geom_line(color = "#1f77b4", linewidth = 1.2) +
  geom_vline(xintercept = top_hit$scan_num, linetype = "dashed", color = "red", linewidth = 1.2) +
  labs(title = paste("EIC for m/z", target_mz),
       x = "Scan Number", y = "Intensity") +
  theme_minimal(base_size = 14)
ggsave(eic_plot, file = file.path(results_dir, paste0(file_name, "eic_scan.png")))

# ---- Spectrum Plot ----

# Prepare data for plotting
spectrum_df <- data.frame(mz = best_spectrum@mz,
                          intensity = best_spectrum@intensity)
# Plot with annotation
spectrum_plot <- ggplot(spectrum_df, aes(x = mz, y = intensity)) +
  geom_segment(aes(xend = mz, yend = 0), color = "black") +
  geom_vline(xintercept = target_mz, color = "blue", linetype = "dashed", linewidth = 1.2) +
  labs(title = paste("Mass spectrum of Scan", top_hit$scan_num),
       x = "m/z", y = "Intensity") +
  annotate("text", x = target_mz - 5,  # adjust x position as needed
           y = max(spectrum_df$intensity) * 1.1,
           label = paste("Target m/z:", target_mz),
           color = "blue", size = 5, hjust = 1) +
  theme_minimal(base_size = 14)
ggsave(spectrum_plot, file = file.path(results_dir, paste0(file_name, "spectrum_plot_", top_hit$scan_num, ".png")))
