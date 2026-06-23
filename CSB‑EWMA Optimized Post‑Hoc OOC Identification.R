# ======================================================
# CSB-EWMA: Post‑Hoc OOC Stream Identification (Optimized)
# Parallel Version with Exact Variance (Precomputed)
# ======================================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(doParallel)
library(foreach)
library(parallel)

set.seed(1234)

# ------------------------------------------------------
# Create organized output folders
# ------------------------------------------------------
dir.create("results", showWarnings = FALSE)
dir.create("results/diagnostics", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/examples", showWarnings = FALSE)

# ------------------------------------------------------
# Global parameters
# ------------------------------------------------------
p0 <- 0.5 # in-control Bernoulli probability
k  <- 10   # streams

distributions <- c("normal")

optimal_combinations <- data.frame(
  lambda = c(0.175, 0.15),
  L      = c(1.375, 1.525),
  ARL0_target = c(370, 500),
  Group = c("ARL0_370", "ARL0_500"),
  stringsAsFactors = FALSE
)
optimal_combinations$combo_label <- sprintf("λ=%.3f, L=%.3f",
                                            optimal_combinations$lambda,
                                            optimal_combinations$L)


# Diagnostic settings
n_sim_diag    <- 2000
max_time_diag <- 2000
alpha_diag    <- 0.05
m_values      <- c(0,1,2)
delta_diag    <- c(0, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4)

# ------------------------------------------------------
# Compute exact CSB-EWMA variance for a single t
# ------------------------------------------------------
var_rt_exact_single <- function(lambda, t) {
  if (t == 0) return(0)
  
  S <- 0 
  
  for (j in 1:(t-1)) {
    weights_j <- (1 - lambda)^(2*t - j - (j+1):t)
    cov_terms <- sqrt(j / ((j+1):t))
    S <- S + 2 * sum(weights_j * cov_terms)
  }
  
  diag_weights <- (1 - lambda)^(2*t - 2*(1:t))
  S <- S + sum(diag_weights)
  
  return(lambda^2 * S)
}

# ------------------------------------------------------
# Precompute Var(r_t) – exact up to convergence, then 1
# ------------------------------------------------------
precompute_variance <- function(lambda, max_t = max_time_diag, converge_t = 1000) {
  varvec <- numeric(max_t)
  for (t in 1:min(converge_t, max_t)) {
    varvec[t] <- var_rt_exact_single(lambda, t)
  }
  if (max_t > converge_t) {
    # After convergence, variance is essentially 1
    varvec[(converge_t+1):max_t] <- 1
  }
  return(varvec)
}

# ------------------------------------------------------
# Distribution helpers
# ------------------------------------------------------
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n, -0.5, 0.5)
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

generate_mixed_continuous <- function(distribution, k, shift_vec, p0 = 0.5) {
  out <- numeric(k)
  for (i in 1:k) {
    shift <- shift_vec[i]
    
    if (distribution == "normal") {
      threshold_in <- qnorm(p0)
      p1 <- min(max(p0 + shift, 0.01), 0.99)
      mean_shift <- threshold_in - qnorm(1 - p1)
      out[i] <- rnorm(1, mean = mean_shift, sd = 1)
    } else {
      stop("Only Normal distribution enabled at this stage.")
    }
  }
  return(out)
}

dichotomize_data <- function(data, distribution, p0 = 0.5) {
  if (distribution == "normal") th <- qnorm(p0)
  else stop("Only Normal allowed now.")
  
  as.integer(data > th)
}

# ------------------------------------------------------
# CSB‑EWMA diagnostic simulation 
# ------------------------------------------------------
simulate_one_diagnostic <- function(lambda, L, distribution, delta, m,
                                    var_cache, max_time = 600,
                                    p0 = 0.5, k = 10, alpha = 0.05) {
  
  ooc_idx <- sample(1:k, size = m, replace = FALSE)
  is_ooc <- rep(FALSE, k)
  is_ooc[ooc_idx] <- TRUE
  shift_vec <- ifelse(is_ooc, delta, 0)
  
  mu0      <- k * p0
  sigma2_0 <- k * p0 * (1 - p0)
  
  cum_sum  <- 0
  r_prev   <- 0
  r_hist   <- numeric(max_time)
  bin_hist <- matrix(0, nrow = k, ncol = max_time)
  
  # ---------------------------------------------
  # Simulation loop – use while for explicit break
  # ---------------------------------------------
  t <- 1
  signal_t <- NA
  
  while (t <= max_time) {
    # generate continuous, dichotomize
    cont  <- generate_mixed_continuous(distribution, k, shift_vec, p0)
    bin   <- dichotomize_data(cont, distribution, p0)
    bin_hist[, t] <- bin
    
    # update cumulative statistics
    C_t <- sum(bin)
    cum_sum <- cum_sum + C_t
    
    W_t <- (cum_sum - mu0 * t) / sqrt(t * sigma2_0)
    r_t <- lambda * W_t + (1 - lambda) * r_prev
    r_hist[t] <- r_t
    
    # exact control limits from cached variance
    v_t    <- var_cache[t]
    UCL_t  <-  L * sqrt(v_t)
    LCL_t  <- -L * sqrt(v_t)
    
    if (r_t > UCL_t || r_t < LCL_t) {
      signal_t <- t
      break
    }
    
    r_prev <- r_t
    t <- t + 1
  }
  
  # If no signal, return NULL
  if (is.na(signal_t)) return(NULL)
  
  T_sig <- signal_t
  bin_mat <- bin_hist[, 1:T_sig, drop = FALSE]
  successes <- rowSums(bin_mat)
  
  # ---------------------------------------------
  # Post‑hoc p‑value diagnostics
  # ---------------------------------------------
  pvals <- sapply(1:k, function(i) {
    1 - pbinom(successes[i] - 1, size = T_sig, prob = p0)
  })
  
  # multiple testing corrections
  p_bonf <- p.adjust(pvals, "bonferroni")
  p_holm <- p.adjust(pvals, "holm")
  p_BH   <- p.adjust(pvals, "BH")
  
  # decision rules
  flag_bonf <- p_bonf < alpha
  flag_holm <- p_holm < alpha
  flag_BH   <- p_BH   < alpha
  flag_raw  <- pvals  < alpha
  
  # --- Compute FWER and FDR  ---
  TP_raw  <- sum(flag_raw  & is_ooc)
  FP_raw  <- sum(flag_raw  & !is_ooc)
  TP_bonf <- sum(flag_bonf & is_ooc)
  FP_bonf <- sum(flag_bonf & !is_ooc)
  TP_holm <- sum(flag_holm & is_ooc)
  FP_holm <- sum(flag_holm & !is_ooc)
  TP_BH   <- sum(flag_BH   & is_ooc)
  FP_BH   <- sum(flag_BH   & !is_ooc)
  
  FWER_raw  <- as.numeric(FP_raw  > 0)
  FWER_bonf <- as.numeric(FP_bonf > 0)
  FWER_holm <- as.numeric(FP_holm > 0)
  FWER_bh   <- as.numeric(FP_BH   > 0)
  
  FDR_raw  <- ifelse((TP_raw  + FP_raw)  > 0, FP_raw  / (TP_raw  + FP_raw),  0)
  FDR_bonf <- ifelse((TP_bonf + FP_bonf) > 0, FP_bonf / (TP_bonf + FP_bonf), 0)
  FDR_holm <- ifelse((TP_holm + FP_holm) > 0, FP_holm / (TP_holm + FP_holm), 0)
  FDR_bh   <- ifelse((TP_BH   + FP_BH)   > 0, FP_BH   / (TP_BH   + FP_BH),   0)
  
  # ---------------------------------------------
  # Return metrics
  # ---------------------------------------------
  metrics <- data.frame(
    signal_time = T_sig,
    
    TP_bonf = TP_bonf,
    FP_bonf = FP_bonf,
    FN_bonf = sum(!flag_bonf & is_ooc),
    TN_bonf = sum(!flag_bonf & !is_ooc),
    
    TP_holm = TP_holm,
    FP_holm = FP_holm,
    FN_holm = sum(!flag_holm & is_ooc),
    TN_holm = sum(!flag_holm & !is_ooc),
    
    TP_BH = TP_BH,
    FP_BH = FP_BH,
    FN_BH = sum(!flag_BH & is_ooc),
    TN_BH = sum(!flag_BH & !is_ooc),
    
    TP_raw = TP_raw,
    FP_raw = FP_raw,
    FN_raw = sum(!flag_raw & is_ooc),
    TN_raw = sum(!flag_raw & !is_ooc),
    
    # FWER columns
    FWER_raw  = FWER_raw,
    FWER_bonf = FWER_bonf,
    FWER_holm = FWER_holm,
    FWER_bh   = FWER_bh,
    
    # FDR columns
    FDR_raw  = FDR_raw,
    FDR_bonf = FDR_bonf,
    FDR_holm = FDR_holm,
    FDR_bh   = FDR_bh
  )
  
  list(
    metrics = metrics,
    history = data.frame(t = 1:T_sig, r = r_hist[1:T_sig]),
    flagged = data.frame(
      stream = 1:k,
      successes = successes,
      prop = successes / T_sig,
      p_raw = pvals,
      bonf_sig = flag_bonf,
      holm_sig = flag_holm,
      bh_sig = flag_BH,
      raw_sig = flag_raw,
      truth = is_ooc
    )
  )
}

#============================================================================================
# Parallel simulation core
#============================================================================================

ncores <- parallel::detectCores() - 1
cl <- makeCluster(ncores)
registerDoParallel(cl)

cat("\n=== PARALLEL DIAGNOSTICS STARTED ===\n")
start_time <- Sys.time()

diagnostic_results <- data.frame()
example_store <- list()

for (dist in distributions) {
  cat("\n--- Distribution:", dist, "---\n")
  
  for (i in 1:nrow(optimal_combinations)) {
    
    lambda <- optimal_combinations$lambda[i]
    L      <- optimal_combinations$L[i]
    group  <- optimal_combinations$Group[i]
    target <- optimal_combinations$ARL0_target[i]
    combo  <- optimal_combinations$combo_label[i]
    
    # Precompute exact variance – exact up to t=500, then 1
    var_cache <- precompute_variance(lambda, max_t = max_time_diag, converge_t = 500)
    
    for (delta in delta_diag) {
      for (m in m_values) {
        
        cat(sprintf("  λ=%.3f, L=%.3f, δ=%.2f, m=%d ... ",
                    lambda, L, delta, m))
        flush.console()
        
        # ----------------------------------------------------
        # Parallel simulation with error-protecting tryCatch
        # ----------------------------------------------------
        sim_list <- foreach(s = 1:n_sim_diag, .packages = c("dplyr")) %dopar% {
          tryCatch(
            simulate_one_diagnostic(lambda, L, dist, delta, m,
                                    var_cache = var_cache,
                                    max_time = max_time_diag,
                                    p0 = p0, k = k, alpha = alpha_diag),
            error = function(e) NULL
          )
        }
        
        # ----------------------------------------------------
        # Robust Filtering
        # ----------------------------------------------------
        sim_list <- Filter(function(x) {
          !is.null(x) &&
            is.list(x) &&
            !is.null(x$metrics) &&
            is.data.frame(x$metrics) &&
            ncol(x$metrics) > 0
        }, sim_list)
        
        n_signals <- length(sim_list)
        
        if (n_signals == 0) {
          cat("no valid signals\n")
          next
        }
        
        # ----------------------------------------------------
        # Safe rbind for metrics
        # ----------------------------------------------------
        metrics_df <- tryCatch(
          do.call(rbind, lapply(sim_list, function(x) x$metrics)),
          error = function(e) {
            message("Error binding metrics_df: ", e$message)
            return(NULL)
          }
        )
        
        if (is.null(metrics_df)) {
          cat("metrics_df failed, skipping condition\n")
          next
        }
        
        # ----------------------------------------------------
        # Calculate Performance Metrics (including FWER & FDR)
        # ----------------------------------------------------
        safe_mean <- function(num, den) {
          mean(ifelse(den > 0, num / den, 0), na.rm = TRUE)
        }
        
        perf <- data.frame(
          Distribution = dist,
          lambda = lambda,
          L = L,
          Group = group,
          ARL0_target = target,
          combo_label = combo,
          delta = delta,
          m = m,
          n_signals = n_signals,
          
          FP_bonf = mean(metrics_df$FP_bonf),
          TN_bonf = mean(metrics_df$TN_bonf),
          FP_holm = mean(metrics_df$FP_holm),
          TN_holm = mean(metrics_df$TN_holm),
          FP_BH = mean(metrics_df$FP_BH),
          TN_BH = mean(metrics_df$TN_BH),
          FP_raw = mean(metrics_df$FP_raw),
          TN_raw = mean(metrics_df$TN_raw),
          
          # Sensitivity, Specificity, FDR
          Sensitivity_Bonf = safe_mean(metrics_df$TP_bonf, metrics_df$TP_bonf + metrics_df$FN_bonf),
          Specificity_Bonf = safe_mean(metrics_df$TN_bonf, metrics_df$TN_bonf + metrics_df$FP_bonf),
          FDR_Bonf        = safe_mean(metrics_df$FP_bonf, metrics_df$TP_bonf + metrics_df$FP_bonf),
          
          Sensitivity_Holm = safe_mean(metrics_df$TP_holm, metrics_df$TP_holm + metrics_df$FN_holm),
          Specificity_Holm = safe_mean(metrics_df$TN_holm, metrics_df$TN_holm + metrics_df$FP_holm),
          FDR_Holm        = safe_mean(metrics_df$FP_holm, metrics_df$TP_holm + metrics_df$FP_holm),
          
          Sensitivity_BH = safe_mean(metrics_df$TP_BH, metrics_df$TP_BH + metrics_df$FN_BH),
          Specificity_BH = safe_mean(metrics_df$TN_BH, metrics_df$TN_BH + metrics_df$FP_BH),
          FDR_BH         = safe_mean(metrics_df$FP_BH, metrics_df$TP_BH + metrics_df$FP_BH),
          
          Sensitivity_raw = safe_mean(metrics_df$TP_raw, metrics_df$TP_raw + metrics_df$FN_raw),
          Specificity_raw = safe_mean(metrics_df$TN_raw, metrics_df$TN_raw + metrics_df$FP_raw),
          FDR_raw         = safe_mean(metrics_df$FP_raw, metrics_df$TP_raw + metrics_df$FP_raw),
          
          # FWER
          FWER_bonf = mean(metrics_df$FWER_bonf, na.rm = TRUE),
          FWER_holm = mean(metrics_df$FWER_holm, na.rm = TRUE),
          FWER_bh   = mean(metrics_df$FWER_bh,   na.rm = TRUE),
          FWER_raw  = mean(metrics_df$FWER_raw,  na.rm = TRUE)
        )
        
        diagnostic_results <- rbind(diagnostic_results, perf)
        
        if (length(example_store) < 10)
          example_store[[length(example_store) + 1]] <- sim_list[[1]]
        
        cat(sprintf("done (signals: %d)\n", n_signals))
      }
    }
  }
}

# ================================================
# Stop cluster + report time
# ================================================
stopCluster(cl)
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")
cat(sprintf("\n=== PARALLEL PROCESS COMPLETE (%.2f minutes) ===\n", elapsed))


# ==========================================================
# 5. SAVE RESULTS & SUMMARY TABLES
# ==========================================================
write.csv(diagnostic_results,
          "results/diagnostics/Diagnostic_Performance_Focused.csv",
          row.names = FALSE)

summary_diag <- diagnostic_results %>%
  group_by(ARL0_target, delta, m) %>%
  summarise(across(starts_with("Sensitivity_") |
                     starts_with("Specificity_") |
                     starts_with("FDR_") |
                     starts_with("FWER_"),
                   ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

write.csv(summary_diag,
          "results/diagnostics/Diagnostic_Summary_Focused.csv",
          row.names = FALSE)

# ==========================================================
# TYPE I ERROR SUMMARY TABLE 
# ==========================================================
type1_table <- diagnostic_results %>%
  mutate(
    Type1_Bonf = FP_bonf / (FP_bonf + TN_bonf),
    Type1_Holm = FP_holm / (FP_holm + TN_holm),
    Type1_BH   = FP_BH   / (FP_BH   + TN_BH),
    Type1_raw  = FP_raw  / (FP_raw  + TN_raw)
  ) %>%
  group_by(ARL0_target, delta, m) %>%
  summarise(
    Bonferroni = mean(Type1_Bonf, na.rm = TRUE),
    Holm       = mean(Type1_Holm, na.rm = TRUE),
    BH         = mean(Type1_BH, na.rm = TRUE),
    Raw_p      = mean(Type1_raw, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(type1_table,
          "results/diagnostics/Type1_Error_Table.csv",
          row.names = FALSE)
print(type1_table)

# ==========================================================
# 6. PLOTS Setup (Sensitivity, Specificity, FDR, FWER)
# ==========================================================

# ---- Long data for Sensitivity, Specificity, FDR ----
diag_long <- diagnostic_results %>%
  pivot_longer(
    cols = c(
      Sensitivity_Bonf, Sensitivity_Holm, Sensitivity_BH, Sensitivity_raw,
      Specificity_Bonf, Specificity_Holm, Specificity_BH, Specificity_raw,
      FDR_Bonf, FDR_Holm, FDR_BH, FDR_raw
    ),
    names_to = c("metric", "method"),
    names_pattern = "(.*)_(.*)",
    values_to = "value"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("Bonf", "Holm", "BH", "raw"),
                    labels = c("Bonferroni", "Holm", "BH", "Raw p < 0.05")),
    metric = factor(metric,
                    levels = c("Sensitivity", "Specificity", "FDR"),
                    labels = c("Sensitivity", "Specificity", "FDR")),
    ARL0_label = factor(ARL0_target,
                        levels = unique(ARL0_target),
                        labels = paste0("ARL[0] == ", unique(ARL0_target)))
  )

# ---- Long data for FWER ----
fwer_long <- diagnostic_results %>%
  pivot_longer(
    cols = c(FWER_bonf, FWER_holm, FWER_bh, FWER_raw),
    names_to = "method",
    values_to = "FWER"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("FWER_bonf", "FWER_holm", "FWER_bh", "FWER_raw"),
                    labels = c("Bonferroni", "Holm", "BH", "Raw p < 0.05")),
    ARL0_label = factor(ARL0_target,
                        levels = unique(ARL0_target),
                        labels = paste0("ARL[0] == ", unique(ARL0_target)))
  )

### Custom facet labels for m
label_m <- c(
  "0" = "m == 0~stream",
  "1" = "m == 1~stream",
  "2" = "m == 2~streams",
  "3" = "m == 3~streams"
)

#============================================================
#  PLOTS: Sensitivity, Specificity, FDR 
#============================================================

# Sensitivity plot
p_sens <- ggplot(
  diag_long %>% filter(metric == "Sensitivity"),
  aes(x = delta, y = value, color = method, group = method)
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  labs(
    title = "Sensitivity Across Methods and Shift Magnitudes",
    subtitle = expression("Detection power across different ARL"[0]~"targets and number of streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Sensitivity",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave("results/plots/Sensitivity_LinePlot.png", p_sens, width = 14, height = 7)

# Specificity plot
p_spec <- ggplot(
  diag_long %>% filter(metric == "Specificity"),
  aes(x = delta, y = value, color = method, group = method)
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  labs(
    title = "Specificity Across Methods and Shift Magnitudes",
    subtitle = expression("True negative rate across different ARL"[0]~"targets and number of streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Specificity",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave("results/plots/Specificity_LinePlot.png", p_spec, width = 14, height = 7)

# FDR plot
p_fdr <- ggplot(
  diag_long %>% filter(metric == "FDR"),
  aes(x = delta, y = value, color = method, group = method)
) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  labs(
    title = "False Discovery Rate Across Methods and Shift Magnitudes",
    subtitle = expression("FDR trends across different ARL"[0]~"targets and number of streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "FDR",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave("results/plots/FDR_LinePlot.png", p_fdr, width = 14, height = 7)

#============================================================
#  Heatmap for Sensitivity/Specificity/FDR
#============================================================
p_heat <- ggplot(
  diag_long,
  aes(x = factor(delta), y = method, fill = value)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.5, color = "white", fontface = "bold") +
  facet_grid(
    metric ~ ARL0_label + m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  scale_fill_viridis_c(option = "plasma", name = "Value") +
  labs(
    title = "Performance Heatmap Across Methods",
    subtitle = expression("Sensitivity, Specificity, and FDR by ARL"[0]~"target and number of streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(2, "cm")
  )

ggsave("results/plots/Performance_Heatmap.png", p_heat, width = 16, height = 9)

#============================================================
#  Seperate FDR Heatmap
#============================================================
p_fdr_heat <- ggplot(
  diag_long %>% filter(metric == "FDR"),
  aes(x = factor(delta), y = method, fill = value)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", value)), size = 3.5, color = "white", fontface = "bold") +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  scale_fill_viridis_c(option = "plasma", name = "FDR") +
  labs(
    title = "False Discovery Rate Heatmap",
    subtitle = expression("FDR across methods by ARL"[0]~"target and number of streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(2, "cm")
  )
ggsave("results/plots/FDR_Heatmap.png", p_fdr_heat, width = 16, height = 9)
#============================================================
# NEW: FWER LINE PLOTS
#============================================================
p_fwer <- ggplot(fwer_long,
                 aes(x = delta, y = FWER, color = method, group = method)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  labs(
    title = "Family-wise Error Rate (FWER) Across Methods",
    subtitle = expression("Probability of at least one false positive among in-control streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "FWER",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave("results/plots/FWER_LinePlot.png", p_fwer, width = 14, height = 7)

#============================================================
# NEW: FWER HEATMAP
#============================================================
p_fwer_heat <- ggplot(fwer_long,
                      aes(x = factor(delta), y = method, fill = FWER)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", FWER)), color = "white", size = 3.5, fontface = "bold") +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  scale_fill_viridis_c(option = "magma", name = "FWER") +
  labs(
    title = "Family-Wise Error Rate (FWER) Heatmap",
    subtitle = expression("Probability of at least one false positive across methods"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(2, "cm")
  )

ggsave("results/plots/FWER_Heatmap.png", p_fwer_heat, width = 16, height = 9)

#============================================================
# Type I Error LINE PLOT
#============================================================
type1_long <- type1_table %>%
  pivot_longer(
    cols = c(Bonferroni, Holm, BH, Raw_p),
    names_to = "Method",
    values_to = "Type1"
  ) %>%
  mutate(
    Method = factor(Method,
                    levels = c("Bonferroni","Holm","BH","Raw_p"),
                    labels = c("Bonferroni","Holm","BH","Raw p < 0.05")),
    ARL0_label = factor(ARL0_target,
                        levels = unique(ARL0_target),
                        labels = paste0("ARL[0] == ", unique(ARL0_target)))
  )

p_type1 <- ggplot(type1_long,
                  aes(x = delta, y = Type1, color = Method, group = Method)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  labs(
    title = "Type I Error Across Multiple-Testing Methods",
    subtitle = expression("False positive rate among in-control streams"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Type I Error Rate",
    color = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    panel.spacing = unit(1, "lines")
  )

ggsave("results/plots/Type1_LinePlot.png", p_type1, width = 14, height = 7)


#============================================================
# Type I Error HEATMAP
#============================================================
p_type1_heat <- ggplot(type1_long,
                       aes(x = factor(delta), y = Method, fill = Type1)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", Type1)), color = "white", size = 3.5, fontface = "bold") +
  facet_grid(
    ARL0_label ~ m,
    labeller = labeller(
      ARL0_label = label_parsed,
      m = as_labeller(label_m, label_parsed)
    )
  ) +
  scale_fill_viridis_c(option = "magma", name = "Type I Error") +
  labs(
    title = "Type I Error Heatmap Across Methods",
    subtitle = expression("False positive rate across ARL"[0]~"targets"),
    x = expression("Shift Magnitude (" * delta * ")"),
    y = "Method"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "gray90", color = NA),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.width = unit(2, "cm")
  )

ggsave("results/plots/Type1_Heatmap.png", p_type1_heat, width = 16, height = 9)

cat("\nAll plots saved to 'results/plots/'\n")


# Load required libraries

library(dplyr)
library(tidyr)
library(writexl)

# Assume diagnostic_results is already loaded (from the simulation)
# Create separate tables for each metric

# Sensitivity table
sens_table <- diagnostic_results %>%
  select(ARL0_target, m, delta, starts_with("Sensitivity")) %>%
  pivot_longer(cols = starts_with("Sensitivity"),
               names_to = "Method", values_to = "Sensitivity") %>%
  mutate(Method = gsub("Sensitivity_", "", Method)) %>%
  pivot_wider(names_from = Method, values_from = Sensitivity) %>%
  arrange(ARL0_target, m, delta)

write_xlsx(sens_table, "results/diagnostics/Sensitivity_Table.xlsx")

# Specificity table
spec_table <- diagnostic_results %>%
  select(ARL0_target, m, delta, starts_with("Specificity")) %>%
  pivot_longer(cols = starts_with("Specificity"),
               names_to = "Method", values_to = "Specificity") %>%
  mutate(Method = gsub("Specificity_", "", Method)) %>%
  pivot_wider(names_from = Method, values_from = Specificity) %>%
  arrange(ARL0_target, m, delta)

write_xlsx(spec_table, "results/diagnostics/Specificity_Table.xlsx")

# FDR table
fdr_table <- diagnostic_results %>%
  select(ARL0_target, m, delta, starts_with("FDR")) %>%
  pivot_longer(cols = starts_with("FDR"),
               names_to = "Method", values_to = "FDR") %>%
  mutate(Method = gsub("FDR_", "", Method)) %>%
  pivot_wider(names_from = Method, values_from = FDR) %>%
  arrange(ARL0_target, m, delta)

write_xlsx(fdr_table, "results/diagnostics/FDR_Table.xlsx")

# FWER table (already computed as columns FWER_bonf, FWER_holm, FWER_bh, FWER_raw)
fwer_table <- diagnostic_results %>%
  select(ARL0_target, m, delta, FWER_bonf, FWER_holm, FWER_bh, FWER_raw) %>%
  rename(Bonferroni = FWER_bonf, Holm = FWER_holm, BH = FWER_bh, Raw = FWER_raw) %>%
  arrange(ARL0_target, m, delta)

write_xlsx(fwer_table, "results/diagnostics/FWER_Table.xlsx")

# Type I error table (use type1_table already created)
# type1_table is already in the correct format; save it
write_xlsx(type1_table, "results/diagnostics/Type1_Table.xlsx")