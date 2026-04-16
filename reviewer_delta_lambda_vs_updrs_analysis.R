## Reviewer extension of the motor latent-variable analysis.
## Refits the 33-item OFF-state motor factor model on T1-T4 rows to obtain
## lambda_motor change scores at the manuscript's baseline and six-month visits.

get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  matches <- grep(file_arg, args, value = TRUE)

  if (length(matches) == 0) {
    stop("Could not determine script path from commandArgs().")
  }

  normalizePath(sub(file_arg, "", matches[1]), winslash = "/", mustWork = TRUE)
}

script_path <- get_script_path()
root_dir <- dirname(script_path)
output_dir <- file.path(root_dir, "reviewer_delta_lambda_vs_updrs_output")

.libPaths(c(file.path(root_dir, "Rlibs"), .libPaths()))

suppressMessages(suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(psych)
}))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

load(file.path(root_dir, "tidyEPOdata.Rdata"))

if (!"side" %in% names(subTB) && ncol(subTB) >= 7) {
  colnames(subTB)[7] <- "side"
}

group_lookup <- subTimeTB %>%
  mutate(
    group = case_when(
      time == 2 & Dose == 0 ~ "Placebo",
      time == 2 & Dose == 5 ~ "NeuroEPO",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  distinct(ID, group)

analysis_data <- cognitionTB %>%
  inner_join(motorTB, by = c("ID", "time")) %>%
  filter(time %in% c(1, 2, 3, 4))

motor_matrix <- analysis_data %>%
  select(speech:dyskinesia) %>%
  as.matrix()

fa_fit <- suppressMessages(psych::fa(
  motor_matrix,
  nfactors = 3,
  fm = "minres",
  rotate = "oblimin",
  scores = "regression"
))

if (is.null(fa_fit$scores) || nrow(fa_fit$scores) != nrow(analysis_data)) {
  stop("psych::fa did not return factor scores for all analysis rows.")
}

factor_scores <- as.data.frame(fa_fit$scores)
raw_total <- rowSums(motor_matrix)
factor_fit_corrs <- vapply(
  factor_scores,
  function(score_column) {
    suppressWarnings(cor(score_column, raw_total, use = "complete.obs"))
  },
  numeric(1)
)

fa_model_path <- file.path(output_dir, "delta_lambda_motor_factor_scores.rds")
saveRDS(
  list(
    fa_spec = list(
      nfactors = 3,
      fm = "minres",
      rotate = "oblimin",
      scores = "regression"
    ),
    analysis_rows = analysis_data %>% select(ID, time),
    fa_fit = fa_fit,
    factor_scores = factor_scores
  ),
  fa_model_path
)

anchor_sign <- 1
if (!is.na(factor_fit_corrs[[1]]) && factor_fit_corrs[[1]] < 0) {
  anchor_sign <- -1
}

row_level <- analysis_data %>%
  transmute(
    ID = ID,
    time = time,
    UPDRS_III_OFF_total = raw_total,
    lambda_motor_raw = factor_scores$MR1,
    lambda_motor_anchored = anchor_sign * factor_scores$MR1
  ) %>%
  left_join(group_lookup, by = "ID")

t1_rows <- row_level %>%
  filter(time == 1) %>%
  transmute(
    ID,
    group,
    lambda_motor_T1 = lambda_motor_anchored,
    UPDRS_III_OFF_total_T1 = UPDRS_III_OFF_total
  )

t4_rows <- row_level %>%
  filter(time == 4) %>%
  transmute(
    ID,
    group_t4 = group,
    lambda_motor_T4 = lambda_motor_anchored,
    UPDRS_III_OFF_total_T4 = UPDRS_III_OFF_total
  )

delta_rows <- t1_rows %>%
  inner_join(t4_rows, by = "ID") %>%
  mutate(
    group = coalesce(group, group_t4),
    group = if_else(is.na(group), "Unknown", group),
    delta_lambda_motor_T4_minus_T1 = lambda_motor_T4 - lambda_motor_T1,
    delta_UPDRS_III_OFF_total_T4_minus_T1 = UPDRS_III_OFF_total_T4 - UPDRS_III_OFF_total_T1
  ) %>%
  select(
    ID,
    group,
    lambda_motor_T1,
    lambda_motor_T4,
    delta_lambda_motor_T4_minus_T1,
    UPDRS_III_OFF_total_T1,
    UPDRS_III_OFF_total_T4,
    delta_UPDRS_III_OFF_total_T4_minus_T1
  ) %>%
  arrange(ID)

if (nrow(delta_rows) < 3) {
  stop("Need at least three participants with both T1 and T4 to compute the requested correlations.")
}

delta_rows <- delta_rows %>%
  mutate(group = factor(group, levels = c("Placebo", "NeuroEPO", "Unknown")))

pearson <- cor.test(
  delta_rows$delta_lambda_motor_T4_minus_T1,
  delta_rows$delta_UPDRS_III_OFF_total_T4_minus_T1,
  method = "pearson"
)

spearman <- cor.test(
  delta_rows$delta_lambda_motor_T4_minus_T1,
  delta_rows$delta_UPDRS_III_OFF_total_T4_minus_T1,
  method = "spearman",
  exact = FALSE
)

regression_fit <- lm(
  delta_UPDRS_III_OFF_total_T4_minus_T1 ~ delta_lambda_motor_T4_minus_T1,
  data = delta_rows
)

bootstrap_spearman_ci <- function(data, n_boot = 5000, conf_level = 0.95, seed = 123) {
  set.seed(seed)

  bootstrap_estimates <- replicate(n_boot, {
    sample_idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
    suppressWarnings(cor(
      data$delta_lambda_motor_T4_minus_T1[sample_idx],
      data$delta_UPDRS_III_OFF_total_T4_minus_T1[sample_idx],
      method = "spearman"
    ))
  })

  alpha <- (1 - conf_level) / 2
  unname(quantile(bootstrap_estimates, probs = c(alpha, 1 - alpha), na.rm = TRUE))
}

pearson_ci <- unname(pearson$conf.int)
spearman_ci <- bootstrap_spearman_ci(delta_rows)

format_ci <- function(ci_values, digits = 3) {
  paste0("[", paste(sprintf(paste0("%.", digits, "f"), ci_values), collapse = ", "), "]")
}

format_plot_pvalue <- function(p_value) {
  sprintf("%.6f", p_value)
}

annotation_text <- paste(
  paste0("n = ", nrow(delta_rows)),
  paste0(
    "Pearson r = ",
    sprintf("%.3f", unname(pearson$estimate)),
    "; p = ",
    format_plot_pvalue(pearson$p.value),
    "\n95% CI = ",
    format_ci(pearson_ci)
  ),
  paste0(
    "Spearman rho = ",
    sprintf("%.3f", unname(spearman$estimate)),
    "; p = ",
    format_plot_pvalue(spearman$p.value),
    "\nBootstrap 95% CI = ",
    format_ci(spearman_ci)
  ),
  sep = "\n"
)

x_range <- range(delta_rows$delta_lambda_motor_T4_minus_T1)
y_range <- range(delta_rows$delta_UPDRS_III_OFF_total_T4_minus_T1)
x_padding <- 0.04 * diff(x_range)
y_padding <- 0.04 * diff(y_range)

annotation_x <- x_range[1] + x_padding
annotation_y <- y_range[2] - y_padding

plot_obj <- ggplot(
  delta_rows,
  aes(x = delta_lambda_motor_T4_minus_T1, y = delta_UPDRS_III_OFF_total_T4_minus_T1)
) +
  geom_hline(yintercept = 0, color = "#b0b0b0", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "#b0b0b0", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 3, color = "#1f5f8b") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#c44e52", fill = "#c44e52", alpha = 0.15, linewidth = 0.9) +
  annotate(
    "label",
    x = annotation_x,
    y = annotation_y,
    hjust = 0,
    vjust = 1,
    label = annotation_text,
    size = 3.5,
    linewidth = 0.2,
    fill = "white",
    color = "black"
  ) +
  labs(
    title = "Participant-Level Change in Latent Motor Score\nvs UPDRS-III OFF Total",
    x = "Delta lambda_motor (T4 - T1)",
    y = "Delta UPDRS-III OFF total (T4 - T1)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(16, 14, 10, 10)
  )

group_plot_obj <- ggplot(
  delta_rows,
  aes(
    x = delta_lambda_motor_T4_minus_T1,
    y = delta_UPDRS_III_OFF_total_T4_minus_T1,
    color = group
  )
) +
  geom_hline(yintercept = 0, color = "#b0b0b0", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = 0, color = "#b0b0b0", linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 3) +
  geom_smooth(
    data = delta_rows,
    mapping = aes(
      x = delta_lambda_motor_T4_minus_T1,
      y = delta_UPDRS_III_OFF_total_T4_minus_T1
    ),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#2f2f2f",
    fill = "#2f2f2f",
    alpha = 0.15,
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  annotate(
    "label",
    x = annotation_x,
    y = annotation_y,
    hjust = 0,
    vjust = 1,
    label = annotation_text,
    size = 3.5,
    linewidth = 0.2,
    fill = "white",
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      Placebo = "#1f5f8b",
      NeuroEPO = "#c44e52",
      Unknown = "#6c757d"
    ),
    drop = TRUE
  ) +
  labs(
    title = "Participant-Level Change in Latent Motor Score\nvs UPDRS-III OFF Total by Group",
    x = "Delta lambda_motor (T4 - T1)",
    y = "Delta UPDRS-III OFF total (T4 - T1)",
    color = "Group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(16, 14, 10, 10)
  )

png_path <- file.path(output_dir, "delta_lambda_motor_vs_delta_updrs_T4_T1.png")
pdf_path <- file.path(output_dir, "delta_lambda_motor_vs_delta_updrs_T4_T1.pdf")
group_png_path <- file.path(output_dir, "delta_lambda_motor_vs_delta_updrs_T4_T1_by_group.png")
group_pdf_path <- file.path(output_dir, "delta_lambda_motor_vs_delta_updrs_T4_T1_by_group.pdf")
csv_path <- file.path(output_dir, "participant_delta_lambda_motor_vs_updrs_T4_T1.csv")
summary_path <- file.path(output_dir, "delta_lambda_motor_vs_delta_updrs_summary.txt")

ggsave(png_path, plot = plot_obj, width = 8, height = 6, dpi = 200)
ggsave(pdf_path, plot = plot_obj, width = 8, height = 6)
ggsave(group_png_path, plot = group_plot_obj, width = 8, height = 6, dpi = 200)
ggsave(group_pdf_path, plot = group_plot_obj, width = 8, height = 6)
write.csv(delta_rows, csv_path, row.names = FALSE)

summary_lines <- c(
  "Reviewer delta-lambda analysis",
  "Latent-score method: psych::fa fit directly to the row-level 33-item motor matrix on the joined cognition+motor table used by msmIPWTcommented.R, extended to include T4.",
  "Factor-analysis specification: nfactors = 3, fm = minres, rotate = oblimin, scores = regression.",
  "Original msmIPWTcommented.R subsets to times 1-3, so this reviewer extension refits the same 33-item, 3-factor motor model on times 1-4 to obtain lambda_motor at T4.",
  "First-factor sign anchored so larger lambda_motor corresponds to worse raw UPDRS-III total.",
  "Under that anchor, negative deltas indicate improvement for both lambda_motor and raw UPDRS-III total.",
  paste0(
    "Fit-data Pearson correlations between latent factors and raw total (T1-T4): ",
    paste(sprintf("%.4f", factor_fit_corrs), collapse = ", ")
  ),
  paste0("Participants with both T1 and T4: ", nrow(delta_rows)),
  paste0(
    "Pearson r, 95% CI, p: ",
    sprintf("%.6f", unname(pearson$estimate)),
    ", ",
    format_ci(pearson_ci, digits = 6),
    ", ",
    sprintf("%.6f", pearson$p.value)
  ),
  paste0(
    "Spearman rho, bootstrap 95% CI, p: ",
    sprintf("%.6f", unname(spearman$estimate)),
    ", ",
    format_ci(spearman_ci, digits = 6),
    ", ",
    sprintf("%.6f", spearman$p.value)
  ),
  paste0(
    "Regression line: delta_updrs = ",
    sprintf("%.6f", coef(regression_fit)[1]),
    " + ",
    sprintf("%.6f", coef(regression_fit)[2]),
    " * delta_lambda_motor"
  ),
  paste0("Saved FA scoring object: ", fa_model_path)
)

writeLines(summary_lines, con = summary_path)
cat(paste(summary_lines, collapse = "\n"), "\n", sep = "")
cat("Saved plot: ", png_path, "\n", sep = "")
cat("Saved plot: ", pdf_path, "\n", sep = "")
cat("Saved plot: ", group_png_path, "\n", sep = "")
cat("Saved plot: ", group_pdf_path, "\n", sep = "")
cat("Saved table: ", csv_path, "\n", sep = "")