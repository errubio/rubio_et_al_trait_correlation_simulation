library(MASS)
library(phytools)
library(mvMORPH)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caper)
library(tibble)
library(ape)

n_species <- 500
n_reps    <- 1000
n_indiv   <- 1000

eps_sign  <- 0.05 #avoid noise complication
clamp_cor <- function(r, lo = -0.95, hi = 0.95) pmax(pmin(r, hi), lo) #avoid >1
sign_eps  <- function(x, eps = 0.05) ifelse(!is.finite(x) | abs(x) < eps, 0, sign(x))

#scenarios (correlation and tree depth)
corr_grid <- c(-0.6, -0.3, 0, 0.3, 0.6)
tree_scales_bm <- c(1, 10)

scenario_grid <- tidyr::crossing(
  among_root_cor = corr_grid,
  within_root_cor = corr_grid,
  tree_scale = tree_scales_bm
) %>%
  dplyr::mutate(
    among_root_cor = as.numeric(among_root_cor),
    within_root_cor = as.numeric(within_root_cor),
    tree_scale = as.integer(tree_scale),
    scenario_id = sprintf("A%+.1f_W%+.1f_T%d",
                          among_root_cor, within_root_cor, tree_scale)
  )

output_dir <- file.path(getwd(), "output", "BM")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)


#run one replicate
run_one_rep <- function(rep, scenario_grid, n_species, n_indiv,
                        tree_seed_base = 12345, sim_seed_base = 500000) {
  
  rep_results <- vector("list", nrow(scenario_grid))
  
  #set variance values
  root_means      <- c(bite = 0, sprint = 0)
  root_var_bite   <- 0.071
  root_var_sprint <- 0.183
  
  #create tree
  set.seed(tree_seed_base + rep)
  base_tree <- pbtree(n = n_species, scale = 1)
  
  H <- max(phytools::nodeHeights(base_tree))
  base_tree$edge.length <- base_tree$edge.length / H
  
  stopifnot(ape::is.ultrametric(base_tree))
  tip_names <- base_tree$tip.label
  
  
  for (s in seq_len(nrow(scenario_grid))) {
    
    sc <- scenario_grid[s, ]
    
    scenario_id     <- sc$scenario_id
    among_root_cor  <- sc$among_root_cor
    within_root_cor <- sc$within_root_cor
    tree_scale      <- sc$tree_scale
    
    tree <- base_tree
    tree$edge.length <- base_tree$edge.length * tree_scale
    
    stopifnot(ape::is.ultrametric(tree))
    
    set.seed(sim_seed_base + rep * 1000 + s)
    
    
    #among-species parameters
    rA <- clamp_cor(among_root_cor)
    covA <- rA * sqrt(root_var_bite * root_var_sprint)
    
    root_vcv_means <- matrix(
      c(root_var_bite, covA,
        covA, root_var_sprint),
      nrow = 2, byrow = TRUE,
      dimnames = list(c("bite", "sprint"), c("bite", "sprint"))
    )
    
    trait_means_bm <- mvSIM(
      tree = tree, nsim = 1, model = "BM1",
      param = list(theta = root_means, sigma = root_vcv_means)
    )
    
    trait_means_bm <- as.matrix(trait_means_bm)
    colnames(trait_means_bm) <- c("bite", "sprint")
    rownames(trait_means_bm) <- tip_names
    
    
    #within-species parameters
    rW <- clamp_cor(within_root_cor, lo = -0.999999, hi = 0.999999)
    
    theta_int <- c(
      log_var_bite = log(root_var_bite),
      log_var_sprint = log(root_var_sprint),
      z_cor = atanh(rW)
    )
    
    sigma_int <- diag(c(0.01, 0.01, 0.1))
    rownames(sigma_int) <- colnames(sigma_int) <- names(theta_int)
    
    vcv_sim <- mvSIM(
      tree, nsim = 1, model = "BM1",
      param = list(theta = theta_int, sigma = sigma_int)
    )
    
    vcv_sim <- as.matrix(vcv_sim)
    colnames(vcv_sim) <- names(theta_int)
    rownames(vcv_sim) <- tip_names
    
    var_bite <- exp(vcv_sim[, "log_var_bite"])
    var_sprint <- exp(vcv_sim[, "log_var_sprint"])
    
    target_cor <- tanh(vcv_sim[, "z_cor"])
    target_cor <- clamp_cor(target_cor, lo = -0.999999, hi = 0.999999)
    
    covar <- target_cor * sqrt(var_bite * var_sprint)
    
    
    interspecific_data <- data.frame(
      species = tip_names,
      mean_bite = trait_means_bm[tip_names, "bite"],
      mean_sprint = trait_means_bm[tip_names, "sprint"],
      var_bite = var_bite,
      var_sprint = var_sprint,
      covar = covar
    )
    
    
    #pgls (estimate lambda, fix to 1 if fail)
    rownames(interspecific_data) <- interspecific_data$species
    
    comp_data <- comparative.data(
      tree, interspecific_data,
      names.col = "species", vcv = TRUE
    )
    
    if (var(interspecific_data$mean_bite, na.rm = TRUE) < 1e-12 ||
        var(interspecific_data$mean_sprint, na.rm = TRUE) < 1e-12) {
      
      pgls_ok <- FALSE
      pgls_err <- "Near-zero variance in mean_bite or mean_sprint"
      
      pgls_fit <- pgls(mean_sprint ~ mean_bite, data = comp_data, lambda = 1)
      pgls_sum <- summary(pgls_fit)
      
    } else {
      
      pgls_ok <- TRUE
      pgls_err <- NA_character_
      
      pgls_fit <- tryCatch(
        pgls(mean_sprint ~ mean_bite, data = comp_data, lambda = "ML"),
        error = function(e) {
          pgls_ok <<- FALSE
          pgls_err <<- conditionMessage(e)
          pgls(mean_sprint ~ mean_bite, data = comp_data, lambda = 1)
        }
      )
      
      pgls_sum <- summary(pgls_fit)
    }
    
    pgls_lambda <- NA_real_
    
    if (!is.null(pgls_fit$param) && "lambda" %in% names(pgls_fit$param)) {
      pgls_lambda <- as.numeric(pgls_fit$param[["lambda"]])
    }
    
    pgls_lambda_mode <- ifelse(pgls_ok, "ML", "Fixed(1)")
    
    among_species_slope <- unname(coef(pgls_fit)["mean_bite"])
    among_species_p <- pgls_sum$coefficients["mean_bite", "Pr(>|t|)"]
    among_species_R2 <- pgls_sum$r.squared
    
    
    #among-species correlation and covariance results for rep
    among_species_cor <- cor(
      interspecific_data$mean_bite,
      interspecific_data$mean_sprint,
      use = "pairwise.complete.obs"
    )
    
    among_species_cov <- cov(
      interspecific_data$mean_bite,
      interspecific_data$mean_sprint,
      use = "pairwise.complete.obs"
    )
    
    
    #simulate within-species data
    within_cor <- rep(NA_real_, n_species)
    
    for (i in seq_len(n_species)) {
      
      mu <- c(
        bite = interspecific_data$mean_bite[i],
        sprint = interspecific_data$mean_sprint[i]
      )
      
      vb <- interspecific_data$var_bite[i]
      vs <- interspecific_data$var_sprint[i]
      cov_ws <- interspecific_data$covar[i]
      
      clamp_cov <- sign(cov_ws) *
        min(abs(cov_ws), sqrt(vb * vs) * (1 - 1e-6))
      
      Sigma <- matrix(
        c(vb, clamp_cov,
          clamp_cov, vs),
        nrow = 2, byrow = TRUE
      )
      
      sim_traits <- MASS::mvrnorm(n = n_indiv, mu = mu, Sigma = Sigma)
      
      sd_bite <- sd(sim_traits[, 1])
      sd_sprint <- sd(sim_traits[, 2])
      
      if (sd_bite >= 1e-6 && sd_sprint >= 1e-6) {
        within_cor[i] <- cor(sim_traits[, 1], sim_traits[, 2])
      }
    }
    
    
    within_species_mean_cor <- tanh(
      mean(
        atanh(clamp_cor(within_cor[is.finite(within_cor)],
                        lo = -0.99, hi = 0.99)),
        na.rm = TRUE
      )
    )
    
    within_species_sd_cor <- sd(within_cor, na.rm = TRUE)
    
    
    #sign flip analysis per rep
    sign_flip_rep <-
      (sign_eps(within_species_mean_cor, eps_sign) !=
         sign_eps(among_species_cor, eps_sign)) &
      (sign_eps(within_species_mean_cor, eps_sign) != 0 |
         sign_eps(among_species_cor, eps_sign) != 0)
    
    
    tip_sign_match <-
      (sign_eps(within_cor, eps_sign) ==
         sign_eps(among_species_cor, eps_sign)) &
      (sign_eps(within_cor, eps_sign) != 0 |
         sign_eps(among_species_cor, eps_sign) != 0)
    
    species_match_rate <- mean(tip_sign_match, na.rm = TRUE)
    
    
    rep_results[[s]] <- tibble(
      rep = rep,
      pgls_ml_ok = pgls_ok,
      pgls_error = pgls_err,
      pgls_lambda = as.numeric(pgls_lambda),
      pgls_lambda_mode = pgls_lambda_mode,
      among_species_slope = as.numeric(among_species_slope),
      among_species_p = as.numeric(among_species_p),
      among_species_R2 = as.numeric(among_species_R2),
      among_species_cor = as.numeric(among_species_cor),
      among_species_cov = as.numeric(among_species_cov),
      within_species_mean_cor = as.numeric(within_species_mean_cor),
      within_species_sd_cor = as.numeric(within_species_sd_cor),
      tree_height = max(phytools::nodeHeights(tree)),
      scenario_id = scenario_id,
      among_root_cor = among_root_cor,
      within_root_cor = within_root_cor,
      tree_scale = tree_scale,
      sign_flip_rep = sign_flip_rep,
      species_match_rate = species_match_rate
    )
  }
  
  bind_rows(rep_results)
}


#run scenarios
run_start_time <- Sys.time()

for (rep in seq_len(n_reps)) {
  
  checkpoint_file <- file.path(
    checkpoint_dir,
    sprintf("BM_rep_%04d.rds", rep)
  )
  
  if (file.exists(checkpoint_file)) {
    cat("Replicate", rep, "of", n_reps, "already completed - skipping\n")
    next
  }
  
  rep_start_time <- Sys.time()
  
  cat("\nStarting replicate", rep, "of", n_reps,
      "at", format(Sys.time()), "\n")
  
  rep_results <- run_one_rep(
    rep = rep,
    scenario_grid = scenario_grid,
    n_species = n_species,
    n_indiv = n_indiv
  )
  
  temp_file <- paste0(checkpoint_file, ".tmp")
  saveRDS(rep_results, temp_file)
  
  if (file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
  }
  
  file.rename(temp_file, checkpoint_file)
  
  rep_elapsed <- as.numeric(
    difftime(Sys.time(), rep_start_time, units = "mins")
  )
  
  completed_files <- length(
    list.files(
      checkpoint_dir,
      pattern = "^BM_rep_[0-9]{4}\\.rds$"
    )
  )
  
  elapsed_hours <- as.numeric(
    difftime(Sys.time(), run_start_time, units = "hours")
  )
  
  cat("Completed replicate", rep, "of", n_reps, "\n")
  cat("Replicate time:", round(rep_elapsed, 2), "minutes\n")
  cat("Completed checkpoints:", completed_files, "of", n_reps, "\n")
  cat("Elapsed this session:", round(elapsed_hours, 2), "hours\n")
  
  rm(rep_results)
  gc()
}


checkpoint_files <- list.files(
  checkpoint_dir,
  pattern = "^BM_rep_[0-9]{4}\\.rds$",
  full.names = TRUE
)

checkpoint_files <- sort(checkpoint_files)

results_all <- bind_rows(lapply(checkpoint_files, readRDS))


#scenarios summary
scenario_summary_df <- results_all %>%
  group_by(
    scenario_id,
    among_root_cor,
    within_root_cor,
    tree_scale
  ) %>%
  summarise(
    n = n(),
    ml_success_rate = mean(pgls_ml_ok, na.rm = TRUE),
    among_species_cor_mean = mean(among_species_cor, na.rm = TRUE),
    within_species_cor_mean = mean(within_species_mean_cor, na.rm = TRUE),
    mean_diff = mean(within_species_mean_cor - among_species_cor, na.rm = TRUE),
    flip_rate = mean(sign_flip_rep, na.rm = TRUE),
    track_r = cor(among_species_cor, within_species_mean_cor, use = "complete.obs"),
    mean_tree_height = mean(tree_height, na.rm = TRUE),
    species_match_rate = mean(species_match_rate, na.rm = TRUE),
    .groups = "drop"
  )


saveRDS(results_all, file.path(output_dir, "BM_results_final.rds"))
saveRDS(scenario_summary_df, file.path(output_dir, "BM_scenario_summary_final.rds"))

write.csv(results_all, file.path(output_dir, "BM_results_final.csv"), row.names = FALSE)
write.csv(scenario_summary_df, file.path(output_dir, "BM_scenario_summary_final.csv"), row.names = FALSE)



#analyses

results_for_lm <- results_all %>%
  mutate(
    abs_diff = abs(within_species_mean_cor - among_species_cor),
    
    sign_relation = case_when(
      among_root_cor == 0 | within_root_cor == 0 ~ "Zero-involved",
      among_root_cor * within_root_cor > 0 ~ "Same sign",
      among_root_cor * within_root_cor < 0 ~ "Opposite sign"
    ),
    
    sign_relation = factor(
      sign_relation,
      levels = c("Same sign", "Zero-involved", "Opposite sign")
    ),
    
    tree_scale = factor(
      tree_scale,
      levels = c(1, 10)
    )
  )


#full mismatch magnitude model
mismatch_model <- lm(
  abs_diff ~ sign_relation * tree_scale,
  data = results_for_lm
)

bm_anova <- anova(mismatch_model)
bm_anova
summary(mismatch_model)

#partial eta squared 
bm_anova_effects <- bm_anova[!is.na(bm_anova$`F value`), , drop = FALSE]
bm_eta_p2 <- (bm_anova_effects$`F value` * bm_anova_effects$Df) /
  ((bm_anova_effects$`F value` * bm_anova_effects$Df) + df.residual(mismatch_model))
bm_eta_p2


#full sign flip model
flip_model_results <- glm(
  sign_flip_rep ~ sign_relation * tree_scale,
  data = results_for_lm,
  family = binomial
)

summary(flip_model_results)
anova(flip_model_results, test = "Chisq")


#mismatch magnitude summary
mismatch_summary <- results_for_lm %>%
  group_by(sign_relation, tree_scale) %>%
  summarise(
    n = n(),
    mean_mismatch = mean(abs_diff, na.rm = TRUE),
    sd_mismatch = sd(abs_diff, na.rm = TRUE),
    median_mismatch = median(abs_diff, na.rm = TRUE),
    iqr_low = quantile(abs_diff, 0.25, na.rm = TRUE),
    iqr_high = quantile(abs_diff, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(mean_mismatch, sd_mismatch, median_mismatch, iqr_low, iqr_high),
      ~ round(.x, 3)
    )
  )

mismatch_summary


#replicate sign flip probability
flip_prob_summary <- results_for_lm %>%
  group_by(sign_relation, tree_scale) %>%
  summarise(
    n = n(),
    n_flips = sum(sign_flip_rep, na.rm = TRUE),
    flip_rate = mean(sign_flip_rep, na.rm = TRUE),
    flip_percent = round(100 * flip_rate, 1),
    .groups = "drop"
  )

flip_prob_summary

#matched scenario sign flip summary
matched_flip_summary <- results_all %>%
  filter(
    among_root_cor == within_root_cor,
    among_root_cor != 0
  ) %>%
  group_by(tree_scale, among_root_cor) %>%
  summarise(
    n = n(),
    n_flips = sum(sign_flip_rep, na.rm = TRUE),
    flip_rate = mean(sign_flip_rep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    flip_percent = round(100 * flip_rate, 1)
  )

matched_flip_summary

#phylogenetic signal summary
lambda_summary <- results_all %>%
  group_by(tree_scale) %>%
  summarise(
    n = n(),
    mean_lambda = mean(pgls_lambda, na.rm = TRUE),
    sd_lambda = sd(pgls_lambda, na.rm = TRUE),
    median_lambda = median(pgls_lambda, na.rm = TRUE),
    iqr_low = quantile(pgls_lambda, 0.25, na.rm = TRUE),
    iqr_high = quantile(pgls_lambda, 0.75, na.rm = TRUE),
    ml_success_rate = mean(pgls_lambda_mode == "ML", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        mean_lambda,
        sd_lambda,
        median_lambda,
        iqr_low,
        iqr_high,
        ml_success_rate
      ),
      ~ round(.x, 3)
    ),
    
    ml_success_percent = round(100 * ml_success_rate, 1)
  )

lambda_summary

#sign flip plot
flip_plot_data <- results_for_lm %>%
  mutate(
    tree_scale = factor(
      as.numeric(as.character(tree_scale)),
      levels = c(1, 10),
      labels = c("Tree depth = 1", "Tree depth = 10")
    ),
    
    sign_relation = factor(
      sign_relation,
      levels = c("Same sign", "Zero-involved", "Opposite sign")
    )
  ) %>%
  group_by(sign_relation, tree_scale) %>%
  summarise(
    n = n(),
    n_flips = sum(sign_flip_rep, na.rm = TRUE),
    flip_percent = 100 * mean(sign_flip_rep, na.rm = TRUE),
    .groups = "drop"
  )

flip_plot <- ggplot(
  flip_plot_data,
  aes(
    x = sign_relation,
    y = flip_percent,
    fill = sign_relation
  )
) +
  geom_col(width = 0.7) +
  facet_wrap(~ tree_scale, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Same sign" = "#4C7434",
      "Zero-involved" = "#2F7F73",
      "Opposite sign" = "#1F4E79"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Correlation Scenario",
    y = "Sign-Flip Frequency (%)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 8.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    panel.grid.minor = element_blank()
  )

flip_plot

ggsave(
  filename = file.path(output_dir, "BM sign flip plot.png"),
  plot = flip_plot,
  width = 8,
  height = 5,
  dpi = 600
)

#mismatch magnitude plot
mismatch_plot_data <- results_for_lm %>%
  mutate(
    tree_scale = factor(
      as.numeric(as.character(tree_scale)),
      levels = c(1, 10),
      labels = c("Tree depth = 1", "Tree depth = 10")
    ),
    
    sign_relation = factor(
      sign_relation,
      levels = c("Same sign", "Zero-involved", "Opposite sign")
    )
  )

mismatch_plot <- ggplot(
  mismatch_plot_data,
  aes(
    x = sign_relation,
    y = abs_diff,
    fill = sign_relation
  )
) +
  geom_boxplot(
    width = 0.7,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  facet_wrap(~ tree_scale, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Same sign" = "#4C7434",
      "Zero-involved" = "#2F7F73",
      "Opposite sign" = "#1F4E79"
    )
  ) +
  coord_cartesian(
    ylim = c(0, 2)
  ) +
  labs(
    x = "Correlation Scenario",
    y = "Correlation Mismatch Magnitude"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(
      face = "bold",
      size = 11,
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      face = "bold",
      size = 11,
      margin = margin(r = 10)
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 8.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    panel.grid.minor = element_blank()
  )

mismatch_plot

ggsave(
  filename = file.path(output_dir, "BM mismatch plot.png"),
  plot = mismatch_plot,
  width = 8,
  height = 5,
  dpi = 600
)
