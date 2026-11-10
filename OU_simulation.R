library(MASS)
library(phytools)
library(mvMORPH)
library(dplyr)
library(tidyr)
library(ggplot2)
library(caper)
library(tibble)

n_species <- 500 
n_reps <- 1000 
n_indiv <- 1000 
regime_names <- c("A","B","C") 
tree_scales <- c(1, 10) #tree depth

output_dir <- file.path(getwd(), "output", "OU")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#helpers
paint_random_regimes <- function(tree, min_tips = 15, max_tries = 200) {
  nTips <- Ntip(tree)
  internal_nodes <- (nTips + 2):(nTips + Nnode(tree))
  
  for (k in 1:max_tries) {
    nodes <- sample(internal_nodes, 3, replace = FALSE)
    
    tree_p <- paintSubTree(tree, nodes[1], state = "A", anc.state = "A", stem = TRUE)
    tree_p <- paintSubTree(tree_p, nodes[2], state = "B", anc.state = "A", stem = TRUE)
    tree_p <- paintSubTree(tree_p, nodes[3], state = "C", anc.state = "A", stem = TRUE)
    
    tip_states <- getStates(tree_p, type = "tips")
    tab <- table(tip_states)
    
    if (all(c("A","B","C") %in% names(tab)) && all(tab[c("A","B","C")] >= min_tips)) {
      return(list(tree_simmap = tree_p, tip_states = tip_states))
    }
  }
  
  stop("could not generate regimes with enough tips")
}

rescale_simmap <- function(tree, target_depth) {
  current_depth <- max(nodeHeights(tree))
  scale_factor <- target_depth / current_depth

  tree$edge.length <- tree$edge.length * scale_factor

  if (!is.null(tree$maps)) {
    tree$maps <- lapply(tree$maps, function(x) x * scale_factor)
  }

  if (!is.null(tree$mapped.edge)) {
    tree$mapped.edge <- tree$mapped.edge * scale_factor
  }

  tree
}

clamp_cor <- function(r, lo = -0.99, hi = 0.99) {
  pmax(pmin(r, hi), lo)
}

sign3 <- function(x, eps = 0.05) {
  dplyr::case_when(
    !is.finite(x) ~ NA_integer_,
    abs(x) < eps  ~ 0L,
    x > 0 ~ 1L,
    TRUE ~ -1L
  )
}

make_theta_vcv <- function(within_structure, root_cor,
                           aligned_spread = 0.10,
                           var_bite = 0.071,
                           var_sprint = 0.183) {
  
  root_cor <- clamp_cor(root_cor, lo = -0.99, hi = 0.99)
  base_sign <- ifelse(root_cor < 0, -1, ifelse(root_cor > 0, 1, 0))
  
  if (within_structure == "aligned") {
    rA <- root_cor - aligned_spread
    rB <- root_cor
    rC <- root_cor + aligned_spread
  } else if (within_structure == "mismatched") {
    if (base_sign == 0) {
      rA <- -0.20
      rB <-  0.00
      rC <-  0.20
    } else {
      rA <- -root_cor - 0.05 * base_sign
      rB <- -root_cor
      rC <- -root_cor + 0.05 * base_sign
    }
  } else {
    rA <- -0.30
    rB <- -0.10
    rC <-  0.30
  }
  
  rA <- clamp_cor(rA, lo = -0.95, hi = 0.95)
  rB <- clamp_cor(rB, lo = -0.95, hi = 0.95)
  rC <- clamp_cor(rC, lo = -0.95, hi = 0.95)
  
  theta_vcv <- t(rbind(
    A = c(log_var_bite = log(var_bite),
          log_var_sprint = log(var_sprint),
          correlation_z = atanh(rA)),
    B = c(log_var_bite = log(var_bite),
          log_var_sprint = log(var_sprint),
          correlation_z = atanh(rB)),
    C = c(log_var_bite = log(var_bite),
          log_var_sprint = log(var_sprint),
          correlation_z = atanh(rC))
  ))
  
  t(theta_vcv)
}

#within-species VCV parameters
alpha_vcv <- diag(c(0.1, 0.1, 0.8))
sigma_vcv <- diag(c(0.01, 0.01, 0.005))
rownames(alpha_vcv) <- colnames(alpha_vcv) <- c("log_var_bite","log_var_sprint","correlation_z")
rownames(sigma_vcv) <- colnames(sigma_vcv) <- c("log_var_bite","log_var_sprint","correlation_z")

#attraction strength scenarios 
alpha_scenarios <- tibble(
  alpha_scenario = c("Weak attraction", "Intermediate attraction", "Strong attraction"),
  alpha_value = c(0.05, 0.50, 1.00)
)

root_cor_levels <- c(-0.6, -0.3, 0, 0.3, 0.6)
within_structure_levels <- c("aligned", "mismatched", "regime-divergent")

scenario_grid <- tidyr::crossing(
  root_cor = root_cor_levels,
  within_structure = within_structure_levels,
  alpha_scenarios
) %>%
  mutate(
    scenario_id = row_number(),
    root_cor = as.numeric(root_cor),
    within_structure = as.character(within_structure),
    alpha_scenario = as.character(alpha_scenario),
    alpha_value = as.numeric(alpha_value)
  )

n_sc <- nrow(scenario_grid)
rep_results <- vector("list", n_reps * length(tree_scales) * n_sc)
idx <- 0L
example_trees <- list(`1` = NULL, `10` = NULL)

for (rep in 1:n_reps) {
  
  set.seed(123 + rep * 1000)
  base_tree <- pbtree(n = n_species, scale = 1)
  
  painted <- tryCatch(
    paint_random_regimes(base_tree),
    error = function(e) NULL
  )
  
  if (is.null(painted)) next
  
  base_tree_simmap <- painted$tree_simmap
  tip_regimes <- painted$tip_states
  
  for (T in tree_scales) {
    
    tree_simmap <- rescale_simmap(base_tree_simmap, target_depth = T)
    
    if (is.null(example_trees[[as.character(T)]])) {
      example_trees[[as.character(T)]] <- list(
        tree = tree_simmap,
        tip_regimes = tip_regimes,
        rep = rep,
        tree_scale = T
      )
    }
    
    for (s in 1:n_sc) {
      idx <- idx + 1L
      sc <- scenario_grid[s, ]
      
      scenario_id_s <- pull(sc, scenario_id)
      root_cor_s <- pull(sc, root_cor)
      within_structure_s <- pull(sc, within_structure)
      alpha_scenario_s <- pull(sc, alpha_scenario)
      alpha_value_s <- pull(sc, alpha_value)
      
      root_cor_s <- clamp_cor(root_cor_s, lo = -0.95, hi = 0.95)
      
      #species mean trait optima
      theta_traits <- matrix(
        c(-0.01, 0.00, 0.01,
          -0.01, 0.00, 0.01),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("bite","sprint"), regime_names)
      )
      
      #species mean alpha
      alpha_traits <- diag(c(alpha_value_s, alpha_value_s))
      rownames(alpha_traits) <- colnames(alpha_traits) <- c("bite","sprint")
      
      #species mean covariance
      var_bite_macro <- 0.071
      var_sprint_macro <- 0.183
      cov_macro <- root_cor_s * sqrt(var_bite_macro * var_sprint_macro)
      
      sigma_traits <- matrix(
        c(var_bite_macro, cov_macro,
          cov_macro, var_sprint_macro),
        nrow = 2, byrow = TRUE,
        dimnames = list(c("bite","sprint"), c("bite","sprint"))
      )
      
      #simulate species mean trait values 
      trait_sim <- mvSIM(
        tree_simmap, nsim = 1, model = "OUM",
        param = list(theta = theta_traits, alpha = alpha_traits, sigma = sigma_traits)
      )
      
      trait_sim <- as.matrix(trait_sim)
      trait_means <- trait_sim[tree_simmap$tip.label, , drop = FALSE]
      colnames(trait_means) <- c("bite","sprint")
      
      bite_means <- trait_means[, "bite"]
      sprint_means <- trait_means[, "sprint"]
      
      names(bite_means) <- tree_simmap$tip.label
      names(sprint_means) <- tree_simmap$tip.label
      
      #estimate phylo signal for species means
      lambda_bite_hat <- tryCatch(
        suppressWarnings(
          phytools::phylosig(tree_simmap, bite_means, method = "lambda", test = FALSE)
        ),
        error = function(e) NA_real_
      )
      if (is.list(lambda_bite_hat)) lambda_bite_hat <- lambda_bite_hat$lambda
      lambda_bite_hat <- as.numeric(lambda_bite_hat)[1]
      
      lambda_sprint_hat <- tryCatch(
        suppressWarnings(
          phytools::phylosig(tree_simmap, sprint_means, method = "lambda", test = FALSE)
        ),
        error = function(e) NA_real_
      )
      if (is.list(lambda_sprint_hat)) lambda_sprint_hat <- lambda_sprint_hat$lambda
      lambda_sprint_hat <- as.numeric(lambda_sprint_hat)[1]
      
      #evolve within-species VCV parameters
      theta_vcv_rep <- make_theta_vcv(
        within_structure = within_structure_s,
        root_cor = root_cor_s,
        aligned_spread = 0.10
      )
      
      vcv_sim_raw <- mvSIM(
        tree_simmap, nsim = 1, model = "OUM",
        param = list(theta = theta_vcv_rep, alpha = alpha_vcv, sigma = sigma_vcv)
      )
      
      vcv_sim_raw <- as.matrix(vcv_sim_raw)
      colnames(vcv_sim_raw) <- colnames(theta_vcv_rep)
      
      #back-transform VCV parameters
      vcv_params <- tibble(species = rownames(vcv_sim_raw)) %>%
        mutate(
          regime = tip_regimes[species],
          var_bite = exp(vcv_sim_raw[, "log_var_bite"]),
          var_sprint = exp(vcv_sim_raw[, "log_var_sprint"]),
          correlation = clamp_cor(tanh(vcv_sim_raw[, "correlation_z"]), lo = -0.95, hi = 0.95),
          covar = correlation * sqrt(var_bite * var_sprint)
        ) %>%
        dplyr::select(species, regime, var_bite, var_sprint, correlation, covar)
      
      #simulate individuals (from each species' mean vector and VCV matrix)
      indiv_data <- vector("list", Ntip(tree_simmap))
      
      for (i in seq_len(Ntip(tree_simmap))) {
        sp <- tree_simmap$tip.label[i]
        
        mu <- c(
          bite = trait_means[sp, "bite"],
          sprint = trait_means[sp, "sprint"]
        )
        
        var_bite_i <- vcv_params$var_bite[match(sp, vcv_params$species)]
        var_sprint_i <- vcv_params$var_sprint[match(sp, vcv_params$species)]
        covar_i <- vcv_params$covar[match(sp, vcv_params$species)]
        
        covar_i <- sign(covar_i) * min(
          abs(covar_i),
          sqrt(var_bite_i * var_sprint_i) * (1 - 1e-6)
        )
        
        Sigma <- matrix(
          c(var_bite_i, covar_i,
            covar_i, var_sprint_i),
          nrow = 2, byrow = TRUE
        )
        
        sim_traits <- MASS::mvrnorm(n = n_indiv, mu = mu, Sigma = Sigma)
        colnames(sim_traits) <- c("bite","sprint")
        
        indiv_data[[i]] <- data.frame(
          rep = rep,
          tree_scale = T,
          scenario_id = scenario_id_s,
          alpha_scenario = alpha_scenario_s,
          alpha_value = alpha_value_s,
          root_cor = root_cor_s,
          within_structure = within_structure_s,
          species = sp,
          bite = sim_traits[, "bite"],
          sprint = sim_traits[, "sprint"]
        )
      }
      
      intra_data <- bind_rows(indiv_data)
      
      #within-species correlations and slopes
      within_by_species <- intra_data %>%
        group_by(
          rep, tree_scale, scenario_id, alpha_scenario, alpha_value,
          root_cor, within_structure, species
        ) %>%
        summarise(
          within_cor = if (sd(bite) < 1e-6 || sd(sprint) < 1e-6) {
            NA_real_
          } else {
            cor(bite, sprint)
          },
          .groups = "drop"
        )
      
      #mean within-species correlation (fisher z)
      within_mean_cor <- within_by_species %>%
        filter(is.finite(within_cor)) %>%
        mutate(z = atanh(clamp_cor(within_cor, lo = -0.999999, hi = 0.999999))) %>%
        summarise(rbar = tanh(mean(z, na.rm = TRUE)), .groups = "drop") %>%
        pull(rbar)
      
      if (length(within_mean_cor) == 0) within_mean_cor <- NA_real_
      
      #species means from evolved trait values
      species_means <- data.frame(
        species = tree_simmap$tip.label,
        bite_mean = trait_means[tree_simmap$tip.label, "bite"],
        sprint_mean = trait_means[tree_simmap$tip.label, "sprint"]
      )
      
      #among-species correlation
      among_cor <- with(species_means, cor(bite_mean, sprint_mean))
      
      #pgls among species
      comp_dat <- tryCatch(
        caper::comparative.data(
          phy = tree_simmap,
          data = as.data.frame(species_means),
          names.col = "species",
          vcv = TRUE,
          na.omit = FALSE,
          warn.dropped = FALSE
        ),
        error = function(e) NULL
      )
      
      pgls_fit <- tryCatch(
        if (is.null(comp_dat)) NULL else {
          caper::pgls(sprint_mean ~ bite_mean, data = comp_dat, lambda = "ML")
        },
        error = function(e) NULL
      )
      
      pgls_lambda <- if (is.null(pgls_fit)) {
        NA_real_
      } else {
        as.numeric(pgls_fit$param["lambda"])
      }
      
      rep_results[[idx]] <- tibble(
        rep = rep,
        tree_scale = T,
        scenario_id = scenario_id_s,
        alpha_scenario = alpha_scenario_s,
        alpha_value = alpha_value_s,
        root_cor = root_cor_s,
        within_structure = within_structure_s,
        lambda_bite_hat = lambda_bite_hat,
        lambda_sprint_hat = lambda_sprint_hat,
        among_cor = as.numeric(among_cor),
        within_mean_cor = as.numeric(within_mean_cor),
        pgls_lambda = as.numeric(pgls_lambda),
        nA = as.integer(sum(tip_regimes == "A")),
        nB = as.integer(sum(tip_regimes == "B")),
        nC = as.integer(sum(tip_regimes == "C"))
      )
    }
  }
}

results_ou <- bind_rows(rep_results)
saveRDS(results_ou, file.path(output_dir, "OU_results_final.rds"))

cols_num <- c(
  "rep","tree_scale","scenario_id","alpha_value",
  "root_cor","lambda_bite_hat","lambda_sprint_hat",
  "among_cor","within_mean_cor","pgls_lambda",
  "nA","nB","nC"
)

for (nm in cols_num) {
  if (nm %in% names(results_ou) && is.list(results_ou[[nm]])) {
    results_ou[[nm]] <- vapply(results_ou[[nm]], function(z) {
      if (is.null(z) || length(z) == 0) return(NA_real_)
      as.numeric(z)[1]
    }, numeric(1))
  }
}

cols_chr <- c("alpha_scenario","within_structure")
for (nm in cols_chr) {
  if (nm %in% names(results_ou) && is.list(results_ou[[nm]])) {
    results_ou[[nm]] <- vapply(results_ou[[nm]], function(z) {
      if (is.null(z) || length(z) == 0) return(NA_character_)
      as.character(z)[1]
    }, character(1))
  }
}

results_clean <- results_ou %>%
  filter(
    is.finite(within_mean_cor),
    is.finite(among_cor)
  ) %>%
  mutate(
    abs_diff = abs(within_mean_cor - among_cor),
    diff_signed = within_mean_cor - among_cor,
    sign_within = sign3(within_mean_cor, eps = 0.05),
    sign_among = sign3(among_cor, eps = 0.05),
    sign_match = (sign_within == sign_among),
    sign_flip = !sign_match,
    tree_scale = factor(tree_scale, levels = tree_scales),
    root_cor_f = factor(root_cor, levels = root_cor_levels),
    within_structure = factor(within_structure, levels = within_structure_levels),
    alpha_scenario = factor(
      alpha_scenario,
      levels = c("Weak attraction", "Intermediate attraction", "Strong attraction")
    )
  ) %>%
  filter(!is.na(sign_match))

model_data <- results_clean


#analyses

results_clean <- results_ou %>%
  filter(
    is.finite(within_mean_cor),
    is.finite(among_cor)
  ) %>%
  mutate(
    abs_diff = abs(within_mean_cor - among_cor),
    diff_signed = within_mean_cor - among_cor,
    sign_within = sign3(within_mean_cor, eps = 0.05),
    sign_among = sign3(among_cor, eps = 0.05),
    sign_match = (sign_within == sign_among),
    sign_flip = !sign_match,
    tree_scale = factor(tree_scale, levels = tree_scales),
    root_cor_f = factor(root_cor, levels = root_cor_levels),
    within_structure = factor(within_structure, levels = within_structure_levels),
    alpha_scenario = factor(
      alpha_scenario,
      levels = c("Weak attraction", "Intermediate attraction", "Strong attraction")
    )
  ) %>%
  filter(!is.na(sign_match))

model_data <- results_clean

#mismatch magnitude
mismatch_model <- lm(
  abs_diff ~ within_structure * alpha_scenario + root_cor_f + tree_scale,
  data = model_data
)

ou_anova <- anova(mismatch_model)
ou_anova
summary(mismatch_model)

#partial eta squared 
ou_anova_effects <- ou_anova[!is.na(ou_anova$`F value`), , drop = FALSE]
ou_eta_p2 <- (ou_anova_effects$`F value` * ou_anova_effects$Df) /
  ((ou_anova_effects$`F value` * ou_anova_effects$Df) + df.residual(mismatch_model))
ou_eta_p2

#replicate sign flips
flip_model <- glm(
  sign_flip ~ within_structure * alpha_scenario + root_cor_f + tree_scale,
  data = model_data,
  family = binomial
)

summary(flip_model)
anova(flip_model, test = "Chisq")

#phylogenetic signal summary 
lambda_summary <- model_data %>%
  group_by(alpha_scenario, tree_scale) %>%
  summarise(
    n = n(),
    mean_lambda_bite = mean(lambda_bite_hat, na.rm = TRUE),
    sd_lambda_bite = sd(lambda_bite_hat, na.rm = TRUE),
    mean_lambda_sprint = mean(lambda_sprint_hat, na.rm = TRUE),
    sd_lambda_sprint = sd(lambda_sprint_hat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(mean_lambda_bite, sd_lambda_bite,
        mean_lambda_sprint, sd_lambda_sprint),
      ~ round(.x, 3)
    )
  )

lambda_summary

#sign flip summary by root correlation
flip_prob_summary <- model_data %>%
  group_by(within_structure, root_cor_f, alpha_scenario, tree_scale) %>%
  summarise(
    n = n(),
    n_flips = sum(sign_flip, na.rm = TRUE),
    flip_rate = mean(sign_flip, na.rm = TRUE),
    flip_percent = round(100 * flip_rate, 1),
    .groups = "drop"
  )

flip_prob_summary

#mismatch summary 
mismatch_summary <- model_data %>%
  group_by(within_structure, alpha_scenario, tree_scale) %>%
  summarise(
    n = n(),
    mean_mismatch = mean(abs_diff, na.rm = TRUE),
    sd_mismatch = sd(abs_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    across(c(mean_mismatch, sd_mismatch), ~ round(.x, 3))
  )

mismatch_summary


#sign flip plot
flip_plot_data <- model_data %>%
  mutate(
    tree_scale = factor(as.numeric(as.character(tree_scale)),
                        levels = c(1, 10),
                        labels = c("Tree depth = 1", "Tree depth = 10")),
    alpha_scenario = dplyr::recode(as.character(alpha_scenario),
                                   "Weak attraction" = "Weak",
                                   "Intermediate attraction" = "Intermediate",
                                   "Strong attraction" = "Strong"),
    alpha_scenario = factor(alpha_scenario,
                            levels = c("Weak", "Intermediate", "Strong")),
    within_structure_plot = dplyr::recode(as.character(within_structure),
                                          "aligned" = "Aligned",
                                          "mismatched" = "Mismatched",
                                          "regime-divergent" = "Regime-\ndivergent"),
    within_structure_plot = factor(within_structure_plot,
                                   levels = c("Aligned", "Regime-\ndivergent", "Mismatched"))
  ) %>%
  group_by(within_structure_plot, alpha_scenario, tree_scale) %>%
  summarise(
    n = n(),
    n_flips = sum(sign_flip, na.rm = TRUE),
    flip_percent = 100 * mean(sign_flip, na.rm = TRUE),
    .groups = "drop"
  )

flip_plot <- ggplot(
  flip_plot_data,
  aes(x = within_structure_plot, y = flip_percent, fill = alpha_scenario)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~ tree_scale, nrow = 1) +
  scale_fill_manual(values = c(
    "Weak" = "#4C7434",
    "Intermediate" = "#3E8C7A",
    "Strong" = "#1F4E79"
  )) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = "Intraspecific Covariance Scenario",
    y = "Sign-Flip Frequency (%)",
    fill = "Attraction Strength"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 10, hjust = 0.5, size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(face = "bold", size = 8),
    legend.text = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_dir, "OU sign flip plot.png"), plot = flip_plot, width = 8, height = 5)

flip_plot


#mismatch magnitude (all scenarios)
mismatch_plot_data <- model_data %>%
  mutate(
    tree_scale = factor(as.numeric(as.character(tree_scale)),
                        levels = c(1, 10),
                        labels = c("Tree depth = 1", "Tree depth = 10")),
    alpha_scenario = dplyr::recode(as.character(alpha_scenario),
                                   "Weak attraction" = "Weak",
                                   "Intermediate attraction" = "Intermediate",
                                   "Strong attraction" = "Strong"),
    alpha_scenario = factor(alpha_scenario, levels = c("Weak", "Intermediate", "Strong")),
    within_structure_plot = dplyr::recode(as.character(within_structure),
                                          "aligned" = "Aligned",
                                          "mismatched" = "Mismatched",
                                          "regime-divergent" = "Regime-divergent"),
    within_structure_plot = factor(within_structure_plot,
                                   levels = c("Aligned", "Regime-divergent", "Mismatched"))
  )

mismatch_plot <- ggplot(
  mismatch_plot_data,
  aes(x = within_structure_plot, y = abs_diff, fill = alpha_scenario)
) +
  geom_boxplot(position = position_dodge(width = 0.8),
               width = 0.7, outlier.shape = NA, alpha = 0.8) +
  facet_wrap(~ tree_scale, nrow = 1) +
  scale_x_discrete(labels = c(
    "Aligned" = "Aligned",
    "Regime-divergent" = "Regime-\ndivergent",
    "Mismatched" = "Mismatched"
  )) +
  scale_fill_manual(values = c(
    "Weak" = "#4C7434",
    "Intermediate" = "#2F7F73",
    "Strong" = "#1F4E79"
  )) +
  coord_cartesian(ylim = c(0, 2)) +
  labs(
    x = "Intraspecific Covariance Scenario",
    y = "Correlation Mismatch Magnitude",
    fill = "Attraction Strength"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 10)),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8.5, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 7.5),
    legend.key.size = grid::unit(0.4, "cm"),
    legend.spacing.x = grid::unit(0.1, "cm"),
    panel.grid.minor = element_blank()
  )

mismatch_plot

ggsave(file.path(output_dir, "ou mismatch plot.png"),
       plot = mismatch_plot,
       width = 8,
       height = 5,
       dpi = 600)
