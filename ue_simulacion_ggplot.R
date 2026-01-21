############################################################
# UE (Universo Emergente) - SIMULACION MEJORADA EN R
# Con tidyverse y visualizaciones ggplot2 avanzadas
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
})

set.seed(42)

# === TEMA PERSONALIZADO PARA GRAFICOS ===
theme_ue <- function() {
  theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.background = element_rect(fill = "white", color = NA)
    )
}

# === 1. GENERACION DE DATOS SINTETICOS ===

generate_synthetic_data <- function(n_times = 500,
                                     window_size = 30,
                                     noise_sd = 0.3) {

  window_starts <- c(50, 200, 380)
  windows_R <- map(window_starts, ~.x:(.x + window_size - 1))

  is_R <- rep(FALSE, n_times)
  for (w in windows_R) is_R[w] <- TRUE

  # Generar señales
  signal1 <- rnorm(n_times, mean = 0, sd = noise_sd)
  signal2 <- rnorm(n_times, mean = 0, sd = noise_sd)
  signal3 <- rnorm(n_times, mean = 0, sd = noise_sd * 0.5)  # Canal adicional

  # Modificar en ventanas R

  for (w in windows_R) {
    signal1[w] <- signal1[w] + 2.0
    signal2[w] <- signal2[w] + sin(seq(0, 2*pi, length.out = length(w))) * 0.8
    signal3[w] <- signal3[w] + 0.5 * cos(seq(0, 4*pi, length.out = length(w)))
  }

  data <- tibble(
    t = 1:n_times,
    signal1 = signal1,
    signal2 = signal2,
    signal3 = signal3,
    is_R = is_R,
    region = ifelse(is_R, "R (Evento Raro)", "Normal")
  )

  list(data = data, windows_R = windows_R, n_times = n_times)
}

# === 2. ERS: EVENTO-RARO-SOSTENIDO ===

build_ers <- function(data, predicate_fn, persistence_tau = 10) {
  mask <- map_lgl(1:nrow(data$data), ~predicate_fn(data$data, .x))
  rle_result <- rle(mask)

  windows_R <- list()
  idx <- 1
  pos <- 1

  for (i in seq_along(rle_result$lengths)) {
    len <- rle_result$lengths[i]
    val <- rle_result$values[i]
    if (val && len >= persistence_tau) {
      windows_R[[idx]] <- pos:(pos + len - 1)
      idx <- idx + 1
    }
    pos <- pos + len
  }

  all_R_times <- unlist(windows_R)
  notR_times <- setdiff(1:nrow(data$data), all_R_times)

  windows_notR <- map(windows_R, function(w) {
    w_len <- length(w)
    if (length(notR_times) >= w_len) {
      start_idx <- sample(1:(length(notR_times) - w_len + 1), 1)
      notR_times[start_idx:(start_idx + w_len - 1)]
    }
  })

  list(name = "ERS", windows_R = windows_R, windows_notR = windows_notR,
       persistence_tau = persistence_tau)
}

# === 3. REPRESENTACION: EXTRAER MICROESTADOS ===

extract_microstates <- function(data, ers) {
  extract_from_windows <- function(windows, label_R, start_id) {
    window_id <- start_id
    states_list <- list()

    for (w in windows) {
      for (t_idx in w) {
        row <- data$data[t_idx, ]
        states_list[[length(states_list) + 1]] <- tibble(
          t = t_idx,
          x1 = row$signal1,
          x2 = row$signal2,
          x3 = row$signal3,
          x1_sq = row$signal1^2,
          x2_sq = row$signal2^2,
          x1x2 = row$signal1 * row$signal2,
          label_R = label_R,
          window_id = window_id
        )
      }
      window_id <- window_id + 1
    }
    bind_rows(states_list)
  }

  states_R <- extract_from_windows(ers$windows_R, TRUE, 1)
  states_notR <- extract_from_windows(ers$windows_notR, FALSE, length(ers$windows_R) + 1)

  bind_rows(states_R, states_notR)
}

# === 4. ENERGIA / COSTE ===

fit_energy_classifier <- function(states) {
  model <- glm(
    label_R ~ x1 + x2 + x3 + x1_sq + x2_sq + x1x2,
    data = states,
    family = binomial(link = "logit")
  )

  energy_fn <- function(row) {
    new_data <- tibble(
      x1 = row$x1, x2 = row$x2, x3 = row$x3,
      x1_sq = row$x1^2, x2_sq = row$x2^2, x1x2 = row$x1 * row$x2
    )
    pR <- predict(model, newdata = new_data, type = "response")
    pR <- pmax(pmin(pR, 1 - 1e-9), 1e-9)
    -log(pR / (1 - pR))
  }

  predictions <- predict(model, type = "response") > 0.5

  list(model = model, energy = energy_fn,
       accuracy = mean(predictions == states$label_R))
}

# === 5. CANDADOS (LOCKS) ===

detect_locks <- function(states, min_score = 0.1) {
  states_R <- filter(states, label_R == TRUE)
  states_notR <- filter(states, label_R == FALSE)

  # Reglas candidatas
  rules <- list(
    list(name = "x1 > Q25", pred = function(s) s$x1 > quantile(states_R$x1, 0.25)),
    list(name = "x1 > Q50", pred = function(s) s$x1 > quantile(states_R$x1, 0.50)),
    list(name = "x1_sq > Q25", pred = function(s) s$x1_sq > quantile(states_R$x1_sq, 0.25)),
    list(name = "x1 + x2 > 1.5", pred = function(s) (s$x1 + s$x2) > 1.5),
    list(name = "x1 + x2 > 2.0", pred = function(s) (s$x1 + s$x2) > 2.0),
    list(name = "|x2| < 0.5 & x1 > 1", pred = function(s) abs(s$x2) < 0.5 & s$x1 > 1)
  )

  score_rule <- function(rule) {
    eval_R <- map_lgl(1:nrow(states_R), ~rule$pred(states_R[.x, ]))
    eval_notR <- map_lgl(1:nrow(states_notR), ~rule$pred(states_notR[.x, ]))

    support_R <- mean(eval_R)
    support_notR <- mean(eval_notR)
    specificity <- log((support_R + 1e-9) / (support_notR + 1e-9))
    score <- support_R * specificity - 0.1 * (1 - support_R)

    tibble(name = rule$name, support_R = support_R, support_notR = support_notR,
           specificity = specificity, score = score)
  }

  locks <- map_dfr(rules, score_rule) %>%
    filter(score >= min_score) %>%
    arrange(desc(score))

  locks
}

# === 6. COARSE-GRAINING Y MSM ===

build_msm <- function(states_R, n_clusters = 4) {
  X <- select(states_R, x1, x2) %>% as.matrix()
  km <- kmeans(X, centers = n_clusters, nstart = 20)

  states_R <- mutate(states_R, macrostate = km$cluster)
  states_R <- arrange(states_R, t)

  seq_macro <- states_R$macrostate
  P <- matrix(0, nrow = n_clusters, ncol = n_clusters)

  for (i in 1:(length(seq_macro) - 1)) {
    P[seq_macro[i], seq_macro[i + 1]] <- P[seq_macro[i], seq_macro[i + 1]] + 1
  }

  row_sums <- rowSums(P)
  row_sums[row_sums == 0] <- 1
  P <- P / row_sums

  eigen_result <- eigen(P)
  eigenvalues <- Re(eigen_result$values)
  sorted_ev <- sort(abs(eigenvalues), decreasing = TRUE)

  spectral_gap <- sorted_ev[1] - sorted_ev[2]
  lambda2 <- sorted_ev[2]
  tau_relax <- if (lambda2 > 0 && lambda2 < 1) -1 / log(lambda2) else Inf

  list(
    P = P,
    states = states_R,
    centroids = km$centers,
    eigenvalues = eigenvalues,
    spectral_gap = spectral_gap,
    tau_relax = tau_relax,
    metastable_ok = mean(diag(P)) > 0.4
  )
}

# === 7. VISUALIZACIONES MEJORADAS ===

plot_timeseries_enhanced <- function(data, ers) {
  # Preparar datos de ventanas R
  windows_df <- tibble(
    xmin = map_dbl(ers$windows_R, min),
    xmax = map_dbl(ers$windows_R, max),
    window = paste("R", seq_along(ers$windows_R))
  )

  # Datos en formato largo para las señales
  data_long <- data$data %>%
    pivot_longer(cols = c(signal1, signal2, signal3),
                 names_to = "canal", values_to = "valor") %>%
    mutate(canal = case_when(
      canal == "signal1" ~ "Canal 1 (Principal)",
      canal == "signal2" ~ "Canal 2 (Oscilatorio)",
      canal == "signal3" ~ "Canal 3 (Secundario)"
    ))

  ggplot() +
    geom_rect(data = windows_df,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "#E74C3C", alpha = 0.15) +
    geom_line(data = data_long, aes(x = t, y = valor, color = canal),
              linewidth = 0.6, alpha = 0.8) +
    geom_text(data = windows_df, aes(x = (xmin + xmax)/2, y = 3, label = window),
              color = "#C0392B", fontface = "bold", size = 3) +
    scale_color_manual(values = c("#3498DB", "#27AE60", "#9B59B6")) +
    labs(
      title = "Serie Temporal - Modelo Universo Emergente",
      subtitle = "Ventanas R (Evento Raro Sostenido) marcadas en rojo",
      x = "Tiempo", y = "Amplitud", color = "Canal"
    ) +
    theme_ue() +
    coord_cartesian(ylim = c(-2, 3.5))
}

plot_energy_enhanced <- function(states, energy_clf) {
  states <- states %>%
    rowwise() %>%
    mutate(energy = energy_clf$energy(cur_data())) %>%
    ungroup() %>%
    mutate(condicion = ifelse(label_R, "R (Evento Raro)", "No-R (Control)"))

  p1 <- ggplot(states, aes(x = energy, fill = condicion)) +
    geom_density(alpha = 0.6, color = NA) +
    geom_rug(aes(color = condicion), alpha = 0.3) +
    scale_fill_manual(values = c("R (Evento Raro)" = "#E74C3C", "No-R (Control)" = "#3498DB")) +
    scale_color_manual(values = c("R (Evento Raro)" = "#E74C3C", "No-R (Control)" = "#3498DB")) +
    labs(
      title = "Distribución de Energía",
      subtitle = paste("E(x) = -log(odds) | Accuracy:", round(energy_clf$accuracy, 3)),
      x = "Energía E(x)", y = "Densidad"
    ) +
    theme_ue() +
    guides(color = "none")

  p2 <- ggplot(states, aes(x = condicion, y = energy, fill = condicion)) +
    geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    geom_jitter(width = 0.1, alpha = 0.2, size = 0.8) +
    scale_fill_manual(values = c("R (Evento Raro)" = "#E74C3C", "No-R (Control)" = "#3498DB")) +
    labs(x = "", y = "Energía E(x)") +
    theme_ue() +
    theme(legend.position = "none")

  grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
}

plot_msm_enhanced <- function(msm) {
  n <- nrow(msm$P)

  P_df <- as_tibble(msm$P, .name_repair = "unique") %>%
    mutate(from = factor(1:n)) %>%
    pivot_longer(cols = -from, names_to = "to_col", values_to = "prob") %>%
    mutate(to = factor(as.integer(gsub("\\.\\.\\.", "", to_col))))

  p1 <- ggplot(P_df, aes(x = to, y = fct_rev(from), fill = prob)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", prob)), color = "white",
              fontface = "bold", size = 4) +
    scale_fill_gradient2(low = "#2C3E50", mid = "#E74C3C", high = "#F39C12",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(
      title = "Matriz de Transición MSM",
      subtitle = paste("Gap espectral:", round(msm$spectral_gap, 3),
                       "| τ relajación:", round(msm$tau_relax, 2),
                       "| Metaestable:", msm$metastable_ok),
      x = "Estado destino (Φ)", y = "Estado origen (Φ)", fill = "P(trans)"
    ) +
    theme_ue() +
    coord_equal()

  # Gráfico de autovalores
  ev_df <- tibble(
    idx = 1:length(msm$eigenvalues),
    valor = abs(msm$eigenvalues)
  )

  p2 <- ggplot(ev_df, aes(x = factor(idx), y = valor)) +
    geom_col(fill = "#3498DB", alpha = 0.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = "Espectro MSM", x = "Índice", y = "|λ|") +
    theme_ue()

  grid.arrange(p1, p2, ncol = 2, widths = c(2, 1))
}

plot_locks_enhanced <- function(locks) {
  if (nrow(locks) == 0) {
    return(ggplot() +
             annotate("text", x = 0.5, y = 0.5, label = "No se detectaron candados",
                      size = 6) + theme_void())
  }

  locks_long <- locks %>%
    select(name, support_R, support_notR, score) %>%
    pivot_longer(cols = c(support_R, support_notR),
                 names_to = "tipo", values_to = "soporte") %>%
    mutate(tipo = ifelse(tipo == "support_R", "Soporte R", "Soporte No-R"))

  p1 <- ggplot(locks_long, aes(x = reorder(name, -soporte), y = soporte, fill = tipo)) +
    geom_col(position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("Soporte R" = "#E74C3C", "Soporte No-R" = "#3498DB")) +
    labs(
      title = "Candados (Locks) Detectados",
      subtitle = "Reglas que distinguen R de No-R",
      x = "", y = "Soporte (proporción)", fill = ""
    ) +
    theme_ue() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  p2 <- ggplot(locks, aes(x = reorder(name, score), y = score)) +
    geom_col(fill = "#27AE60", alpha = 0.8) +
    geom_text(aes(label = round(score, 2)), hjust = -0.2, size = 3) +
    coord_flip() +
    labs(title = "Score de Candados", x = "", y = "Score") +
    theme_ue()

  grid.arrange(p1, p2, ncol = 2)
}

plot_macrostates_enhanced <- function(msm) {
  states <- msm$states %>%
    mutate(macrostate = factor(macrostate))

  centroids <- as_tibble(msm$centroids) %>%
    mutate(macrostate = factor(1:n()))

  ggplot() +
    stat_ellipse(data = states, aes(x = x1, y = x2, color = macrostate),
                 level = 0.9, linetype = "dashed", linewidth = 0.8) +
    geom_point(data = states, aes(x = x1, y = x2, color = macrostate),
               alpha = 0.6, size = 2) +
    geom_point(data = centroids, aes(x = x1, y = x2),
               shape = 4, size = 5, stroke = 2, color = "black") +
    geom_label(data = centroids, aes(x = x1, y = x2, label = paste("Φ", macrostate)),
               nudge_y = 0.3, size = 3, fontface = "bold") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = "Macroestados (Coarse-Graining)",
      subtitle = paste("Clustering k-means |", nrow(centroids), "macroestados | τ =",
                       round(msm$tau_relax, 2)),
      x = "x1 (Canal principal)", y = "x2 (Canal oscilatorio)",
      color = "Macroestado Φ"
    ) +
    theme_ue()
}

plot_phase_space <- function(states) {
  states <- mutate(states, condicion = ifelse(label_R, "R", "No-R"))

  ggplot(states, aes(x = x1, y = x2, color = condicion)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_density_2d(linewidth = 0.5) +
    scale_color_manual(values = c("R" = "#E74C3C", "No-R" = "#3498DB")) +
    facet_wrap(~condicion) +
    labs(
      title = "Espacio de Fases",
      subtitle = "Distribución de microestados en el plano (x1, x2)",
      x = "x1", y = "x2"
    ) +
    theme_ue() +
    theme(legend.position = "none")
}

# === 8. PIPELINE PRINCIPAL ===

run_ue_simulation <- function() {
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║     SIMULACIÓN UNIVERSO EMERGENTE (UE) - VERSIÓN MEJORADA   ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n\n")

  # 1. Generar datos
  cat("▶ 1. Generando datos sintéticos...\n")
  data <- generate_synthetic_data(n_times = 500)
  cat("    ├─ Tiempos:", data$n_times, "\n")
  cat("    └─ Ventanas R:", length(data$windows_R), "\n\n")

  # 2. Construir ERS
  cat("▶ 2. Construyendo ERS (Evento-Raro-Sostenido)...\n")
  predicate <- function(df, t) df$is_R[t]
  ers <- build_ers(data, predicate, persistence_tau = 10)
  cat("    ├─ Ventanas R detectadas:", length(ers$windows_R), "\n")
  cat("    └─ Ventanas control:", length(ers$windows_notR), "\n\n")

  # 3. Extraer microestados
  cat("▶ 3. Extrayendo microestados...\n")
  states <- extract_microstates(data, ers)
  cat("    ├─ Total microestados:", nrow(states), "\n")
  cat("    ├─ Estados R:", sum(states$label_R), "\n")
  cat("    └─ Estados No-R:", sum(!states$label_R), "\n\n")

  # 4. Clasificador de energía
  cat("▶ 4. Ajustando clasificador de energía...\n")
  energy_clf <- fit_energy_classifier(states)
  cat("    └─ Accuracy:", round(energy_clf$accuracy, 3), "\n\n")

  # 5. Detectar candados
  cat("▶ 5. Detectando candados (locks)...\n")
  locks <- detect_locks(states, min_score = 0.1)
  cat("    └─ Candados detectados:", nrow(locks), "\n")
  if (nrow(locks) > 0) {
    for (i in 1:min(3, nrow(locks))) {
      cat("       •", locks$name[i], "| Score:", round(locks$score[i], 2), "\n")
    }
  }
  cat("\n")

  # 6. MSM
  cat("▶ 6. Construyendo macroestados y MSM...\n")
  states_R <- filter(states, label_R == TRUE)
  msm <- build_msm(states_R, n_clusters = 4)
  cat("    ├─ Macroestados:", length(unique(msm$states$macrostate)), "\n")
  cat("    ├─ Gap espectral:", round(msm$spectral_gap, 3), "\n")
  cat("    ├─ τ relajación:", round(msm$tau_relax, 2), "\n")
  cat("    └─ Metaestable:", msm$metastable_ok, "\n\n")

  # 7. Generar visualizaciones
  cat("▶ 7. Generando visualizaciones...\n")

  ggsave("plot_01_timeseries.png", plot_timeseries_enhanced(data, ers),
         width = 12, height = 5, dpi = 150)
  cat("    ├─ plot_01_timeseries.png\n")

  png("plot_02_energy.png", width = 1200, height = 500, res = 150)
  plot_energy_enhanced(states, energy_clf)
  dev.off()
  cat("    ├─ plot_02_energy.png\n")

  png("plot_03_msm.png", width = 1200, height = 500, res = 150)
  plot_msm_enhanced(msm)
  dev.off()
  cat("    ├─ plot_03_msm.png\n")

  png("plot_04_locks.png", width = 1200, height = 500, res = 150)
  plot_locks_enhanced(locks)
  dev.off()
  cat("    ├─ plot_04_locks.png\n")

  ggsave("plot_05_macrostates.png", plot_macrostates_enhanced(msm),
         width = 10, height = 7, dpi = 150)
  cat("    ├─ plot_05_macrostates.png\n")

  ggsave("plot_06_phase_space.png", plot_phase_space(states),
         width = 10, height = 5, dpi = 150)
  cat("    └─ plot_06_phase_space.png\n\n")

  # 8. Resumen
  status <- ifelse(nrow(locks) > 0 && msm$metastable_ok, "PASS", "FAIL")

  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║                    RESUMEN DOMINIO UE                        ║\n")
  cat("╠══════════════════════════════════════════════════════════════╣\n")
  cat(sprintf("║  Nombre: %-51s ║\n", "D_ERS"))
  cat(sprintf("║  Status: %-51s ║\n", status))
  cat(sprintf("║  Candados encontrados: %-38s ║\n", nrow(locks)))
  cat(sprintf("║  Metaestabilidad: %-43s ║\n", msm$metastable_ok))
  cat(sprintf("║  Accuracy clasificador: %-37s ║\n", round(energy_clf$accuracy, 3)))
  cat("╚══════════════════════════════════════════════════════════════╝\n")

  invisible(list(
    data = data, ers = ers, states = states,
    energy_classifier = energy_clf, locks = locks, msm = msm, status = status
  ))
}

# === EJECUTAR ===
cat("\n🚀 Iniciando simulación UE mejorada...\n")
results <- run_ue_simulation()
cat("\n✅ Simulación completada. Gráficos guardados.\n\n")
