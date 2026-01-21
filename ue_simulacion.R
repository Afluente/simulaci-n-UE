############################################################
# UE (Universo Emergente) - SIMULACION EN R
# Demo minima funcional con datos sinteticos
# Version: Base R (sin dependencias externas)
############################################################

set.seed(42)

# === 1. GENERACION DE DATOS SINTETICOS ===

generate_synthetic_data <- function(n_times = 500,
                                     n_windows_R = 3,
                                     window_size = 30,
                                     noise_sd = 0.3) {

  times <- 1:n_times

  # Inicializar canales
  signal1 <- rnorm(n_times, mean = 0, sd = noise_sd)
  signal2 <- rnorm(n_times, mean = 0, sd = noise_sd)

  # Definir ventanas R (posiciones fijas para reproducibilidad)
  window_starts <- c(50, 200, 380)
  windows_R <- lapply(window_starts, function(s) s:(s + window_size - 1))

  # Crear mascara R
  is_R <- rep(FALSE, n_times)
  for (w in windows_R) {
    is_R[w] <- TRUE
  }

  # Modificar senales en ventanas R
  for (w in windows_R) {
    signal1[w] <- signal1[w] + 2.0  # shift de media
    signal2[w] <- signal2[w] + sin(seq(0, 2*pi, length.out = length(w))) * 0.8
  }

  # Crear dataframe
  data <- data.frame(
    t = times,
    signal1 = signal1,
    signal2 = signal2,
    is_R = is_R
  )

  list(
    data = data,
    windows_R = windows_R,
    n_times = n_times
  )
}

# === 2. ERS: EVENTO-RARO-SOSTENIDO ===

build_ers <- function(data, predicate_fn, persistence_tau = 10) {
  # Aplica predicado a cada tiempo
  mask <- sapply(1:nrow(data$data), function(t) predicate_fn(data$data, t))

  # Encontrar segmentos consecutivos TRUE
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

  # Crear ventanas de control (notR) emparejadas
  all_R_times <- unlist(windows_R)
  notR_times <- setdiff(1:nrow(data$data), all_R_times)

  windows_notR <- list()
  for (i in seq_along(windows_R)) {
    w_len <- length(windows_R[[i]])
    if (length(notR_times) >= w_len) {
      start_idx <- sample(1:(length(notR_times) - w_len + 1), 1)
      windows_notR[[i]] <- notR_times[start_idx:(start_idx + w_len - 1)]
    }
  }

  list(
    name = "ERS",
    windows_R = windows_R,
    windows_notR = windows_notR,
    persistence_tau = persistence_tau
  )
}

# === 3. REPRESENTACION: EXTRAER MICROESTADOS ===

represent <- function(data_row) {
  c(
    x1 = data_row$signal1,
    x2 = data_row$signal2,
    x1_sq = data_row$signal1^2,
    x2_sq = data_row$signal2^2,
    x1x2 = data_row$signal1 * data_row$signal2
  )
}

extract_microstates <- function(data, ers) {
  states <- data.frame(
    t = integer(),
    x1 = numeric(),
    x2 = numeric(),
    x1_sq = numeric(),
    x2_sq = numeric(),
    x1x2 = numeric(),
    label_R = logical(),
    window_id = integer()
  )

  window_id <- 1

  # Extraer estados de ventanas R
  for (w in ers$windows_R) {
    for (t_idx in w) {
      features <- represent(data$data[t_idx, ])
      states <- rbind(states, data.frame(
        t = t_idx,
        x1 = features["x1"],
        x2 = features["x2"],
        x1_sq = features["x1_sq"],
        x2_sq = features["x2_sq"],
        x1x2 = features["x1x2"],
        label_R = TRUE,
        window_id = window_id
      ))
    }
    window_id <- window_id + 1
  }

  # Extraer estados de ventanas notR
  for (w in ers$windows_notR) {
    for (t_idx in w) {
      features <- represent(data$data[t_idx, ])
      states <- rbind(states, data.frame(
        t = t_idx,
        x1 = features["x1"],
        x2 = features["x2"],
        x1_sq = features["x1_sq"],
        x2_sq = features["x2_sq"],
        x1x2 = features["x1x2"],
        label_R = FALSE,
        window_id = window_id
      ))
    }
    window_id <- window_id + 1
  }

  rownames(states) <- NULL
  states
}

# === 4. ENERGIA / COSTE ===

fit_energy_classifier <- function(states) {
  # Clasificador logistico para distinguir R vs notR
  model <- glm(
    label_R ~ x1 + x2 + x1_sq + x2_sq + x1x2,
    data = states,
    family = binomial(link = "logit")
  )

  # Funcion de energia: E(x) = -log(odds)
  energy_fn <- function(x_vec) {
    new_data <- data.frame(
      x1 = x_vec[1],
      x2 = x_vec[2],
      x1_sq = x_vec[1]^2,
      x2_sq = x_vec[2]^2,
      x1x2 = x_vec[1] * x_vec[2]
    )
    pR <- predict(model, newdata = new_data, type = "response")
    pR <- pmax(pmin(pR, 1 - 1e-9), 1e-9)  # clamp
    odds <- pR / (1 - pR)
    -log(odds)
  }

  predictions <- predict(model, type = "response") > 0.5
  accuracy <- mean(predictions == states$label_R)

  list(
    model = model,
    energy = energy_fn,
    accuracy = accuracy
  )
}

build_cost <- function(energy_fn, epsilon = 1.0) {
  function(x_vec) {
    epsilon * energy_fn(x_vec)
  }
}

# === 5. CANDADOS (LOCKS) ===

propose_lock_hypotheses <- function(states_R, states_notR) {
  rules <- list()

  # Regla 1: x1 > umbral
  threshold_x1 <- quantile(states_R$x1, 0.25)
  rules[[1]] <- list(
    name = "x1 > threshold",
    predicate = function(x) x[1] > threshold_x1,
    threshold = threshold_x1
  )

  # Regla 2: x1^2 > umbral
  threshold_x1sq <- quantile(states_R$x1_sq, 0.25)
  rules[[2]] <- list(
    name = "x1_sq > threshold",
    predicate = function(x) x[3] > threshold_x1sq,
    threshold = threshold_x1sq
  )

  # Regla 3: combinacion lineal
  rules[[3]] <- list(
    name = "x1 + x2 > 1.5",
    predicate = function(x) (x[1] + x[2]) > 1.5,
    threshold = 1.5
  )

  rules
}

score_lock <- function(rule, states_R, states_notR) {
  eval_R <- sapply(1:nrow(states_R), function(i) {
    rule$predicate(c(states_R$x1[i], states_R$x2[i], states_R$x1_sq[i]))
  })

  eval_notR <- sapply(1:nrow(states_notR), function(i) {
    rule$predicate(c(states_notR$x1[i], states_notR$x2[i], states_notR$x1_sq[i]))
  })

  support_R <- mean(eval_R)
  support_notR <- mean(eval_notR)

  specificity <- log((support_R + 1e-9) / (support_notR + 1e-9))
  tau_persist <- sum(eval_R)
  score <- support_R * specificity - 0.1 * (1 - support_R)

  list(
    name = rule$name,
    support_R = support_R,
    support_notR = support_notR,
    specificity = specificity,
    tau_persist = tau_persist,
    score = score
  )
}

detect_locks <- function(states, ers, min_score = 0.5) {
  states_R <- states[states$label_R == TRUE, ]
  states_notR <- states[states$label_R == FALSE, ]

  rules <- propose_lock_hypotheses(states_R, states_notR)
  locks <- lapply(rules, function(r) score_lock(r, states_R, states_notR))

  # Filtrar por score minimo
  locks <- Filter(function(l) l$score >= min_score, locks)

  # Ordenar por score descendente
  scores <- sapply(locks, function(l) l$score)
  locks <- locks[order(scores, decreasing = TRUE)]

  locks
}

# === 6. COARSE-GRAINING Y MSM ===

build_macrostates <- function(states_R, n_clusters = 3) {
  X <- as.matrix(states_R[, c("x1", "x2")])
  km <- kmeans(X, centers = n_clusters, nstart = 10)

  states_R$macrostate <- km$cluster

  list(
    states = states_R,
    centroids = km$centers,
    n_clusters = n_clusters
  )
}

fit_msm <- function(macro_result, lag = 1) {
  states <- macro_result$states[order(macro_result$states$t), ]
  n_clusters <- macro_result$n_clusters

  seq_macro <- states$macrostate

  # Estimar matriz de transicion
  P <- matrix(0, nrow = n_clusters, ncol = n_clusters)

  for (i in 1:(length(seq_macro) - lag)) {
    from <- seq_macro[i]
    to <- seq_macro[i + lag]
    P[from, to] <- P[from, to] + 1
  }

  # Normalizar filas
  row_sums <- rowSums(P)
  row_sums[row_sums == 0] <- 1
  P <- P / row_sums

  # Calcular autovalores
  eigen_result <- eigen(P)
  eigenvalues <- Re(eigen_result$values)

  sorted_ev <- sort(abs(eigenvalues), decreasing = TRUE)
  spectral_gap <- sorted_ev[1] - sorted_ev[2]

  lambda2 <- sorted_ev[2]
  tau_relax <- if (lambda2 > 0 && lambda2 < 1) -lag / log(lambda2) else Inf

  diag_mean <- mean(diag(P))
  metastable_ok <- diag_mean > 0.5

  list(
    P = P,
    eigenvalues = eigenvalues,
    spectral_gap = spectral_gap,
    tau_relax = tau_relax,
    metastable_ok = metastable_ok,
    states = macro_result$states,
    centroids = macro_result$centroids
  )
}

# === 7. VISUALIZACIONES (Base R) ===

plot_timeseries <- function(data, ers, filename = "plot_timeseries.png") {
  png(filename, width = 1000, height = 400)

  par(mar = c(4, 4, 3, 1))
  plot(data$data$t, data$data$signal1, type = "l", col = "blue",
       xlab = "Tiempo", ylab = "Valor",
       main = "Serie Temporal con Ventanas R (Evento Raro Sostenido)",
       ylim = range(c(data$data$signal1, data$data$signal2)))
  lines(data$data$t, data$data$signal2, col = "darkgreen")

  # Marcar ventanas R
  for (w in ers$windows_R) {
    rect(min(w), par("usr")[3], max(w), par("usr")[4],
         col = rgb(1, 0, 0, 0.2), border = NA)
  }

  legend("topright", legend = c("Signal 1", "Signal 2", "Ventana R"),
         col = c("blue", "darkgreen", rgb(1, 0, 0, 0.3)),
         lty = c(1, 1, NA), pch = c(NA, NA, 15), bg = "white")

  dev.off()
  cat("  Guardado:", filename, "\n")
}

plot_energy_distribution <- function(states, energy_clf, filename = "plot_energy.png") {
  # Calcular energia para cada estado
  energies <- sapply(1:nrow(states), function(i) {
    energy_clf$energy(c(states$x1[i], states$x2[i]))
  })

  states$energy <- energies

  png(filename, width = 800, height = 500)

  par(mar = c(4, 4, 3, 1))

  energy_R <- states$energy[states$label_R]
  energy_notR <- states$energy[!states$label_R]

  hist(energy_notR, col = rgb(0, 0, 1, 0.5), border = "blue",
       main = paste("Distribucion de Energia: R vs No-R\nAccuracy:",
                    round(energy_clf$accuracy, 3)),
       xlab = "Energia E(x)", ylab = "Frecuencia",
       xlim = range(c(energy_R, energy_notR)),
       breaks = 20)
  hist(energy_R, col = rgb(1, 0, 0, 0.5), border = "red", add = TRUE, breaks = 20)

  legend("topright", legend = c("R", "No-R"),
         fill = c(rgb(1, 0, 0, 0.5), rgb(0, 0, 1, 0.5)), bg = "white")

  dev.off()
  cat("  Guardado:", filename, "\n")
}

plot_msm_transitions <- function(msm, filename = "plot_msm.png") {
  png(filename, width = 600, height = 500)

  par(mar = c(4, 4, 3, 1))

  n <- nrow(msm$P)
  image(1:n, 1:n, t(msm$P)[, n:1], col = heat.colors(20),
        xlab = "Estado destino", ylab = "Estado origen",
        main = paste("Matriz de Transicion MSM\nGap espectral:",
                     round(msm$spectral_gap, 3),
                     "| Metaestable:", msm$metastable_ok),
        axes = FALSE)
  axis(1, at = 1:n)
  axis(2, at = 1:n, labels = n:1)

  # Agregar valores
  for (i in 1:n) {
    for (j in 1:n) {
      text(j, n - i + 1, round(msm$P[i, j], 2), cex = 1.2)
    }
  }

  dev.off()
  cat("  Guardado:", filename, "\n")
}

plot_locks <- function(locks, filename = "plot_locks.png") {
  png(filename, width = 800, height = 500)

  if (length(locks) == 0) {
    plot.new()
    text(0.5, 0.5, "No se detectaron candados", cex = 1.5)
    dev.off()
    return()
  }

  par(mar = c(8, 4, 3, 1))

  names_locks <- sapply(locks, function(l) l$name)
  support_R <- sapply(locks, function(l) l$support_R)
  support_notR <- sapply(locks, function(l) l$support_notR)

  n <- length(locks)
  x <- barplot(rbind(support_R, support_notR),
               beside = TRUE,
               col = c("red", "blue"),
               names.arg = rep("", n),
               main = "Candados Detectados\nComparacion de soporte R vs No-R",
               ylab = "Soporte (proporcion)",
               ylim = c(0, 1))

  text(colMeans(x), par("usr")[3] - 0.05, labels = names_locks,
       srt = 45, adj = 1, xpd = TRUE, cex = 0.9)

  legend("topright", legend = c("Soporte R", "Soporte No-R"),
         fill = c("red", "blue"), bg = "white")

  dev.off()
  cat("  Guardado:", filename, "\n")
}

plot_macrostates <- function(msm, filename = "plot_macrostates.png") {
  png(filename, width = 800, height = 600)

  par(mar = c(4, 4, 3, 1))

  colors <- c("red", "blue", "green", "purple", "orange")[1:max(msm$states$macrostate)]

  plot(msm$states$x1, msm$states$x2,
       col = colors[msm$states$macrostate],
       pch = 19, cex = 1.2,
       xlab = "x1", ylab = "x2",
       main = paste("Macroestados (Coarse-Graining)\nTau relajacion:",
                    round(msm$tau_relax, 2)))

  # Marcar centroides
  points(msm$centroids[, 1], msm$centroids[, 2],
         pch = 4, cex = 3, lwd = 3, col = "black")

  legend("topright",
         legend = c(paste("Macroestado", 1:nrow(msm$centroids)), "Centroide"),
         col = c(colors[1:nrow(msm$centroids)], "black"),
         pch = c(rep(19, nrow(msm$centroids)), 4),
         pt.cex = c(rep(1.2, nrow(msm$centroids)), 2),
         bg = "white")

  dev.off()
  cat("  Guardado:", filename, "\n")
}

# === 8. PIPELINE PRINCIPAL ===

run_ue_simulation <- function() {
  cat("========================================\n")
  cat("   SIMULACION UNIVERSO EMERGENTE (UE)  \n")
  cat("========================================\n\n")

  # 1. Generar datos
  cat("1. Generando datos sinteticos...\n")
  data <- generate_synthetic_data(n_times = 500)
  cat("   - Tiempos:", data$n_times, "\n")
  cat("   - Ventanas R:", length(data$windows_R), "\n\n")

  # 2. Construir ERS
  cat("2. Construyendo ERS (Evento-Raro-Sostenido)...\n")
  predicate <- function(df, t) df$is_R[t]
  ers <- build_ers(data, predicate, persistence_tau = 10)
  cat("   - Ventanas R detectadas:", length(ers$windows_R), "\n")
  cat("   - Ventanas control:", length(ers$windows_notR), "\n\n")

  # 3. Extraer microestados
  cat("3. Extrayendo microestados...\n")
  states <- extract_microstates(data, ers)
  cat("   - Total microestados:", nrow(states), "\n")
  cat("   - Estados R:", sum(states$label_R), "\n")
  cat("   - Estados No-R:", sum(!states$label_R), "\n\n")

  # 4. Ajustar clasificador de energia
  cat("4. Ajustando clasificador de energia...\n")
  energy_clf <- fit_energy_classifier(states)
  cat("   - Accuracy:", round(energy_clf$accuracy, 3), "\n\n")

  # 5. Detectar candados
  cat("5. Detectando candados (locks)...\n")
  locks <- detect_locks(states, ers, min_score = 0.1)
  cat("   - Candados detectados:", length(locks), "\n")
  for (l in locks) {
    cat("     *", l$name, "| Score:", round(l$score, 3),
        "| Especificidad:", round(l$specificity, 3), "\n")
  }
  cat("\n")

  # 6. Coarse-graining y MSM
  cat("6. Construyendo macroestados y MSM...\n")
  states_R <- states[states$label_R == TRUE, ]
  macro <- build_macrostates(states_R, n_clusters = 3)
  msm <- fit_msm(macro, lag = 1)
  cat("   - Macroestados:", length(unique(msm$states$macrostate)), "\n")
  cat("   - Gap espectral:", round(msm$spectral_gap, 3), "\n")
  cat("   - Tau relajacion:", round(msm$tau_relax, 2), "\n")
  cat("   - Metaestable:", msm$metastable_ok, "\n\n")

  # 7. Generar visualizaciones
  cat("7. Generando visualizaciones...\n")
  plot_timeseries(data, ers)
  plot_energy_distribution(states, energy_clf)
  plot_msm_transitions(msm)
  plot_locks(locks)
  plot_macrostates(msm)
  cat("\n")

  # 8. Resumen
  cat("========================================\n")
  cat("           RESUMEN DOMINIO UE          \n")
  cat("========================================\n")
  cat("Nombre: D_ERS\n")
  status <- ifelse(length(locks) > 0 && msm$metastable_ok, "PASS", "FAIL")
  cat("Status:", status, "\n")
  cat("Candados encontrados:", length(locks), "\n")
  cat("Metaestabilidad:", msm$metastable_ok, "\n")
  cat("========================================\n")

  # Retornar resultados
  list(
    data = data,
    ers = ers,
    states = states,
    energy_classifier = energy_clf,
    locks = locks,
    msm = msm,
    status = status
  )
}

# === EJECUTAR SIMULACION ===
cat("\nIniciando simulacion UE...\n\n")
results <- run_ue_simulation()

cat("\nSimulacion completada. Resultados en 'results'.\n")
cat("Graficos guardados en el directorio actual.\n")
