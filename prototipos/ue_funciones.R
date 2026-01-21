# =============================================================================
# UE Framework - Funciones Modulares
# =============================================================================
# Cada concepto del framework es una funcion independiente y composable.
#
# Uso:
#   source("ue_funciones.R")
#   resultado <- ue_pipeline(mi_potencial, mi_observable, params)
# =============================================================================

library(stats)

# =============================================================================
# 1. DINAMICA Y SIMULACION
# =============================================================================

#' Crear funcion de paso Langevin
#' @param grad_V funcion gradiente del potencial
#' @param dt paso temporal
#' @param sigma nivel de ruido
#' @return funcion que ejecuta un paso
ue_crear_dinamica <- function(grad_V, dt = 0.01, sigma = 0.4) {
  function(x) {
    x - grad_V(x) * dt + sigma * sqrt(dt) * rnorm(length(x))
  }
}

#' Simular una trayectoria
#' @param x0 condicion inicial
#' @param n_pasos longitud de trayectoria
#' @param step_fn funcion de paso (de ue_crear_dinamica)
#' @param limites vector c(min, max) para limitar espacio
#' @return vector con la trayectoria
ue_simular_trayectoria <- function(x0, n_pasos, step_fn, limites = c(-3, 3)) {
  x <- numeric(n_pasos)
  x[1] <- x0
  for (t in 2:n_pasos) {
    x[t] <- step_fn(x[t-1])
    x[t] <- pmax(limites[1], pmin(limites[2], x[t]))
  }
  x
}

#' Simular multiples trayectorias
#' @param n_traj numero de trayectorias
#' @param n_pasos pasos por trayectoria
#' @param step_fn funcion de paso
#' @param x0_sampler funcion que genera condiciones iniciales
#' @param limites limites del espacio
#' @param seed semilla aleatoria
#' @return lista de trayectorias
ue_simular_ensemble <- function(n_traj, n_pasos, step_fn,
                                 x0_sampler = function() runif(1, -2, 2),
                                 limites = c(-3, 3), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  lapply(1:n_traj, function(i) {
    ue_simular_trayectoria(x0_sampler(), n_pasos, step_fn, limites)
  })
}

# =============================================================================
# 2. EVENTO RARO SOSTENIDO (ERS)
# =============================================================================

#' Verificar si una trayectoria cumple la condicion ERS
#' @param traj vector de la trayectoria
#' @param observable funcion θ(x)
#' @param theta_star umbral
#' @param tau_star persistencia minima
#' @param above TRUE si θ >= θ*, FALSE si θ <= θ*
#' @return TRUE/FALSE
ue_cumple_ERS <- function(traj, observable, theta_star, tau_star, above = TRUE) {
  theta_t <- observable(traj)

  if (above) {
    condicion <- theta_t >= theta_star
  } else {
    condicion <- theta_t <= theta_star
  }

  if (!any(condicion)) return(FALSE)

  rachas <- rle(condicion)
  duraciones <- rachas$lengths[rachas$values == TRUE]

  if (length(duraciones) == 0) return(FALSE)

  max(duraciones) >= tau_star
}

#' Filtrar trayectorias que cumplen R
#' @param trayectorias lista de trayectorias
#' @param observable funcion θ(x)
#' @param theta_star umbral
#' @param tau_star persistencia minima
#' @param above direccion de la condicion
#' @return lista con indices y trayectorias condicionadas
ue_filtrar_ERS <- function(trayectorias, observable, theta_star, tau_star, above = TRUE) {
  cumple <- sapply(trayectorias, ue_cumple_ERS,
                   observable = observable,
                   theta_star = theta_star,
                   tau_star = tau_star,
                   above = above)

  list(
    indices = which(cumple),
    trayectorias = trayectorias[cumple],
    n_total = length(trayectorias),
    n_condicionadas = sum(cumple),
    tasa = mean(cumple)
  )
}

# =============================================================================
# 3. PROBABILIDAD CONDICIONADA
# =============================================================================

#' Calcular distribucion P(x|R) discretizada
#' @param trayectorias_R lista de trayectorias condicionadas
#' @param n_bins numero de bins
#' @param rango vector c(min, max)
#' @return lista con mids, frecuencias, breaks
ue_distribucion_condicionada <- function(trayectorias_R, n_bins = 100, rango = c(-3, 3)) {
  x_cond <- unlist(trayectorias_R)

  breaks <- seq(rango[1], rango[2], length.out = n_bins + 1)
  mids <- (breaks[-1] + breaks[-length(breaks)]) / 2

  h <- hist(x_cond, breaks = breaks, plot = FALSE)
  freq <- h$counts / sum(h$counts)

  list(
    x = mids,
    prob = freq,
    breaks = breaks,
    n_muestras = length(x_cond),
    media = mean(x_cond),
    varianza = var(x_cond)
  )
}

# =============================================================================
# 4. PAISAJE DE COSTO F_R
# =============================================================================

#' Construir paisaje de costo F_R = -epsilon * log P(x|R)
#' @param distribucion resultado de ue_distribucion_condicionada
#' @param epsilon escala (NULL para auto: 1/log(n))
#' @return lista con x, F_R, epsilon
ue_paisaje_costo <- function(distribucion, epsilon = NULL) {
  prob <- distribucion$prob

  # Evitar log(0)
  prob[prob == 0] <- min(prob[prob > 0]) / 100

  # Epsilon automatico si no se especifica
  if (is.null(epsilon)) {
    epsilon <- 1 / log(distribucion$n_muestras)
  }

  F_R <- -epsilon * log(prob)
  F_R <- F_R - min(F_R)  # Normalizar minimo a 0

  list(
    x = distribucion$x,
    F_R = F_R,
    epsilon = epsilon,
    prob = prob
  )
}

# =============================================================================
# 5. CLASES METAESTABLES
# =============================================================================

#' Encontrar minimos locales del paisaje
#' @param paisaje resultado de ue_paisaje_costo
#' @param min_separacion separacion minima entre minimos
#' @return vector de indices de minimos
ue_encontrar_minimos <- function(paisaje, min_separacion = 0.3) {
  F_R <- paisaje$F_R
  x <- paisaje$x
  n <- length(F_R)

  # Buscar minimos locales
  minima <- c()
  for (i in 3:(n-2)) {
    if (F_R[i] < F_R[i-1] && F_R[i] < F_R[i+1] &&
        F_R[i] < F_R[i-2] && F_R[i] < F_R[i+2]) {
      minima <- c(minima, i)
    }
  }

  # Filtrar minimos muy cercanos
  if (length(minima) > 1) {
    filtered <- minima[1]
    for (m in minima[-1]) {
      if (abs(x[m] - x[filtered[length(filtered)]]) > min_separacion) {
        filtered <- c(filtered, m)
      }
    }
    minima <- filtered
  }

  minima
}

#' Asignar estados a clases (cuencas)
#' @param x_valores valores de x a clasificar
#' @param centros_clase centros de las clases (x de minimos)
#' @return vector de asignaciones (1, 2, 3, ...)
ue_asignar_clases <- function(x_valores, centros_clase) {
  sapply(x_valores, function(x) which.min(abs(x - centros_clase)))
}

#' Detectar clases metaestables completo
#' @param paisaje resultado de ue_paisaje_costo
#' @param trayectorias_R trayectorias condicionadas
#' @param min_separacion separacion minima entre clases
#' @return lista con centros, asignaciones, pesos
ue_detectar_clases <- function(paisaje, trayectorias_R, min_separacion = 0.3) {
  minimos_idx <- ue_encontrar_minimos(paisaje, min_separacion)
  centros <- paisaje$x[minimos_idx]
  n_clases <- length(centros)

  x_cond <- unlist(trayectorias_R)
  asignaciones <- ue_asignar_clases(x_cond, centros)

  # Calcular pesos P(A_k | R)
  pesos <- as.numeric(table(factor(asignaciones, levels = 1:n_clases))) / length(asignaciones)

  list(
    n_clases = n_clases,
    centros = centros,
    minimos_idx = minimos_idx,
    F_R_minimos = paisaje$F_R[minimos_idx],
    pesos = pesos,
    asignaciones = asignaciones
  )
}

# =============================================================================
# 6. METRICAS DE DOMINANCIA
# =============================================================================

#' Calcular metricas de dominancia
#' @param clases resultado de ue_detectar_clases
#' @return lista con D_H, Delta, S_eff, interpretacion
ue_metricas_dominancia <- function(clases) {
  pesos <- sort(clases$pesos, decreasing = TRUE)
  n <- length(pesos)

  if (n < 2) {
    return(list(
      D_H = Inf,
      Delta = 1,
      S_eff = 0,
      S_max = 0,
      pesos_ordenados = pesos,
      interpretacion = "Una sola clase detectada"
    ))
  }

  p1 <- pesos[1]
  p2 <- pesos[2]

  D_H <- ifelse(p2 > 0, p1 / p2, Inf)
  Delta <- p1 - p2

  # Entropia efectiva
  pesos_pos <- pesos[pesos > 0]
  S_eff <- -sum(pesos_pos * log(pesos_pos))
  S_max <- log(n)

  # Interpretacion
  if (D_H > 10) {
    interpretacion <- "Dominancia FUERTE"
  } else if (D_H > 3) {
    interpretacion <- "Dominancia MODERADA"
  } else if (D_H > 1.5) {
    interpretacion <- "Dominancia DEBIL"
  } else {
    interpretacion <- "EQUILIBRIO entre clases"
  }

  list(
    D_H = D_H,
    Delta = Delta,
    S_eff = S_eff,
    S_max = S_max,
    S_norm = S_eff / S_max,
    pesos_ordenados = pesos,
    interpretacion = interpretacion
  )
}

# =============================================================================
# 7. VALIDACION
# =============================================================================

#' Test de robustez: variar parametros ±20%
#' @param trayectorias todas las trayectorias (no filtradas)
#' @param observable funcion observable
#' @param theta_base umbral base
#' @param tau_base persistencia base
#' @param above direccion condicion
#' @return lista con resultados y veredicto
ue_test_robustez <- function(trayectorias, observable, theta_base, tau_base, above = TRUE) {
  variaciones <- expand.grid(
    theta_mult = c(0.8, 1.0, 1.2),
    tau_mult = c(0.8, 1.0, 1.2)
  )

  resultados <- lapply(1:nrow(variaciones), function(i) {
    theta_var <- theta_base * variaciones$theta_mult[i]
    tau_var <- round(tau_base * variaciones$tau_mult[i])

    filtrado <- ue_filtrar_ERS(trayectorias, observable, theta_var, tau_var, above)

    if (filtrado$n_condicionadas < 50) {
      return(list(theta = theta_var, tau = tau_var, D_H = NA, clase_dom = NA))
    }

    dist <- ue_distribucion_condicionada(filtrado$trayectorias)
    paisaje <- ue_paisaje_costo(dist)
    clases <- ue_detectar_clases(paisaje, filtrado$trayectorias)
    metricas <- ue_metricas_dominancia(clases)

    clase_dom <- which.max(clases$pesos)

    list(theta = theta_var, tau = tau_var, D_H = metricas$D_H, clase_dom = clase_dom)
  })

  D_H_vals <- sapply(resultados, function(r) r$D_H)
  clases_dom <- sapply(resultados, function(r) r$clase_dom)

  # Verificar consistencia
  clases_consistentes <- length(unique(na.omit(clases_dom))) == 1
  D_H_estable <- sd(D_H_vals, na.rm = TRUE) / mean(D_H_vals, na.rm = TRUE) < 0.3

  list(
    resultados = resultados,
    pass = clases_consistentes && D_H_estable,
    mensaje = ifelse(clases_consistentes && D_H_estable,
                     "Robustez verificada", "Robustez NO verificada")
  )
}

#' Test de ablacion: permutar estructura temporal
#' @param trayectorias_R trayectorias condicionadas
#' @param observable funcion observable
#' @param theta_star umbral
#' @param tau_star persistencia
#' @param above direccion
#' @return lista con resultados y veredicto
ue_test_ablacion <- function(trayectorias_R, observable, theta_star, tau_star, above = TRUE) {
  n_orig <- length(trayectorias_R)

  # Permutar cada trayectoria
  traj_permutadas <- lapply(trayectorias_R, sample)

  # Contar cuantas aun cumplen R
  cumple_perm <- sapply(traj_permutadas, ue_cumple_ERS,
                        observable = observable,
                        theta_star = theta_star,
                        tau_star = tau_star,
                        above = above)

  n_cumple_perm <- sum(cumple_perm)
  ratio <- n_cumple_perm / n_orig

  list(
    n_original = n_orig,
    n_permutadas_cumplen = n_cumple_perm,
    ratio = ratio,
    pass = ratio < 0.3,
    mensaje = ifelse(ratio < 0.3,
                     "Ablacion verificada: estructura temporal necesaria",
                     "Ablacion parcial: revisar condicion R")
  )
}

#' Test out-of-sample: train/test split
#' @param trayectorias todas las trayectorias
#' @param observable funcion observable
#' @param theta_star umbral
#' @param tau_star persistencia
#' @param above direccion
#' @param train_ratio proporcion para train
#' @param seed semilla
#' @return lista con resultados y veredicto
ue_test_out_of_sample <- function(trayectorias, observable, theta_star, tau_star,
                                   above = TRUE, train_ratio = 0.7, seed = 123) {
  set.seed(seed)
  n <- length(trayectorias)
  idx_train <- sample(1:n, size = round(train_ratio * n))
  idx_test <- setdiff(1:n, idx_train)

  # Funcion auxiliar para calcular metricas
  calcular <- function(trajs) {
    filtrado <- ue_filtrar_ERS(trajs, observable, theta_star, tau_star, above)
    if (filtrado$n_condicionadas < 30) return(list(D_H = NA, clase_dom = NA))

    dist <- ue_distribucion_condicionada(filtrado$trayectorias)
    paisaje <- ue_paisaje_costo(dist)
    clases <- ue_detectar_clases(paisaje, filtrado$trayectorias)
    metricas <- ue_metricas_dominancia(clases)

    list(D_H = metricas$D_H, clase_dom = which.max(clases$pesos), pesos = clases$pesos)
  }

  res_train <- calcular(trayectorias[idx_train])
  res_test <- calcular(trayectorias[idx_test])

  mismo_clase <- !is.na(res_train$clase_dom) && !is.na(res_test$clase_dom) &&
                 res_train$clase_dom == res_test$clase_dom

  diff_D_H <- abs(res_train$D_H - res_test$D_H) / res_train$D_H
  D_H_similar <- !is.na(diff_D_H) && diff_D_H < 0.2

  list(
    train = res_train,
    test = res_test,
    diff_D_H = diff_D_H,
    pass = mismo_clase && D_H_similar,
    mensaje = ifelse(mismo_clase && D_H_similar,
                     "Out-of-sample verificado",
                     "Out-of-sample NO verificado")
  )
}

#' Ejecutar validacion completa
#' @param trayectorias todas las trayectorias
#' @param trayectorias_R trayectorias condicionadas
#' @param observable funcion observable
#' @param theta_star umbral
#' @param tau_star persistencia
#' @param above direccion
#' @return lista con los tres tests
ue_validacion_completa <- function(trayectorias, trayectorias_R, observable,
                                    theta_star, tau_star, above = TRUE) {
  cat("Ejecutando validacion UE...\n")

  cat("  [1/3] Test de robustez...\n")
  robustez <- ue_test_robustez(trayectorias, observable, theta_star, tau_star, above)

  cat("  [2/3] Test de ablacion...\n")
  ablacion <- ue_test_ablacion(trayectorias_R, observable, theta_star, tau_star, above)

  cat("  [3/3] Test out-of-sample...\n")
  oos <- ue_test_out_of_sample(trayectorias, observable, theta_star, tau_star, above)

  all_pass <- robustez$pass && ablacion$pass && oos$pass

  cat("\n")
  cat("Resultados:\n")
  cat("  Robustez:      ", ifelse(robustez$pass, "PASS", "FAIL"), "\n")
  cat("  Ablacion:      ", ifelse(ablacion$pass, "PASS", "PARCIAL"), "\n")
  cat("  Out-of-Sample: ", ifelse(oos$pass, "PASS", "FAIL"), "\n")
  cat("  ---\n")
  cat("  VALIDACION:    ", ifelse(all_pass, "COMPLETA", "PARCIAL"), "\n")

  list(
    robustez = robustez,
    ablacion = ablacion,
    out_of_sample = oos,
    all_pass = all_pass
  )
}

# =============================================================================
# 8. PIPELINE COMPLETO
# =============================================================================

#' Ejecutar pipeline UE completo
#' @param grad_V gradiente del potencial
#' @param observable funcion observable
#' @param theta_star umbral ERS
#' @param tau_star persistencia ERS
#' @param above direccion condicion
#' @param n_traj numero de trayectorias
#' @param n_pasos pasos por trayectoria
#' @param dt paso temporal
#' @param sigma ruido
#' @param seed semilla
#' @param validar ejecutar validacion
#' @param verbose imprimir progreso
#' @return lista con todos los resultados
ue_pipeline <- function(grad_V, observable,
                         theta_star, tau_star, above = TRUE,
                         n_traj = 2000, n_pasos = 500,
                         dt = 0.01, sigma = 0.4,
                         seed = 42, validar = FALSE, verbose = TRUE) {

  if (verbose) cat("=== UE Pipeline ===\n\n")

  # 1. Crear dinamica
  if (verbose) cat("[1/6] Creando dinamica...\n")
  step_fn <- ue_crear_dinamica(grad_V, dt, sigma)

  # 2. Simular
  if (verbose) cat("[2/6] Simulando", n_traj, "trayectorias...\n")
  trayectorias <- ue_simular_ensemble(n_traj, n_pasos, step_fn, seed = seed)

  # 3. Filtrar ERS
  if (verbose) cat("[3/6] Filtrando por condicion R...\n")
  filtrado <- ue_filtrar_ERS(trayectorias, observable, theta_star, tau_star, above)
  if (verbose) cat("      Condicionadas:", filtrado$n_condicionadas,
                   "(", round(filtrado$tasa * 100, 1), "%)\n")

  if (filtrado$n_condicionadas < 50) {
    warning("Muy pocas trayectorias condicionadas. Ajustar parametros.")
    return(list(error = "Pocas trayectorias condicionadas", filtrado = filtrado))
  }

  # 4. Distribucion condicionada
  if (verbose) cat("[4/6] Calculando P(x|R)...\n")
  distribucion <- ue_distribucion_condicionada(filtrado$trayectorias)

  # 5. Paisaje de costo
  if (verbose) cat("[5/6] Construyendo paisaje F_R...\n")
  paisaje <- ue_paisaje_costo(distribucion)

  # 6. Clases y dominancia
  if (verbose) cat("[6/6] Detectando clases y metricas...\n")
  clases <- ue_detectar_clases(paisaje, filtrado$trayectorias)
  metricas <- ue_metricas_dominancia(clases)

  if (verbose) {
    cat("\n=== Resultados ===\n")
    cat("Clases detectadas:", clases$n_clases, "\n")
    cat("Centros:", round(clases$centros, 3), "\n")
    cat("Pesos P(A_k|R):", round(clases$pesos, 4), "\n")
    cat("D_H =", round(metricas$D_H, 3), "\n")
    cat("Interpretacion:", metricas$interpretacion, "\n")
  }

  resultado <- list(
    trayectorias = trayectorias,
    filtrado = filtrado,
    distribucion = distribucion,
    paisaje = paisaje,
    clases = clases,
    metricas = metricas,
    parametros = list(
      theta_star = theta_star,
      tau_star = tau_star,
      above = above,
      n_traj = n_traj,
      n_pasos = n_pasos,
      dt = dt,
      sigma = sigma
    )
  )

  # Validacion opcional
  if (validar) {
    if (verbose) cat("\n")
    resultado$validacion <- ue_validacion_completa(
      trayectorias, filtrado$trayectorias,
      observable, theta_star, tau_star, above
    )
  }

  resultado
}

# =============================================================================
# 9. VISUALIZACION
# =============================================================================

#' Graficar resultados del pipeline
#' @param resultado salida de ue_pipeline
#' @param potencial funcion V(x) para comparar (opcional
ue_plot <- function(resultado, potencial = NULL) {
  par(mfrow = c(2, 2))

  # 1. Distribucion P(x|R)
  plot(resultado$distribucion$x, resultado$distribucion$prob,
       type = "h", col = "steelblue", lwd = 2,
       main = "P(x|R)", xlab = "x", ylab = "probabilidad")

  # 2. Paisaje F_R con clases
  plot(resultado$paisaje$x, resultado$paisaje$F_R,
       type = "l", lwd = 2, col = "darkgreen",
       main = "Paisaje F_R(x)", xlab = "x", ylab = "F_R")
  points(resultado$clases$centros, resultado$clases$F_R_minimos,
         col = "red", pch = 19, cex = 1.5)

  # 3. Pesos de clase
  if (resultado$clases$n_clases > 0) {
    barplot(resultado$clases$pesos,
            names.arg = paste0("A", 1:resultado$clases$n_clases),
            col = rainbow(resultado$clases$n_clases, alpha = 0.7),
            main = "Pesos P(A_k|R)", ylab = "probabilidad")
  }

  # 4. Comparacion con potencial o metricas
  if (!is.null(potencial)) {
    x_grid <- resultado$paisaje$x
    V_vals <- sapply(x_grid, potencial)
    V_norm <- (V_vals - min(V_vals)) / (max(V_vals) - min(V_vals))
    F_norm <- resultado$paisaje$F_R / max(resultado$paisaje$F_R)

    plot(x_grid, V_norm, type = "l", lwd = 2, col = "darkred",
         main = "V(x) vs F_R(x)", xlab = "x", ylab = "normalizado")
    lines(x_grid, F_norm, lwd = 2, col = "darkgreen")
    legend("topright", c("V(x)", "F_R(x)"), col = c("darkred", "darkgreen"), lwd = 2)
  } else {
    # Metricas
    plot.new()
    text(0.5, 0.7, paste("D_H =", round(resultado$metricas$D_H, 2)), cex = 1.5)
    text(0.5, 0.5, paste("Delta =", round(resultado$metricas$Delta, 3)), cex = 1.2)
    text(0.5, 0.3, resultado$metricas$interpretacion, cex = 1.1, col = "darkblue")
  }

  par(mfrow = c(1, 1))
}

cat("UE Framework cargado. Funciones disponibles:\n")
cat("  - ue_pipeline()           : ejecutar analisis completo\n")
cat("  - ue_plot()               : visualizar resultados\n")
cat("  - ue_validacion_completa(): validar resultados\n")
cat("\nEjemplo:\n")
cat("  grad_V <- function(x) 4*x*(x^2-1)\n")
cat("  obs <- function(x) abs(x)\n")
cat("  res <- ue_pipeline(grad_V, obs, theta_star=0.5, tau_star=25)\n")
cat("  ue_plot(res)\n")
