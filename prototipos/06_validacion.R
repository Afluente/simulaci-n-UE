# =============================================================================
# CONCEPTO 6: Protocolo de Validacion UE
# =============================================================================
#
# Tres tests obligatorios para validar resultados:
#
# 1. ROBUSTEZ: Variar parametros ±20%
#    - θ* → θ* ± δθ
#    - τ* → τ* ± δτ
#    - Verificar: clases dominantes persisten
#
# 2. ABLACION: Romper mecanismos preservando marginales
#    - Permutar temporalmente señales
#    - Aleatorizar conexiones
#    - Verificar: dominancia desaparece o cambia
#
# 3. OUT-OF-SAMPLE: Dividir datos en train/test
#    - Entrenar en 70%, predecir en 30%
#    - Verificar: resultados consistentes
#
# =============================================================================

cat("=== PRUEBA: Protocolo de Validacion UE ===\n\n")

# --- FUNCIONES AUXILIARES ---
set.seed(42)

# Potencial doble pozo asimetrico
V <- function(x) (x^2 - 1)^2 - 0.3 * x
grad_V <- function(x) 4 * x * (x^2 - 1) - 0.3

step <- function(x, dt, sigma) {
  x - grad_V(x) * dt + sigma * sqrt(dt) * rnorm(1)
}

simulate <- function(x0, n, dt, sigma) {
  x <- numeric(n)
  x[1] <- x0
  for (t in 2:n) {
    x[t] <- step(x[t-1], dt, sigma)
    x[t] <- max(-3, min(3, x[t]))
  }
  x
}

cumple_ERS <- function(traj, theta_star, tau_star, observable_fn) {
  theta_t <- observable_fn(traj)
  condicion <- theta_t >= theta_star
  if (!any(condicion)) return(FALSE)
  rachas <- rle(condicion)
  duraciones <- rachas$lengths[rachas$values == TRUE]
  if (length(duraciones) == 0) return(FALSE)
  max(duraciones) >= tau_star
}

# Calcular metricas de una simulacion
run_simulation <- function(n_traj, traj_length, dt, noise,
                            theta_star, tau_star, observable_fn, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  trayectorias <- vector("list", n_traj)
  cumple_R <- logical(n_traj)

  for (i in 1:n_traj) {
    x0 <- runif(1, -2, 2)
    trayectorias[[i]] <- simulate(x0, traj_length, dt, noise)
    cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_star, tau_star, observable_fn)
  }

  x_cond <- unlist(trayectorias[cumple_R])

  if (length(x_cond) < 50) {
    return(list(n_cond = sum(cumple_R), D_H = NA, p_izq = NA, p_der = NA))
  }

  p_izq <- mean(x_cond < 0)
  p_der <- mean(x_cond > 0)
  p_max <- max(p_izq, p_der)
  p_min <- min(p_izq, p_der)
  D_H <- ifelse(p_min > 0, p_max / p_min, Inf)

  list(n_cond = sum(cumple_R), D_H = D_H, p_izq = p_izq, p_der = p_der,
       trayectorias = trayectorias, cumple_R = cumple_R)
}

# --- PARAMETROS BASE ---
n_traj <- 1500
traj_length <- 400
dt <- 0.01
noise <- 0.4
theta_base <- 0.5
tau_base <- 25
observable <- function(x) abs(x)

cat("Parametros base:\n")
cat("  θ* =", theta_base, "\n")
cat("  τ* =", tau_base, "\n\n")

# --- SIMULACION BASE ---
cat("=== SIMULACION BASE ===\n")
base <- run_simulation(n_traj, traj_length, dt, noise,
                        theta_base, tau_base, observable, seed = 42)
cat("Condicionadas:", base$n_cond, "\n")
cat("D_H =", round(base$D_H, 3), "\n")
cat("P(izq|R) =", round(base$p_izq, 3), ", P(der|R) =", round(base$p_der, 3), "\n")
clase_dominante_base <- ifelse(base$p_der > base$p_izq, "derecha", "izquierda")
cat("Clase dominante:", clase_dominante_base, "\n\n")

# ============================================================================
# TEST 1: ROBUSTEZ
# ============================================================================
cat("=== TEST 1: ROBUSTEZ (±20% en parametros) ===\n\n")

variaciones <- expand.grid(
  theta_mult = c(0.8, 1.0, 1.2),
  tau_mult = c(0.8, 1.0, 1.2)
)

robustez_resultados <- data.frame()

for (i in 1:nrow(variaciones)) {
  theta_var <- theta_base * variaciones$theta_mult[i]
  tau_var <- round(tau_base * variaciones$tau_mult[i])

  res <- run_simulation(n_traj, traj_length, dt, noise,
                         theta_var, tau_var, observable, seed = 42 + i)

  clase_dom <- ifelse(is.na(res$p_der), NA,
                      ifelse(res$p_der > res$p_izq, "derecha", "izquierda"))

  robustez_resultados <- rbind(robustez_resultados, data.frame(
    theta = theta_var,
    tau = tau_var,
    n_cond = res$n_cond,
    D_H = round(res$D_H, 2),
    clase_dom = clase_dom
  ))
}

print(robustez_resultados)

# Verificar consistencia
clases_consistentes <- all(robustez_resultados$clase_dom == clase_dominante_base, na.rm = TRUE)
D_H_estable <- sd(robustez_resultados$D_H, na.rm = TRUE) / mean(robustez_resultados$D_H, na.rm = TRUE)

cat("\nResultado:\n")
if (clases_consistentes && D_H_estable < 0.3) {
  cat("  ✓ ROBUSTEZ VERIFICADA\n")
  cat("    - Clase dominante consistente:", clase_dominante_base, "\n")
  cat("    - Variabilidad D_H:", round(D_H_estable * 100, 1), "%\n")
} else {
  cat("  ✗ ROBUSTEZ NO VERIFICADA\n")
  cat("    - Clase dominante cambia o D_H muy variable\n")
}

# ============================================================================
# TEST 2: ABLACION
# ============================================================================
cat("\n=== TEST 2: ABLACION (romper estructura temporal) ===\n\n")

# Obtener trayectorias condicionadas originales
traj_cond_orig <- base$trayectorias[base$cumple_R]

# ABLACION: Permutar temporalmente cada trayectoria
# Esto preserva la distribucion marginal P(x) pero rompe correlaciones temporales
cat("Metodo: Permutar aleatoriamente cada trayectoria en el tiempo\n")
cat("(preserva P(x), rompe estructura temporal)\n\n")

traj_ablacion <- lapply(traj_cond_orig, function(traj) sample(traj))

# Verificar cuantas "cumplirian R" despues de ablacion
# (con la estructura temporal destruida)
n_cumple_ablacion <- sum(sapply(traj_ablacion, function(traj) {
  cumple_ERS(traj, theta_base, tau_base, observable)
}))

cat("Trayectorias originales que cumplen R:", length(traj_cond_orig), "\n")
cat("Trayectorias ablacionadas que 'cumplen R':", n_cumple_ablacion, "\n")

# Calcular D_H de ablacionadas (si hay suficientes)
if (n_cumple_ablacion > 20) {
  x_ablacion <- unlist(traj_ablacion[sapply(traj_ablacion, function(traj) {
    cumple_ERS(traj, theta_base, tau_base, observable)
  })])
  p_izq_abl <- mean(x_ablacion < 0)
  p_der_abl <- mean(x_ablacion > 0)
  D_H_abl <- max(p_izq_abl, p_der_abl) / min(p_izq_abl, p_der_abl)
  cat("D_H ablacion:", round(D_H_abl, 2), "\n")
} else {
  D_H_abl <- NA
  cat("Muy pocas trayectorias ablacionadas cumplen R\n")
}

cat("\nResultado:\n")
if (n_cumple_ablacion < length(traj_cond_orig) * 0.3) {
  cat("  ✓ ABLACION VERIFICADA\n")
  cat("    - La estructura temporal es crucial para R\n")
  cat("    - Sin ella, las trayectorias no cumplen la condicion\n")
} else {
  cat("  ? ABLACION PARCIAL\n")
  cat("    - Algunas trayectorias aun cumplen R tras ablacion\n")
  if (!is.na(D_H_abl) && abs(D_H_abl - base$D_H) > 0.5) {
    cat("    - Pero D_H cambia significativamente\n")
  }
}

# ============================================================================
# TEST 3: OUT-OF-SAMPLE
# ============================================================================
cat("\n=== TEST 3: OUT-OF-SAMPLE (train 70% / test 30%) ===\n\n")

# Dividir trayectorias en train y test
set.seed(123)
n_total <- length(base$trayectorias)
idx_train <- sample(1:n_total, size = round(0.7 * n_total))
idx_test <- setdiff(1:n_total, idx_train)

traj_train <- base$trayectorias[idx_train]
traj_test <- base$trayectorias[idx_test]

cumple_train <- base$cumple_R[idx_train]
cumple_test <- base$cumple_R[idx_test]

# Metricas en train
x_train <- unlist(traj_train[cumple_train])
p_izq_train <- mean(x_train < 0)
p_der_train <- mean(x_train > 0)
D_H_train <- max(p_izq_train, p_der_train) / min(p_izq_train, p_der_train)

# Metricas en test
x_test <- unlist(traj_test[cumple_test])
p_izq_test <- mean(x_test < 0)
p_der_test <- mean(x_test > 0)
D_H_test <- max(p_izq_test, p_der_test) / min(p_izq_test, p_der_test)

cat("TRAIN (70%):\n")
cat("  N condicionadas:", sum(cumple_train), "\n")
cat("  D_H =", round(D_H_train, 3), "\n")
cat("  P(izq|R) =", round(p_izq_train, 3), ", P(der|R) =", round(p_der_train, 3), "\n\n")

cat("TEST (30%):\n")
cat("  N condicionadas:", sum(cumple_test), "\n")
cat("  D_H =", round(D_H_test, 3), "\n")
cat("  P(izq|R) =", round(p_izq_test, 3), ", P(der|R) =", round(p_der_test, 3), "\n\n")

# Comparar
diff_D_H <- abs(D_H_train - D_H_test) / D_H_train
clase_train <- ifelse(p_der_train > p_izq_train, "derecha", "izquierda")
clase_test <- ifelse(p_der_test > p_izq_test, "derecha", "izquierda")

cat("Resultado:\n")
if (clase_train == clase_test && diff_D_H < 0.2) {
  cat("  ✓ OUT-OF-SAMPLE VERIFICADO\n")
  cat("    - Clase dominante consistente:", clase_train, "\n")
  cat("    - Diferencia D_H:", round(diff_D_H * 100, 1), "%\n")
} else {
  cat("  ✗ OUT-OF-SAMPLE NO VERIFICADO\n")
  cat("    - Resultados difieren entre train y test\n")
}

# ============================================================================
# RESUMEN FINAL
# ============================================================================
cat("\n")
cat("=============================================\n")
cat("         RESUMEN DE VALIDACION UE           \n")
cat("=============================================\n\n")

test_robustez <- clases_consistentes && D_H_estable < 0.3
test_ablacion <- n_cumple_ablacion < length(traj_cond_orig) * 0.3
test_oos <- clase_train == clase_test && diff_D_H < 0.2

cat("1. Robustez:      ", ifelse(test_robustez, "✓ PASS", "✗ FAIL"), "\n")
cat("2. Ablacion:      ", ifelse(test_ablacion, "✓ PASS", "? PARCIAL"), "\n")
cat("3. Out-of-Sample: ", ifelse(test_oos, "✓ PASS", "✗ FAIL"), "\n\n")

if (test_robustez && test_ablacion && test_oos) {
  cat("==> VALIDACION COMPLETA: Resultados confiables\n")
} else {
  cat("==> VALIDACION PARCIAL: Revisar tests fallidos\n")
}

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Robustez: D_H vs parametros
plot(1:nrow(robustez_resultados), robustez_resultados$D_H,
     type = "b", pch = 19, col = "steelblue",
     main = "Test Robustez: D_H estable?",
     xlab = "configuracion", ylab = "D_H",
     ylim = c(0, max(robustez_resultados$D_H, na.rm = TRUE) * 1.2))
abline(h = base$D_H, col = "red", lty = 2)
legend("topright", "D_H base", col = "red", lty = 2)

# 2. Ablacion: comparacion
barplot(c(length(traj_cond_orig), n_cumple_ablacion),
        names.arg = c("Original", "Ablacion"),
        col = c("steelblue", "coral"),
        main = "Test Ablacion: R requiere estructura?",
        ylab = "N trayectorias que cumplen R")

# 3. Out-of-sample: comparacion
barplot(matrix(c(p_izq_train, p_der_train, p_izq_test, p_der_test), nrow = 2),
        beside = TRUE, col = c("steelblue", "coral"),
        names.arg = c("Train", "Test"),
        main = "Test OOS: Consistencia?",
        ylab = "P(clase|R)")
legend("topright", c("Izq", "Der"), fill = c("steelblue", "coral"))

# 4. Resumen visual
plot(1, type = "n", xlim = c(0, 10), ylim = c(0, 4),
     axes = FALSE, xlab = "", ylab = "",
     main = "Resumen Validacion")
text(5, 3, paste("Robustez:", ifelse(test_robustez, "PASS", "FAIL")),
     col = ifelse(test_robustez, "darkgreen", "red"), cex = 1.3)
text(5, 2, paste("Ablacion:", ifelse(test_ablacion, "PASS", "PARCIAL")),
     col = ifelse(test_ablacion, "darkgreen", "orange"), cex = 1.3)
text(5, 1, paste("Out-of-Sample:", ifelse(test_oos, "PASS", "FAIL")),
     col = ifelse(test_oos, "darkgreen", "red"), cex = 1.3)

par(mfrow = c(1, 1))

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("El protocolo de validacion UE asegura que:\n")
cat("  1. Los resultados no dependen de ajuste fino de parametros (Robustez)\n")
cat("  2. La estructura del sistema es necesaria, no accidental (Ablacion)\n")
cat("  3. Los resultados generalizan a datos no vistos (Out-of-Sample)\n")
cat("\nSin estos tests, hay riesgo de 'posdiction' o sobreajuste.\n")
