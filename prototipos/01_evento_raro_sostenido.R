# =============================================================================
# CONCEPTO 1: Evento Raro Sostenido (ERS / R)
# =============================================================================
#
# Definicion:
#   R = { θ(t) ≥ θ* durante τ ≥ τ* }
#
# Componentes:
#   - θ(t): observable medible en el tiempo
#   - θ*: umbral (threshold)
#   - τ*: tiempo minimo de persistencia
#
# El ERS filtra trayectorias que mantienen una condicion durante suficiente
# tiempo. No basta con cruzar el umbral brevemente.
# =============================================================================

cat("=== PRUEBA: Evento Raro Sostenido (ERS) ===\n\n")

# --- PARAMETROS ---
set.seed(42)
n_traj <- 1000
traj_length <- 200
dt <- 0.1

# Condicion ERS
theta_threshold <- 0.8   # umbral
tau_persistence <- 20    # pasos minimos

# --- GENERAR TRAYECTORIAS (random walk simple) ---
generate_trajectory <- function(n) {
  cumsum(rnorm(n, 0, 0.1))
}

# --- OBSERVABLE: valor absoluto ---
observable <- function(x) abs(x)

# --- FUNCION: Detectar si cumple ERS ---
cumple_ERS <- function(traj, theta_star, tau_star) {
  theta_t <- observable(traj)
  condicion <- theta_t >= theta_star

  # Buscar racha continua >= tau_star
  if (!any(condicion)) return(FALSE)

  rachas <- rle(condicion)
  duraciones_cumplidas <- rachas$lengths[rachas$values == TRUE]

  if (length(duraciones_cumplidas) == 0) return(FALSE)

  max(duraciones_cumplidas) >= tau_star
}

# --- SIMULAR Y CLASIFICAR ---
cat("Simulando", n_traj, "trayectorias...\n")
cat("Condicion R: |x(t)| >=", theta_threshold, "durante >=", tau_persistence, "pasos\n\n")

trayectorias <- vector("list", n_traj)
cumple_R <- logical(n_traj)

for (i in 1:n_traj) {
  trayectorias[[i]] <- generate_trajectory(traj_length)
  cumple_R[i] <- cumple_ERS(trayectorias[[i]], theta_threshold, tau_persistence)
}

# --- RESULTADOS ---
n_condicionadas <- sum(cumple_R)
porcentaje <- 100 * n_condicionadas / n_traj

cat("Resultados:\n")
cat("  Trayectorias totales:", n_traj, "\n")
cat("  Cumplen R:", n_condicionadas, "(", round(porcentaje, 2), "%)\n")
cat("  No cumplen R:", n_traj - n_condicionadas, "\n\n")

# --- VISUALIZACION ---
par(mfrow = c(2, 2))

# 1. Ejemplo de trayectoria que NO cumple R
idx_no <- which(!cumple_R)[1]
plot(trayectorias[[idx_no]], type = "l", col = "gray50",
     main = "Trayectoria que NO cumple R",
     xlab = "tiempo", ylab = "x(t)")
abline(h = c(-theta_threshold, theta_threshold), col = "red", lty = 2)
legend("topright", "umbral θ*", col = "red", lty = 2, cex = 0.8)

# 2. Ejemplo de trayectoria que SI cumple R
idx_si <- which(cumple_R)[1]
plot(trayectorias[[idx_si]], type = "l", col = "steelblue", lwd = 1.5,
     main = "Trayectoria que SI cumple R",
     xlab = "tiempo", ylab = "x(t)")
abline(h = c(-theta_threshold, theta_threshold), col = "red", lty = 2)

# 3. Observable theta(t) de la trayectoria que cumple
theta_ejemplo <- observable(trayectorias[[idx_si]])
plot(theta_ejemplo, type = "l", col = "darkgreen",
     main = "Observable θ(t) = |x(t)|",
     xlab = "tiempo", ylab = "θ(t)")
abline(h = theta_threshold, col = "red", lty = 2)

# Marcar donde cumple la condicion
cumple_en_t <- theta_ejemplo >= theta_threshold
points(which(cumple_en_t), theta_ejemplo[cumple_en_t],
       col = "orange", pch = 16, cex = 0.5)

# 4. Distribucion de duracion maxima sobre umbral
duracion_max <- sapply(trayectorias, function(traj) {
  theta_t <- observable(traj)
  condicion <- theta_t >= theta_threshold
  if (!any(condicion)) return(0)
  rachas <- rle(condicion)
  duraciones <- rachas$lengths[rachas$values == TRUE]
  if (length(duraciones) == 0) return(0)
  max(duraciones)
})

hist(duracion_max, breaks = 30, col = "lightblue", border = "white",
     main = "Distribucion de duracion maxima sobre θ*",
     xlab = "duracion maxima (pasos)")
abline(v = tau_persistence, col = "red", lwd = 2)
text(tau_persistence + 5, par("usr")[4] * 0.9, "τ*", col = "red")

par(mfrow = c(1, 1))

# --- CONCLUSION ---
cat("\n--- CONCLUSION ---\n")
cat("El ERS actua como un FILTRO:\n")
cat("  - Descarta fluctuaciones breves\n")
cat("  - Selecciona solo trayectorias con comportamiento SOSTENIDO\n")
cat("  - La rareza viene de exigir persistencia, no solo cruzar el umbral\n")
